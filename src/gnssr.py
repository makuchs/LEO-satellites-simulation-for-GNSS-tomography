"""
GNSS-R (GNSS Reflectometry) Geometry Simulation Module

This module provides a SPICE-based workflow to simulate bistatic GNSS-R geometry between
GNSS transmitters (Tx) and LEO receivers (Rx). It computes an approximate specular point
on the Earth ellipsoid (from SPICE body radii), filters events by a geographic bounding
box, checks visibility from Tx and Rx to the specular point, and writes results to CSV.
"""

import os
import sys
import csv
import time
import argparse
import subprocess
import multiprocessing as mp
import tempfile

import numpy as np
import pandas as pd
import spiceypy as spice
from spiceypy.utils.support_types import SPICEINT_CELL, SPICEDOUBLE_CELL
from project_paths import DE432_FILE, EARTH_BPC_FILE, LSK_FILE, PCK_TPC_FILE

# Populated by _earth_radii_and_flattening(); cleared whenever spice.kclear() runs in this module.
_EARTH_RADII_FLAT_CACHE = None


def _invalidate_earth_geometry_cache():
    global _EARTH_RADII_FLAT_CACHE
    _EARTH_RADII_FLAT_CACHE = None


def load_spice_kernels(leo_path, gnss_path):
    spice.kclear()
    _invalidate_earth_geometry_cache()
    spice.furnsh(str(LSK_FILE))
    spice.furnsh(str(PCK_TPC_FILE))
    spice.furnsh(str(DE432_FILE))
    if EARTH_BPC_FILE.exists():
        spice.furnsh(str(EARTH_BPC_FILE))
    spice.furnsh(leo_path)
    spice.furnsh(gnss_path)


def load_time_kernel():
    spice.kclear()
    _invalidate_earth_geometry_cache()
    spice.furnsh(str(LSK_FILE))


def get_files_from_folder(folder_path):
    return [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if os.path.isfile(os.path.join(folder_path, f)) and f.lower().endswith(".bsp")
    ]


def get_objects_from_spk_file(spk_path):
    ids_cell = spice.spkobj(spk_path, SPICEINT_CELL(10000))
    obj_ids = sorted(int(x) for x in ids_cell)

    out = []
    for obj_id in obj_ids:
        try:
            obj_name = spice.bodc2s(int(obj_id))
        except Exception:
            obj_name = str(obj_id)

        out.append({"spice_id": int(obj_id), "name": str(obj_name)})
    return out


def get_spk_coverage_windows(spk_path, obj_id):
    try:
        cover = spice.spkcov(spk_path, int(obj_id), SPICEDOUBLE_CELL(200000))
        nwin = spice.wncard(cover)
        return [spice.wnfetd(cover, i) for i in range(nwin)]
    except Exception:
        return []


def filter_epochs_by_windows(epochs, windows):
    if len(epochs) == 0 or not windows:
        return np.array([], dtype=float)

    mask = np.zeros(len(epochs), dtype=bool)
    for et_start, et_stop in windows:
        mask |= (epochs >= float(et_start)) & (epochs <= float(et_stop))
    return epochs[mask]


def _earth_radii_and_flattening():
    global _EARTH_RADII_FLAT_CACHE
    if _EARTH_RADII_FLAT_CACHE is not None:
        return _EARTH_RADII_FLAT_CACHE
    radii = spice.bodvrd("EARTH", "RADII", 3)[1]
    re = float(radii[0])
    rp = float(radii[2])
    f = (re - rp) / re
    _EARTH_RADII_FLAT_CACHE = (re, rp, f)
    return _EARTH_RADII_FLAT_CACHE


def _state_in_frame(target_id, et, out_frame="J2000", abcorr="NONE", obs="EARTH"):
    st, _ = spice.spkezr(str(target_id), float(et), out_frame, abcorr, obs)
    pos = np.array(st[:3], dtype=float)
    vel = np.array(st[3:6], dtype=float)
    return pos, vel


def _transform_state(pos, vel, et, from_frame="J2000", to_frame="IAU_EARTH"):
    xform = spice.sxform(from_frame, to_frame, float(et))
    st6 = np.hstack([pos, vel])
    st6_out = spice.mxvg(xform, st6)
    return np.array(st6_out[:3], dtype=float), np.array(st6_out[3:6], dtype=float)


def _ecef_latlon_alt_from_pos(pos_ecef):
    re, _, f = _earth_radii_and_flattening()
    lon, lat, alt = spice.recgeo(pos_ecef, re, f)
    lat_deg = float(spice.convrt(lat, "RADIANS", "DEGREES"))
    lon_deg = float(spice.convrt(lon, "RADIANS", "DEGREES"))
    return lat_deg, lon_deg, float(alt)


def _wrap_lon_deg(lon_deg):
    return (float(lon_deg) + 180.0) % 360.0 - 180.0


def is_in_bbox(lat_deg, lon_deg, bbox):
    min_lat, max_lat, min_lon, max_lon = bbox

    lat = float(lat_deg)
    lon = _wrap_lon_deg(float(lon_deg))
    min_lon = _wrap_lon_deg(float(min_lon))
    max_lon = _wrap_lon_deg(float(max_lon))

    if not (float(min_lat) <= lat <= float(max_lat)):
        return False

    if min_lon <= max_lon:
        return min_lon <= lon <= max_lon

    return (lon >= min_lon) or (lon <= max_lon)


def _surface_point_from_latlon(lat_deg, lon_deg):
    re, _, f = _earth_radii_and_flattening()
    lat = spice.convrt(float(lat_deg), "DEGREES", "RADIANS")
    lon = spice.convrt(float(lon_deg), "DEGREES", "RADIANS")
    return np.array(spice.pgrrec("EARTH", lon, lat, 0.0, re, f), dtype=float)


def _surface_normal_ecef(surfpt_ecef):
    re, rp, _ = _earth_radii_and_flattening()
    return np.array(spice.surfnm(re, re, rp, surfpt_ecef), dtype=float)


def _specular_error(tx_ecef, rx_ecef, s_ecef):
    n = _surface_normal_ecef(s_ecef)

    u = tx_ecef - s_ecef
    v = rx_ecef - s_ecef
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu < 1e-12 or nv < 1e-12:
        return 1e99

    uhat = u / nu
    vhat = v / nv
    w = uhat + vhat
    nw = np.linalg.norm(w)
    if nw < 1e-12:
        return 1e99

    what = w / nw
    e1 = np.linalg.norm(what - n)
    e2 = np.linalg.norm(what + n)
    return float(min(e1, e2))


def find_specular_point(tx_ecef, rx_ecef, initial_latlon=None,
                        coarse_deg=5.0, refine_iters=7, refine_factor=0.5):
    if initial_latlon is None:
        lat_tx, lon_tx, _ = _ecef_latlon_alt_from_pos(tx_ecef)
        lat_rx, lon_rx, _ = _ecef_latlon_alt_from_pos(rx_ecef)

        lon_tx = _wrap_lon_deg(lon_tx)
        lon_rx = _wrap_lon_deg(lon_rx)
        dlon = _wrap_lon_deg(lon_rx - lon_tx)

        lat0 = 0.5 * (lat_tx + lat_rx)
        lon0 = _wrap_lon_deg(lon_tx + 0.5 * dlon)
    else:
        lat0, lon0 = initial_latlon
        lon0 = _wrap_lon_deg(lon0)

    best_lat = float(lat0)
    best_lon = float(lon0)
    best_s = _surface_point_from_latlon(best_lat, best_lon)
    best_err = _specular_error(tx_ecef, rx_ecef, best_s)

    step = float(coarse_deg)

    for _ in range(int(refine_iters)):
        cand_lat = np.array([best_lat - step, best_lat, best_lat + step], dtype=float)
        cand_lat = np.clip(cand_lat, -89.9, 89.9)

        cand_lon = np.array([best_lon - step, best_lon, best_lon + step], dtype=float)
        cand_lon = np.array([_wrap_lon_deg(x) for x in cand_lon], dtype=float)

        for la in cand_lat:
            for lo in cand_lon:
                s = _surface_point_from_latlon(float(la), float(lo))
                err = _specular_error(tx_ecef, rx_ecef, s)
                if err < best_err:
                    best_err = err
                    best_lat = float(la)
                    best_lon = float(lo)
                    best_s = s

        step *= float(refine_factor)

    return best_s, best_lat, best_lon, best_err


def find_specular_point_in_bbox(tx_ecef, rx_ecef, bbox,
                                grid_deg=1.0,
                                refine_iters=7, refine_factor=0.5):
    min_lat, max_lat, min_lon, max_lon = bbox
    min_lat = float(min_lat)
    max_lat = float(max_lat)
    min_lon = float(min_lon)
    max_lon = float(max_lon)

    if min_lat > max_lat:
        min_lat, max_lat = max_lat, min_lat

    if grid_deg <= 0:
        grid_deg = 1.0

    lats = np.arange(min_lat, max_lat + 1e-9, grid_deg, dtype=float)

    if _wrap_lon_deg(min_lon) <= _wrap_lon_deg(max_lon) and abs(max_lon - min_lon) < 180:
        lons = np.arange(min_lon, max_lon + 1e-9, grid_deg, dtype=float)
    else:
        lon_a = np.arange(min_lon, 180.0 + 1e-9, grid_deg, dtype=float)
        lon_b = np.arange(-180.0, max_lon + 1e-9, grid_deg, dtype=float)
        lons = np.concatenate([lon_a, lon_b]).astype(float)

    best_err = 1e99
    best_lat = float(lats[len(lats) // 2]) if len(lats) else 0.0
    best_lon = float(_wrap_lon_deg(lons[len(lons) // 2])) if len(lons) else 0.0
    best_s = _surface_point_from_latlon(best_lat, best_lon)

    for la in lats:
        for lo in lons:
            lo2 = _wrap_lon_deg(float(lo))
            s = _surface_point_from_latlon(float(la), float(lo2))
            err = _specular_error(tx_ecef, rx_ecef, s)
            if err < best_err:
                best_err = err
                best_lat = float(la)
                best_lon = float(lo2)
                best_s = s

    step = float(grid_deg)

    for _ in range(int(refine_iters)):
        cand_lat = np.array([best_lat - step, best_lat, best_lat + step], dtype=float)
        cand_lat = np.clip(cand_lat, -89.9, 89.9)

        cand_lon = np.array([best_lon - step, best_lon, best_lon + step], dtype=float)
        cand_lon = np.array([_wrap_lon_deg(x) for x in cand_lon], dtype=float)

        for la in cand_lat:
            for lo in cand_lon:
                if not is_in_bbox(float(la), float(lo), bbox):
                    continue
                s = _surface_point_from_latlon(float(la), float(lo))
                err = _specular_error(tx_ecef, rx_ecef, s)
                if err < best_err:
                    best_err = err
                    best_lat = float(la)
                    best_lon = float(lo)
                    best_s = s

        step *= float(refine_factor)

    return best_s, best_lat, best_lon, best_err


def _has_line_of_sight_to_surface(observer_ecef, target_surface_ecef, tol_km=0.01):
    re, rp, _ = _earth_radii_and_flattening()

    dvec = target_surface_ecef - observer_ecef
    nd = np.linalg.norm(dvec)
    if nd < 1e-12:
        return False

    d = dvec / nd

    try:
        xpt = spice.surfpt(observer_ecef, d, re, re, rp)
        if xpt is None:
            return False
        return np.linalg.norm(np.array(xpt, dtype=float) - target_surface_ecef) <= float(tol_km)
    except Exception:
        return False


def incidence_angle_deg(tx_ecef, s_ecef):
    n = _surface_normal_ecef(s_ecef)
    v = tx_ecef - s_ecef
    nv = np.linalg.norm(v)
    if nv < 1e-12:
        return float("nan")
    vhat = v / nv
    c = float(np.clip(np.dot(n, vhat), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def bistatic_angle_deg(tx_ecef, rx_ecef, s_ecef):
    u = tx_ecef - s_ecef
    v = rx_ecef - s_ecef
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)
    if nu < 1e-12 or nv < 1e-12:
        return float("nan")
    uhat = u / nu
    vhat = v / nv
    c = float(np.clip(np.dot(uhat, vhat), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def compute_delay_km(tx_ecef, rx_ecef, s_ecef):
    d1 = float(np.linalg.norm(tx_ecef - s_ecef))
    d2 = float(np.linalg.norm(rx_ecef - s_ecef))
    ddir = float(np.linalg.norm(tx_ecef - rx_ecef))
    return (d1 + d2), (d1 + d2 - ddir)


def process_gnssr_epoch(
    leo_id,
    gnss_id,
    et,
    bbox=None,
    specular_err_max=0.15,
    inc_min_deg=None,
    inc_max_deg=75.0,
    bistatic_max_deg=60.0,
    excess_max_km=200.0,
    los_tol_km=0.01,
    bbox_guided_search=False,
    bbox_grid_deg=1.0,
    initial_latlon=None
):
    tx_pos_j2k, tx_vel_j2k = _state_in_frame(gnss_id, et, "J2000", "NONE", "EARTH")
    rx_pos_j2k, rx_vel_j2k = _state_in_frame(leo_id, et, "J2000", "NONE", "EARTH")

    tx_ecef, _ = _transform_state(tx_pos_j2k, tx_vel_j2k, et, "J2000", "IAU_EARTH")
    rx_ecef, _ = _transform_state(rx_pos_j2k, rx_vel_j2k, et, "J2000", "IAU_EARTH")

    if bbox_guided_search and (bbox is not None):
        grid_deg = 0.5 if (bbox_grid_deg is None) else float(bbox_grid_deg)
        s_ecef, _, _, s_err = find_specular_point_in_bbox(
            tx_ecef, rx_ecef, bbox=bbox, grid_deg=grid_deg
        )
    else:
        s_ecef, _, _, s_err = find_specular_point(
            tx_ecef, rx_ecef, initial_latlon=initial_latlon
        )

    if not np.isfinite(s_err) or float(s_err) > float(specular_err_max):
        return None

    s_lat, s_lon, _ = _ecef_latlon_alt_from_pos(s_ecef)
    s_lon = _wrap_lon_deg(s_lon)

    if bbox is not None and (not is_in_bbox(float(s_lat), float(s_lon), bbox)):
        return None

    if not _has_line_of_sight_to_surface(tx_ecef, s_ecef, tol_km=los_tol_km):
        return None
    if not _has_line_of_sight_to_surface(rx_ecef, s_ecef, tol_km=los_tol_km):
        return None

    inc = incidence_angle_deg(tx_ecef, s_ecef)
    bist = bistatic_angle_deg(tx_ecef, rx_ecef, s_ecef)
    path_km, excess_km = compute_delay_km(tx_ecef, rx_ecef, s_ecef)
    grazing_elev_deg = 90.0 - float(inc) if np.isfinite(inc) else float("nan")

    if inc_min_deg is not None and np.isfinite(inc) and float(inc) < float(inc_min_deg):
        return None
    if inc_max_deg is not None and np.isfinite(inc) and float(inc) > float(inc_max_deg):
        return None
    if bistatic_max_deg is not None and np.isfinite(bist) and float(bist) > float(bistatic_max_deg):
        return None
    if excess_max_km is not None and np.isfinite(excess_km) and abs(float(excess_km)) > float(excess_max_km):
        return None

    dt_utc = spice.timout(float(et), "YYYY-MM-DD HR:MN:SC ::UTC")

    return {
        "time_utc": dt_utc,
        "GNSS": abs(int(gnss_id)),
        "LEO": abs(int(leo_id)),
        "spec_lat_deg": float(s_lat),
        "spec_lon_deg": float(s_lon),
        "specular_err": float(s_err),
        "incidence_deg": float(inc),
        "grazing_elev_deg": float(grazing_elev_deg),
        "bistatic_deg": float(bist),
        "path_km": float(path_km),
        "excess_m": float(excess_km) * 1000.0,
    }


def _format_output_df(df):
    round_map = {
        "spec_lat_deg": 5,
        "spec_lon_deg": 5,
        "specular_err": 5,
        "incidence_deg": 3,
        "grazing_elev_deg": 3,
        "bistatic_deg": 3,
        "path_km": 3,
        "excess_m": 1,
    }
    for col, nd in round_map.items():
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").round(nd)
    if "specular_err" in df.columns:
        df.loc[df["specular_err"].abs() < 1e-9, "specular_err"] = 0.0
    return df


# =========================
# Batch / subprocess processing
# =========================

def _init_base_kernels_for_worker():
    spice.kclear()
    _invalidate_earth_geometry_cache()
    spice.furnsh(str(LSK_FILE))
    spice.furnsh(str(PCK_TPC_FILE))
    spice.furnsh(str(DE432_FILE))
    if EARTH_BPC_FILE.exists():
        spice.furnsh(str(EARTH_BPC_FILE))


def run_gnssr_single_leo(
    leo_path,
    gnss_path,
    out_csv_path,
    start_date,
    end_date,
    start_time_utc="00:00:00",
    end_time_utc="23:00:00",
    bbox=None,
    step_seconds=3600,
    specular_err_max=1.0,
    inc_min_deg=None,
    inc_max_deg=None,
    bistatic_max_deg=None,
    excess_max_km=None,
    los_tol_km=0.01,
    bbox_guided_search=True,
    bbox_grid_deg=0.25,
    initial_latlon=None
):
    header = [
        "time_utc",
        "GNSS",
        "LEO",
        "spec_lat_deg",
        "spec_lon_deg",
        "specular_err",
        "incidence_deg",
        "grazing_elev_deg",
        "bistatic_deg",
        "path_km",
        "excess_m",
    ]

    _init_base_kernels_for_worker()
    spice.furnsh(gnss_path)
    spice.furnsh(leo_path)

    gnss_objs = get_objects_from_spk_file(gnss_path)
    leo_objs = get_objects_from_spk_file(leo_path)
    dates = pd.date_range(start=start_date, end=end_date, freq="D")

    gnss_cov = {obj["spice_id"]: get_spk_coverage_windows(gnss_path, obj["spice_id"]) for obj in gnss_objs}
    leo_cov = {obj["spice_id"]: get_spk_coverage_windows(leo_path, obj["spice_id"]) for obj in leo_objs}

    rows_saved = 0
    os.makedirs(os.path.dirname(out_csv_path) or ".", exist_ok=True)

    with open(out_csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)

        for date in dates:
            d = date.strftime("%Y %m %d")
            et1 = spice.str2et(f"{d} {start_time_utc} UTC")
            et2 = spice.str2et(f"{d} {end_time_utc} UTC")
            epochs = np.arange(float(et1), float(et2) + 1e-9, float(step_seconds), dtype=float)

            for gnss_obj in gnss_objs:
                gnss_epochs = filter_epochs_by_windows(epochs, gnss_cov.get(gnss_obj["spice_id"], []))
                if len(gnss_epochs) == 0:
                    continue

                for leo_obj in leo_objs:
                    pair_epochs = filter_epochs_by_windows(gnss_epochs, leo_cov.get(leo_obj["spice_id"], []))
                    if len(pair_epochs) == 0:
                        continue

                    rows = []

                    for et in pair_epochs:
                        try:
                            row = process_gnssr_epoch(
                                leo_id=leo_obj["spice_id"],
                                gnss_id=gnss_obj["spice_id"],
                                et=float(et),
                                bbox=bbox,
                                specular_err_max=specular_err_max,
                                inc_min_deg=inc_min_deg,
                                inc_max_deg=inc_max_deg,
                                bistatic_max_deg=bistatic_max_deg,
                                excess_max_km=excess_max_km,
                                los_tol_km=los_tol_km,
                                bbox_guided_search=bbox_guided_search,
                                bbox_grid_deg=bbox_grid_deg,
                                initial_latlon=initial_latlon
                            )
                            if row is not None:
                                rows.append(row)
                        except Exception:
                            continue

                    if rows:
                        df = pd.DataFrame(rows, columns=header)
                        df = _format_output_df(df)
                        for r in df.itertuples(index=False):
                            w.writerow(list(r))
                        rows_saved += len(rows)

    try:
        spice.kclear()
        _invalidate_earth_geometry_cache()
    except Exception:
        pass

    return rows_saved


def run_gnssr(
    leo_folder,
    gnss_folder,
    csv_file_path,
    start_date,
    end_date,
    start_time_utc="00:00:00",
    end_time_utc="23:00:00",
    bbox=None,
    step_seconds=3600,
    specular_err_max=1.0,
    inc_min_deg=None,
    inc_max_deg=None,
    bistatic_max_deg=None,
    excess_max_km=None,
    los_tol_km=0.01,
    bbox_guided_search=True,
    bbox_grid_deg=0.25,
    initial_latlon=None,
):
    leo_files = sorted(get_files_from_folder(leo_folder))
    gnss_files = sorted(get_files_from_folder(gnss_folder))

    if not leo_files:
        raise FileNotFoundError("No LEO .bsp files found in leo_folder.")
    if not gnss_files:
        raise FileNotFoundError("No GNSS .bsp files found in gnss_folder.")
    if len(gnss_files) != 1:
        raise ValueError("In gnss_folder keep exactly ONE multi-sat GNSS .bsp.")

    gnss_path = gnss_files[0]
    os.makedirs(os.path.dirname(csv_file_path) or ".", exist_ok=True)

    total_rows = 0
    header_written = False

    with tempfile.TemporaryDirectory(prefix="gnssr_", dir=os.path.dirname(csv_file_path) or ".") as tmp_dir:
        for index, leo_path in enumerate(leo_files, start=1):
            leo_base = os.path.splitext(os.path.basename(leo_path))[0]
            tmp_csv_path = os.path.join(tmp_dir, f"{leo_base}.csv")

            print(f"[{index}/{len(leo_files)}] Processing {leo_base}")
            rows_saved = run_gnssr_single_leo(
                leo_path=leo_path,
                gnss_path=gnss_path,
                out_csv_path=tmp_csv_path,
                start_date=start_date,
                end_date=end_date,
                start_time_utc=start_time_utc,
                end_time_utc=end_time_utc,
                bbox=bbox,
                step_seconds=step_seconds,
                specular_err_max=specular_err_max,
                inc_min_deg=inc_min_deg,
                inc_max_deg=inc_max_deg,
                bistatic_max_deg=bistatic_max_deg,
                excess_max_km=excess_max_km,
                los_tol_km=los_tol_km,
                bbox_guided_search=bbox_guided_search,
                bbox_grid_deg=bbox_grid_deg,
                initial_latlon=initial_latlon,
            )

            if rows_saved == 0:
                continue

            with open(tmp_csv_path, "r", encoding="utf-8") as src, open(
                csv_file_path,
                "a" if header_written else "w",
                encoding="utf-8",
                newline="",
            ) as dst:
                for line_number, line in enumerate(src):
                    if header_written and line_number == 0:
                        continue
                    dst.write(line)

            header_written = True
            total_rows += rows_saved

    if not header_written:
        with open(csv_file_path, "w", encoding="utf-8", newline="") as dst:
            dst.write(
                "time_utc,GNSS,LEO,spec_lat_deg,spec_lon_deg,specular_err,incidence_deg,"
                "grazing_elev_deg,bistatic_deg,path_km,excess_m\n"
            )

    print(f"GNSS-R finished | rows={total_rows} | output={csv_file_path}")
    return total_rows


def run_gnssr_subprocess_per_leo(
    leo_folder,
    gnss_folder,
    output_folder,
    start_date,
    end_date,
    start_time_utc="00:00:00",
    end_time_utc="23:00:00",
    bbox=None,
    step_seconds=3600,
    specular_err_max=1.0,
    inc_min_deg=None,
    inc_max_deg=None,
    bistatic_max_deg=None,
    excess_max_km=None,
    los_tol_km=0.01,
    bbox_guided_search=True,
    bbox_grid_deg=0.25,
    initial_latlon=None,
    max_workers=8
):
    leo_files = sorted(get_files_from_folder(leo_folder))
    gnss_files = sorted(get_files_from_folder(gnss_folder))

    if not leo_files:
        raise FileNotFoundError("No LEO .bsp files found in leo_folder.")
    if not gnss_files:
        raise FileNotFoundError("No GNSS .bsp files found in gnss_folder.")
    if len(gnss_files) != 1:
        raise ValueError("In gnss_folder keep exactly ONE multi-sat GNSS .bsp.")

    gnss_path = gnss_files[0]
    os.makedirs(output_folder, exist_ok=True)

    print(f"Subprocess GNSS-R per LEO | max_workers={max_workers}")
    print(f"LEO kernels={len(leo_files)}")
    print(f"GNSS kernel: {os.path.basename(gnss_path)}")
    print(f"Output folder: {output_folder}")

    pending = list(leo_files)
    running = []
    done_count = 0

    while pending or running:
        while pending and len(running) < max_workers:
            leo_path = pending.pop(0)
            leo_base = os.path.splitext(os.path.basename(leo_path))[0]
            out_csv_path = os.path.join(output_folder, f"{leo_base}.csv")

            cmd = [
                sys.executable,
                os.path.abspath(__file__),
                "--worker-single-leo",
                "--leo-path", leo_path,
                "--gnss-path", gnss_path,
                "--out-csv", out_csv_path,
                "--start-date", start_date,
                "--end-date", end_date,
                "--start-time-utc", start_time_utc,
                "--end-time-utc", end_time_utc,
                "--step-seconds", str(step_seconds),
                "--specular-err-max", str(specular_err_max),
                "--los-tol-km", str(los_tol_km),
                "--bbox-guided-search", "1" if bbox_guided_search else "0",
                "--bbox-grid-deg", "None" if bbox_grid_deg is None else str(bbox_grid_deg),
                "--inc-min-deg", "None" if inc_min_deg is None else str(inc_min_deg),
                "--inc-max-deg", "None" if inc_max_deg is None else str(inc_max_deg),
                "--bistatic-max-deg", "None" if bistatic_max_deg is None else str(bistatic_max_deg),
                "--excess-max-km", "None" if excess_max_km is None else str(excess_max_km),
                "--initial-latlon", "None" if initial_latlon is None else f"{initial_latlon[0]},{initial_latlon[1]}",
            ]

            if bbox is not None:
                cmd += ["--bbox", ",".join(str(x) for x in bbox)]
            else:
                cmd += ["--bbox", "None"]

            proc = subprocess.Popen(cmd)
            running.append((proc, leo_path, out_csv_path))

        still_running = []
        for proc, leo_path, out_csv_path in running:
            ret = proc.poll()
            if ret is None:
                still_running.append((proc, leo_path, out_csv_path))
                continue

            done_count += 1
            leo_name = os.path.basename(leo_path)
            if ret == 0:
                print(f"[{done_count}/{len(leo_files)}] {leo_name} done | {out_csv_path}")
            else:
                print(f"[{done_count}/{len(leo_files)}] {leo_name} FAILED (exit={ret}) | {out_csv_path}")

        running = still_running
        time.sleep(1.0)

    print(f"ALL DONE | output_folder={output_folder}")


def _parse_optional_float(value):
    if value is None or value == "None":
        return None
    return float(value)


def _parse_optional_bbox(value):
    if value is None or value == "None":
        return None
    vals = [float(x) for x in value.split(",")]
    if len(vals) != 4:
        raise ValueError("bbox must have 4 comma-separated values")
    return tuple(vals)


def _parse_optional_latlon(value):
    if value is None or value == "None":
        return None
    vals = [float(x) for x in value.split(",")]
    if len(vals) != 2:
        raise ValueError("initial_latlon must have 2 comma-separated values")
    return tuple(vals)


if __name__ == "__main__":
    mp.freeze_support()

    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-single-leo", action="store_true")
    parser.add_argument("--leo-path")
    parser.add_argument("--gnss-path")
    parser.add_argument("--out-csv")
    parser.add_argument("--start-date")
    parser.add_argument("--end-date")
    parser.add_argument("--start-time-utc", default="00:00:00")
    parser.add_argument("--end-time-utc", default="23:00:00")
    parser.add_argument("--bbox", default="None")
    parser.add_argument("--step-seconds", type=int, default=3600)
    parser.add_argument("--specular-err-max", type=float, default=1.0)
    parser.add_argument("--inc-min-deg", default="None")
    parser.add_argument("--inc-max-deg", default="None")
    parser.add_argument("--bistatic-max-deg", default="None")
    parser.add_argument("--excess-max-km", default="None")
    parser.add_argument("--los-tol-km", type=float, default=0.01)
    parser.add_argument("--bbox-guided-search", default="1")
    parser.add_argument("--bbox-grid-deg", default="0.25")
    parser.add_argument("--initial-latlon", default="None")
    args = parser.parse_args()

    if args.worker_single_leo:
        bbox = _parse_optional_bbox(args.bbox)
        initial_latlon = _parse_optional_latlon(args.initial_latlon)

        rows_saved = run_gnssr_single_leo(
            leo_path=args.leo_path,
            gnss_path=args.gnss_path,
            out_csv_path=args.out_csv,
            start_date=args.start_date,
            end_date=args.end_date,
            start_time_utc=args.start_time_utc,
            end_time_utc=args.end_time_utc,
            bbox=bbox,
            step_seconds=args.step_seconds,
            specular_err_max=args.specular_err_max,
            inc_min_deg=_parse_optional_float(args.inc_min_deg),
            inc_max_deg=_parse_optional_float(args.inc_max_deg),
            bistatic_max_deg=_parse_optional_float(args.bistatic_max_deg),
            excess_max_km=_parse_optional_float(args.excess_max_km),
            los_tol_km=args.los_tol_km,
            bbox_guided_search=(args.bbox_guided_search == "1"),
            bbox_grid_deg=_parse_optional_float(args.bbox_grid_deg),
            initial_latlon=initial_latlon
        )
        print(f"WORKER DONE | rows={rows_saved} | {args.out_csv}")
    else:
        parser.error("Use main.py to run GNSS-R. --worker-single-leo is reserved for internal batch workers.")
