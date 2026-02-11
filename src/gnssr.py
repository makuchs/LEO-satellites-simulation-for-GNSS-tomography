"""
GNSS-R (GNSS Reflectometry) Geometry Simulation Module

This module provides a SPICE-based workflow to simulate bistatic GNSS-R geometry between
GNSS transmitters (Tx) and LEO receivers (Rx). It computes an approximate specular point
on the Earth ellipsoid (from SPICE body radii), filters events by a geographic bounding
box, checks visibility from Tx and Rx to the specular point, and writes results to CSV.

Key upgrades in this version:
  - Optional time-of-day window: start_time_utc / end_time_utc (so you can simulate only the .nc window).
  - Optional bbox-guided specular search (grid search inside bbox) to avoid converging to a different local solution.
  - LOS tolerance relaxed (default 10 m) to avoid rejecting valid geometries due to numerical mismatch.
"""

import os
import numpy as np
import pandas as pd
import spiceypy as spice


def load_spice_kernels(leo_path, gnss_path):
    """
    Load base SPICE kernels and two SPK kernels (LEO and GNSS).
    """
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    kernels_dir = os.path.join(base_dir, "kernels")

    spice.kclear()
    spice.furnsh(os.path.join(kernels_dir, "lsk", "naif0012.tls"))
    spice.furnsh(os.path.join(kernels_dir, "pck", "pck00011.tpc"))
    spice.furnsh(os.path.join(kernels_dir, "spk", "de432s.bsp"))

    earth_bpc = os.path.join(kernels_dir, "pck", "earth_000101_241106_240813.bpc")
    if os.path.exists(earth_bpc):
        spice.furnsh(earth_bpc)

    spice.furnsh(leo_path)
    spice.furnsh(gnss_path)


def load_time_kernel():
    """
    Load only the leapseconds kernel (LSK).
    """
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    kernels_dir = os.path.join(base_dir, "kernels")
    spice.kclear()
    spice.furnsh(os.path.join(kernels_dir, "lsk", "naif0012.tls"))


def get_files_from_folder(folder_path):
    """
    Return a list of .bsp files from the given folder as full paths.
    """
    return [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if os.path.isfile(os.path.join(folder_path, f)) and f.lower().endswith(".bsp")
    ]


def get_name_and_id_from_path(path):
    """
    Extract satellite name and ID from a kernel filename "<ID>_<NAME>.bsp" or "<ID>.bsp".
    """
    name_with_id = os.path.basename(path)
    file_name = os.path.splitext(name_with_id)[0]
    parts = file_name.split("_")
    sat_id = parts[0]
    name = parts[1] if len(parts) > 1 else parts[0]
    return name, sat_id


def _earth_radii_and_flattening():
    """
    Return Earth radii (km) and flattening used by SPICE coordinate helpers.
    """
    radii = spice.bodvrd("EARTH", "RADII", 3)[1]
    re = float(radii[0])
    rp = float(radii[2])
    f = (re - rp) / re
    return re, rp, f


def _state_in_frame(target_id, et, out_frame="J2000", abcorr="NONE", obs="EARTH"):
    """
    Get target state (position, velocity) in the requested frame.
    """
    st, _ = spice.spkezr(str(target_id), float(et), out_frame, abcorr, obs)
    pos = np.array(st[:3], dtype=float)
    vel = np.array(st[3:6], dtype=float)
    return pos, vel


def _transform_state(pos, vel, et, from_frame="J2000", to_frame="IAU_EARTH"):
    """
    Transform a 6D state (pos, vel) between frames using SPICE state transformation.
    """
    xform = spice.sxform(from_frame, to_frame, float(et))
    st6 = np.hstack([pos, vel])
    st6_out = spice.mxvg(xform, st6)
    return np.array(st6_out[:3], dtype=float), np.array(st6_out[3:6], dtype=float)


def _ecef_latlon_alt_from_pos(pos_ecef):
    """
    Convert Earth-fixed Cartesian position (km) to geodetic lat/lon/alt.
    """
    re, _, f = _earth_radii_and_flattening()
    lon, lat, alt = spice.recgeo(pos_ecef, re, f)
    lat_deg = float(spice.convrt(lat, "RADIANS", "DEGREES"))
    lon_deg = float(spice.convrt(lon, "RADIANS", "DEGREES"))
    return lat_deg, lon_deg, float(alt)


def _wrap_lon_deg(lon_deg):
    """
    Wrap longitude to [-180, 180).
    """
    return (float(lon_deg) + 180.0) % 360.0 - 180.0


def is_in_bbox(lat_deg, lon_deg, bbox):
    """
    Check if a (lat, lon) point lies inside bbox=(min_lat, max_lat, min_lon, max_lon).
    Handles dateline-crossing bboxes.
    """
    min_lat, max_lat, min_lon, max_lon = bbox

    lat = float(lat_deg)
    lon = _wrap_lon_deg(float(lon_deg))
    min_lon = _wrap_lon_deg(float(min_lon))
    max_lon = _wrap_lon_deg(float(max_lon))

    if not (float(min_lat) <= lat <= float(max_lat)):
        return False

    if min_lon <= max_lon:
        return (min_lon <= lon <= max_lon)

    return (lon >= min_lon) or (lon <= max_lon)


def _surface_point_from_latlon(lat_deg, lon_deg):
    """
    Convert geodetic lat/lon (deg) to a point on the reference ellipsoid (km, Earth-fixed).
    """
    re, _, f = _earth_radii_and_flattening()
    lat = spice.convrt(float(lat_deg), "DEGREES", "RADIANS")
    lon = spice.convrt(float(lon_deg), "DEGREES", "RADIANS")
    return np.array(spice.pgrrec("EARTH", lon, lat, 0.0, re, f), dtype=float)


def _surface_normal_ecef(surfpt_ecef):
    """
    Compute outward surface normal at an ellipsoid surface point (Earth-fixed, km).
    """
    re, rp, _ = _earth_radii_and_flattening()
    return np.array(spice.surfnm(re, re, rp, surfpt_ecef), dtype=float)


def _specular_error(tx_ecef, rx_ecef, s_ecef):
    """
    Evaluate an error measure for the specular reflection condition at a surface point.
    """
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
    """
    Find an approximate specular point on the Earth ellipsoid using a local multi-resolution search
    in geodetic latitude/longitude.
    """
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
    """
    Find an approximate specular point by first searching a grid inside a bbox, then refining locally.
    """
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
    best_lat = float(lats[len(lats)//2]) if len(lats) else 0.0
    best_lon = float(_wrap_lon_deg(lons[len(lons)//2])) if len(lons) else 0.0
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
    """
    Check whether the observer has a clear line of sight to the surface point on the ellipsoid.
    """
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
    """
    Compute incidence angle at the specular point (degrees) relative to the surface normal.
    """
    n = _surface_normal_ecef(s_ecef)
    v = tx_ecef - s_ecef
    nv = np.linalg.norm(v)
    if nv < 1e-12:
        return float("nan")
    vhat = v / nv
    c = float(np.clip(np.dot(n, vhat), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def bistatic_angle_deg(tx_ecef, rx_ecef, s_ecef):
    """
    Compute bistatic angle at the specular point (degrees) between incoming and outgoing rays.
    """
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
    """
    Compute bistatic path length (km) and excess path (km) relative to direct Tx->Rx.
    """
    d1 = float(np.linalg.norm(tx_ecef - s_ecef))
    d2 = float(np.linalg.norm(rx_ecef - s_ecef))
    ddir = float(np.linalg.norm(tx_ecef - rx_ecef))
    return (d1 + d2), (d1 + d2 - ddir)


def process_gnssr_epoch(
    leo_id,
    gnss_id,
    et,
    leo_name,
    gnss_name,
    bbox=None,
    specular_err_max=0.15,
    inc_max_deg=75.0,
    bistatic_max_deg=60.0,
    excess_max_km=200.0,
    los_tol_km=0.01,
    bbox_guided_search=False,
    bbox_grid_deg=1.0,
    initial_latlon=None
):
    """
    Process one epoch and return a result dict if the GNSS-R geometry passes filters.
    """
    tx_pos_j2k, tx_vel_j2k = _state_in_frame(gnss_id, et, "J2000", "NONE", "EARTH")
    rx_pos_j2k, rx_vel_j2k = _state_in_frame(leo_id, et, "J2000", "NONE", "EARTH")

    tx_ecef, _ = _transform_state(tx_pos_j2k, tx_vel_j2k, et, "J2000", "IAU_EARTH")
    rx_ecef, _ = _transform_state(rx_pos_j2k, rx_vel_j2k, et, "J2000", "IAU_EARTH")

    if bbox_guided_search and (bbox is not None):
        s_ecef, _, _, s_err = find_specular_point_in_bbox(
            tx_ecef, rx_ecef, bbox=bbox, grid_deg=bbox_grid_deg
        )
    else:
        s_ecef, _, _, s_err = find_specular_point(
            tx_ecef, rx_ecef, initial_latlon=initial_latlon
        )

    if not np.isfinite(s_err) or float(s_err) > float(specular_err_max):
        return None

    s_lat, s_lon, s_alt_km = _ecef_latlon_alt_from_pos(s_ecef)
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

    if inc_max_deg is not None and np.isfinite(inc) and float(inc) > float(inc_max_deg):
        return None
    if bistatic_max_deg is not None and np.isfinite(bist) and float(bist) > float(bistatic_max_deg):
        return None
    if excess_max_km is not None and np.isfinite(excess_km) and abs(float(excess_km)) > float(excess_max_km):
        return None

    dt_utc = spice.timout(float(et), "YYYY-MM-DD HR:MN:SC ::UTC")

    return {
        "time_utc": dt_utc,
        "gnss_sat": gnss_name,
        "leo_sat": leo_name,
        "spec_lat_deg": float(s_lat),
        "spec_lon_deg": float(s_lon),
        "spec_alt_km": float(s_alt_km),
        "specular_err": float(s_err),
        "incidence_deg": float(inc),
        "bistatic_deg": float(bist),
        "path_km": float(path_km),
        "excess_km": float(excess_km),
    }


def run_gnssr(
    leo_folder,
    gnss_folder,
    csv_file_path,
    start_date,
    end_date,
    start_time_utc="00:00:00",
    end_time_utc="23:59:59",
    bbox=None,
    step_seconds=30,
    specular_err_max=0.15,
    inc_max_deg=75.0,
    bistatic_max_deg=60.0,
    excess_max_km=200.0,
    los_tol_km=0.01,
    bbox_guided_search=False,
    bbox_grid_deg=1.0,
    initial_latlon=None
):
    """
    Run GNSS-R geometry simulation over a date range for all LEO x GNSS kernel pairs and write results to CSV.
    """
    dates = pd.date_range(start=start_date, end=end_date, freq="D")
    leo_files = get_files_from_folder(leo_folder)
    gnss_files = get_files_from_folder(gnss_folder)

    header = [
        "time_utc",
        "gnss_sat",
        "leo_sat",
        "spec_lat_deg",
        "spec_lon_deg",
        "spec_alt_km",
        "specular_err",
        "incidence_deg",
        "bistatic_deg",
        "path_km",
        "excess_km",
    ]

    os.makedirs(os.path.dirname(csv_file_path) or ".", exist_ok=True)

    if not (os.path.exists(csv_file_path) and os.path.getsize(csv_file_path) > 0):
        with open(csv_file_path, "w", newline="") as f:
            f.write(",".join(header) + "\n")

    for date in dates:
        load_time_kernel()
        d = date.strftime("%Y %m %d")
        et1 = spice.str2et(f"{d} {start_time_utc} UTC")
        et2 = spice.str2et(f"{d} {end_time_utc} UTC")
        epochs = np.arange(float(et1), float(et2) + 1e-9, float(step_seconds), dtype=float)

        for gnss_path in gnss_files:
            for leo_path in leo_files:
                try:
                    load_spice_kernels(leo_path, gnss_path)
                    leo_name, leo_id = get_name_and_id_from_path(leo_path)
                    gnss_name, gnss_id = get_name_and_id_from_path(gnss_path)

                    rows = []
                    for et in epochs:
                        row = process_gnssr_epoch(
                            leo_id=leo_id,
                            gnss_id=gnss_id,
                            et=float(et),
                            leo_name=leo_name,
                            gnss_name=gnss_name,
                            bbox=bbox,
                            specular_err_max=specular_err_max,
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

                    if rows:
                        df = pd.DataFrame(rows, columns=header)
                        df.to_csv(csv_file_path, mode="a", header=False, index=False)

                except Exception as e:
                    print(
                        f"GNSS-R error {os.path.basename(leo_path)} {os.path.basename(gnss_path)} "
                        f"{date.strftime('%Y-%m-%d')}: {e}"
                    )
                finally:
                    try:
                        spice.kclear()
                    except Exception:
                        pass


if __name__ == "__main__":
    run_gnssr(
        leo_folder=r"kernels/leo",
        gnss_folder=r"kernels/gnss",
        csv_file_path=r"output/gnssr_results.csv",
        start_date="2023-10-01",
        end_date="2023-10-01",
        start_time_utc="01:20:00",
        end_time_utc="01:25:10",
        bbox=(-4.0, 2.5, 72.5, 83.3),
        step_seconds=1,
        specular_err_max=1.0,
        inc_max_deg=None,
        bistatic_max_deg=None,
        excess_max_km=None,
        los_tol_km=0.01,
        bbox_guided_search=True,
        bbox_grid_deg=0.5,
        initial_latlon=None
    )
