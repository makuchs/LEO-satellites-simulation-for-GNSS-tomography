"""
Occultation Analysis Module

This module loads SPICE kernels, reads satellite files for LEO and GNSS,
and processes radio occultations using SPICE functions.

Dependencies:
  - os, gc
  - numpy, pandas, spiceypy
"""

import os
import numpy as np
import pandas as pd
import spiceypy as spice

def _kernel_coverage_bounds(kernel_path, body_id):
    cover = spice.cell_double(200000)
    spice.spkcov(kernel_path, int(body_id), cover)
    if spice.wncard(cover) == 0:
        return None, None
    start, _ = spice.wnfetd(cover, 0)
    _, stop = spice.wnfetd(cover, spice.wncard(cover) - 1)
    return start, stop

def load_occultation_kernel_pool(leo_paths, gnss_paths):
    """
    Clears the kernel pool and loads leap seconds, planetary and Earth orientation
    kernels, then every LEO and GNSS SPK path given (e.g. .bsp files).

    Call once before iterating many (LEO, GNSS, time window) combinations.
    """
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    KERNELS_DIR = os.path.join(BASE_DIR, "kernels")

    spice.kclear()
    spice.furnsh(os.path.join(KERNELS_DIR, "lsk", "naif0012.tls"))
    spice.furnsh(os.path.join(KERNELS_DIR, "pck", "pck00011.tpc"))
    spice.furnsh(os.path.join(KERNELS_DIR, "spk", "de432s.bsp"))
    spice.furnsh(os.path.join(KERNELS_DIR, "pck", "earth_000101_241106_240813.bpc"))
    for p in leo_paths:
        spice.furnsh(p)
    for p in gnss_paths:
        spice.furnsh(p)


def load_spice_kernels(leo_path, gnss_path):
    """
    Loads the SPICE kernels for a single LEO and GNSS satellite pair (same pool layout
    as load_occultation_kernel_pool).
    """
    load_occultation_kernel_pool([leo_path], [gnss_path])

def load_time_kernel():
    """
    Loads only the time kernel from the repository's kernels folder.
    """
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    KERNELS_DIR = os.path.join(BASE_DIR, "kernels")
    
    spice.kclear()
    spice.furnsh(os.path.join(KERNELS_DIR, "lsk", "naif0012.tls"))

def get_files_from_folder(folder_path):
    """
    Returns a list of full file paths from the specified folder.
    """
    return [os.path.join(folder_path, file) for file in os.listdir(folder_path)
            if os.path.isfile(os.path.join(folder_path, file))]

def get_name_and_id_from_path(path):
    """
    Extracts the satellite name and ID from the file path.
    Assumes filename format: "ID_NAME.ext"
    """
    name_with_id = os.path.basename(path)
    file_name = os.path.splitext(name_with_id)[0]
    parts = file_name.split('_')
    sat_id = parts[0]
    name = parts[1] if len(parts) > 1 else "Unknown"
    return name, sat_id

def lon_lat_in_bbox(lon_deg, lat_deg, bbox):
    """
    True if (lon_deg, lat_deg) lies inside bbox.

    bbox: (min_lat, max_lat, min_lon, max_lon) in degrees, or None to skip the
    geographic filter (still requires finite tangent lon/lat).
    Missing coordinates (e.g. failed tangent computation) yield False.
    """
    if lon_deg is None or lat_deg is None:
        return False
    if bbox is None:
        return True
    min_lat, max_lat, min_lon, max_lon = bbox
    return (min_lat <= lat_deg <= max_lat) and (min_lon <= lon_deg <= max_lon)

def convert_to_lat_long(start, pos_leo, pos_gnss):
    """
    Converts state vectors to geographic latitude and longitude for both LEO and GNSS.
    """
    matrix = spice.pxform('J2000', 'IAU_EARTH', start)
    pos_leo_earth = spice.mxv(matrix, pos_leo)
    pos_gnss_earth = spice.mxv(matrix, pos_gnss)
    radii = spice.bodvrd('EARTH', 'RADII', 3)
    re = radii[1][0]
    rp = radii[1][2]
    f = (re - rp) / re           
    lon_leo, lat_leo, _ = spice.recgeo(pos_leo_earth, re, f)
    lon_gnss, lat_gnss, _ = spice.recgeo(pos_gnss_earth, re, f)
    lon_leo_deg = spice.convrt(lon_leo, 'RADIANS', 'DEGREES')
    lat_leo_deg = spice.convrt(lat_leo, 'RADIANS', 'DEGREES')
    lon_gnss_deg = spice.convrt(lon_gnss, 'RADIANS', 'DEGREES')
    lat_gnss_deg = spice.convrt(lat_gnss, 'RADIANS', 'DEGREES')

    return lon_leo_deg, lat_leo_deg, lon_gnss_deg, lat_gnss_deg

def is_within_fov(vel_vector, los_vector, fov_angle_deg=55):
    """
    Determines if a line-of-sight vector is within a given field-of-view angle.
    Considers both the ram and wake directions.
    """
    vel_norm = spice.vhat(vel_vector)
    los_norm = spice.vhat(los_vector)
    cos_angle = spice.vdot(vel_norm, los_norm)
    angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
    cos_angle_wake = spice.vdot(spice.vminus(vel_norm), los_norm)
    angle_deg_wake = np.degrees(np.arccos(np.clip(cos_angle_wake, -1.0, 1.0)))
    in_fov_ram = abs(angle_deg) <= fov_angle_deg
    in_fov_wake = abs(angle_deg_wake) <= fov_angle_deg
    
    return in_fov_ram or in_fov_wake

def process_occultation(leo_id, gnss_id, timestamp, gnss_name, leo_name,
                        results_list, fshape, target, fframe, abcorr, locus, observer, rayfrm, bbox):
    """
    Processes a single occultation event, extracting position information
    and checking if the event is within the field-of-view and inside the bounding box.
    """
    start_utc = spice.timout(timestamp, "YYYY-MM-DD HR:MN:SC ::UTC")
    pos_leo, _ = spice.spkpos(leo_id, timestamp, 'J2000', 'NONE', 'EARTH')
    pos_gnss, _ = spice.spkpos(gnss_id, timestamp, 'J2000', 'NONE', 'EARTH')
    lon_leo_deg, lat_leo_deg, lon_gnss_deg, lat_gnss_deg = convert_to_lat_long(timestamp, pos_leo, pos_gnss)
    dvec = np.subtract(pos_gnss, pos_leo)
    
    try:
        tangent_point, _, _, _, _, _ = spice.tangpt(
            fshape, target, timestamp, fframe, abcorr, locus, observer, rayfrm, dvec
        )
        _, lon, lat = spice.reclat(tangent_point)
        lon_deg = spice.convrt(lon, 'RADIANS', 'DEGREES')
        lat_deg = spice.convrt(lat, 'RADIANS', 'DEGREES')
    except Exception as e:
        print(f"Error: {e}")
        lon_deg, lat_deg = None, None

    state_leo, _ = spice.spkezr(leo_id, timestamp, 'J2000', 'NONE', 'EARTH')
    pos_leo_fov, vel_leo_fov = state_leo[:3], state_leo[3:6]
    state_gnss, _ = spice.spkezr(gnss_id, timestamp, 'J2000', 'NONE', 'EARTH')
    pos_gnss_fov = state_gnss[:3]
    los_vector = spice.vsub(pos_gnss_fov, pos_leo_fov)
    in_fov = is_within_fov(vel_leo_fov, los_vector)
    in_bbox = lon_lat_in_bbox(lon_deg, lat_deg, bbox)

    if in_fov and in_bbox:
        results_list.append({
            'GNSS': gnss_name,
            'LEO': leo_name,
            'DATE': start_utc,
            'Pos_GNSS': pos_gnss,
            'Pos_LEO': pos_leo,
            'T_lon_deg': lon_deg,
            'T_lat_deg': lat_deg,
            'lon_leo_deg': lon_leo_deg,
            'lat_leo_deg': lat_leo_deg,
            'lon_gnss_deg': lon_gnss_deg,
            'lat_gnss_deg': lat_gnss_deg,
        })

def find_RO_occultations(leo_path, gnss_path, et1, et2, bbox, load_kernels=True):
    """
    Finds radio occultations between a LEO and GNSS satellite within the given ephemeris time window.

    bbox: (min_lat, max_lat, min_lon, max_lon) in degrees for tangent-point filtering,
    or None to disable geographic filtering.

    If load_kernels is True (default), clears the pool and loads this pair's SPKs plus
    common kernels. Set load_kernels=False when the pool was already loaded (e.g. via
    load_occultation_kernel_pool for all LEO/GNSS files).

    Returns a pandas DataFrame with occultation details.
    """
    if load_kernels:
        load_spice_kernels(leo_path, gnss_path)
    leo_name, leo_id = get_name_and_id_from_path(leo_path)
    gnss_name, gnss_id = get_name_and_id_from_path(gnss_path)

    leo_start, leo_stop = _kernel_coverage_bounds(leo_path, leo_id)
    gnss_start, gnss_stop = _kernel_coverage_bounds(gnss_path, gnss_id)
    et1_clamped = max(et1, leo_start, gnss_start)
    et2_clamped = min(et2, leo_stop, gnss_stop)
    if et1_clamped >= et2_clamped:
        return pd.DataFrame([])
    
    rayfrm = "J2000"
    locus = "TANGENT POINT"
    occtype, front, fshape, fframe = 'any', 'EARTH', 'ELLIPSOID', 'IAU_EARTH'
    back, bshape, bframe, observer, abcorr = gnss_id, 'POINT', 'IAU_EARTH', leo_id, 'NONE'
    step, MAXWIN = 30, 50000
    confine = spice.cell_double(2 * MAXWIN)
    spice.wninsd(et1_clamped, et2_clamped, confine)
    result = spice.cell_double(MAXWIN)
    spice.gfoclt(occtype, front, fshape, fframe, back, bshape, bframe, abcorr, observer, step, confine, result)
    target = front
    num_intervals = spice.wncard(result)
    results_list = []
    for i in range(num_intervals):
        start, stop = spice.wnfetd(result, i)
        process_occultation(leo_id, gnss_id, start, gnss_name, leo_name,
                            results_list, fshape, target, fframe, abcorr, locus, observer, rayfrm, bbox)
        process_occultation(leo_id, gnss_id, stop, gnss_name, leo_name,
                            results_list, fshape, target, fframe, abcorr, locus, observer, rayfrm, bbox)
    occultations_df = pd.DataFrame(results_list)
    return occultations_df

def run_occultation(leo_folder, gnss_folder, csv_file_path, start_date, end_date, bbox):
    """
    Runs the occultation analysis over the specified date range using the given folders.

    bbox: (min_lat, max_lat, min_lon, max_lon) in degrees for tangent-point filtering,
    or None to keep all tangent points with valid coordinates (FOV filter still applies).
    """
    dates = pd.date_range(start=start_date, end=end_date, freq='D')

    leo_files = get_files_from_folder(leo_folder)
    gnss_files = get_files_from_folder(gnss_folder)

    if not os.path.exists(csv_file_path):
        with open(csv_file_path, 'w') as f:
            f.write('')

    load_occultation_kernel_pool(leo_files, gnss_files)
    try:
        for date in dates:
            et1 = spice.str2et(date.strftime('%Y %m %d 00:00:00 UTC'))
            et2 = spice.str2et(date.strftime('%Y %m %d 23:59:59 UTC'))
            for gnss_path in gnss_files:
                for leo_path in leo_files:
                    try:
                        print(
                            f"Processing {os.path.basename(leo_path)} and {os.path.basename(gnss_path)} "
                            f"for date {date.strftime('%Y-%m-%d')}..."
                        )
                        temp_df = find_RO_occultations(
                            leo_path, gnss_path, et1, et2, bbox, load_kernels=False
                        )
                        temp_df.to_csv(
                            csv_file_path,
                            mode='a',
                            header=not os.path.getsize(csv_file_path) > 0,
                            index=False,
                        )
                    except Exception as e:
                        print(f"Error processing {gnss_path} on {date.strftime('%Y-%m-%d')}: {e}")
    finally:
        spice.kclear()
                    
if __name__ == '__main__':
    raise SystemExit("Use main.py with MODE='occultation', or import run_occultation() from Python.")
