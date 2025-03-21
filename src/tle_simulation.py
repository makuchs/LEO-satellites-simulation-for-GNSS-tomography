"""
TLE Simulation and SPICE Kernel Generation

This module loads TLE records, propagates satellite states using sgp4, converts state vectors
from TEME to J2000 coordinates, writes an SPK kernel using spiceypy, and optionally writes 
simulated positions converted to ECEF coordinates to a CSV file.

Dependencies:
  - math, datetime, os, re, csv
  - numpy, sgp4, skyfield, astropy, spiceypy
"""

import math
import datetime
import os
import re
import csv
import numpy as np

from sgp4.api import Satrec
import skyfield.sgp4lib as sgp4lib
from astropy import coordinates as coord, units as u
from astropy.time import Time
import spiceypy as spice  

def load_spice_kernels():
    """
    Loads the required SPICE kernels from the repository's kernels folder.
    """
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    KERNELS_DIR = os.path.join(BASE_DIR, "kernels")
    
    spice.kclear()
    spice.furnsh(os.path.join(KERNELS_DIR, "lsk", "naif0012.tls"))
    spice.furnsh(os.path.join(KERNELS_DIR, "pck", "pck00011.tpc"))
    spice.furnsh(os.path.join(KERNELS_DIR, "spk", "de432s.bsp"))
    spice.furnsh(os.path.join(KERNELS_DIR, "pck", "earth_000101_241106_240813.bpc"))

def load_tle_records(filename):
    """
    Reads a TLE text file and returns a list of TLE records.
    Each record is a tuple: (name, line1, line2).
    Assumes the file is organized in groups of three nonempty lines.
    """
    tle_records = []
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    if len(lines) % 3 != 0:
        raise ValueError("TLE file format error: total nonempty lines not a multiple of 3.")
    for i in range(0, len(lines), 3):
        tle_records.append((lines[i], lines[i+1], lines[i+2]))
    return tle_records

def propagate_to_next_whole_second(satellite, jd=None, fr=None):
    """
    Propagates a satellite's position to the next whole second.
    If jd and fr are not provided, uses the satellite's epoch.
    """
    if jd is None or fr is None:
        jd = satellite.jdsatepoch
        fr = satellite.jdsatepochF

    current_seconds = fr * 86400.0
    next_seconds = math.ceil(current_seconds)

    if next_seconds >= 86400:
        next_jd = jd + 1
        next_fr = 0.0
    else:
        next_jd = jd
        next_fr = next_seconds / 86400.0

    e, r, v = satellite.sgp4(next_jd, next_fr)
    return satellite, e, r, v, next_jd, next_fr

def simulate_tle_states(tle_records, timestep_seconds=30):
    """
    Simulate satellite positions for each TLE record.
    
    Returns:
      simulation_times: List of simulation times (Julian Dates)
      simulation_positions: List of TEME positions
      simulation_velocities: List of TEME velocities
    """
    simulation_times = []
    simulation_positions = []
    simulation_velocities = []
    num_records = len(tle_records)
    
    for i in range(num_records):
        name, line1, line2 = tle_records[i]
        print(f"Processing record {i+1}/{num_records}: {name}")
        sat = Satrec.twoline2rv(line1, line2)
        
        sat, err, pos, vel, current_jd, current_fr = propagate_to_next_whole_second(sat)
        current_time = current_jd + current_fr

        if i < num_records - 1:
            next_name, next_line1, next_line2 = tle_records[i+1]
            next_sat = Satrec.twoline2rv(next_line1, next_line2)
            end_time = next_sat.jdsatepoch + next_sat.jdsatepochF
        else:
            end_time = current_time + timestep_seconds / 86400.0
        
        while current_time < end_time:
            e, r, v = sat.sgp4(current_jd, current_fr)
            simulation_times.append(current_time)
            simulation_positions.append(r)
            simulation_velocities.append(v)
            
            current_time += timestep_seconds / 86400.0
            current_jd = math.floor(current_time)
            current_fr = current_time - current_jd
        
        final_jd = math.floor(end_time)
        final_fr = end_time - final_jd
        e, r, v = sat.sgp4(final_jd, final_fr)
        simulation_times.append(end_time)
        simulation_positions.append(r)
        simulation_velocities.append(v)
        
    return simulation_times, simulation_positions, simulation_velocities

def write_positions_csv_ecef(csv_filename, times, positions, velocities, time_offset=19):
    """
    Converts TEME state vectors to ECEF (ITRF) using teme_to_ecef and writes
    the positions to a CSV file. The CSV file will have columns: Date, X, Y, Z.
    """
    j2000 = datetime.datetime(2000, 1, 1, 12, 0, 0)
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Date", "X", "Y", "Z"])
        for t, p, v in zip(times, positions, velocities):
            adjusted_time, p_ecef, _ = teme_to_ecef(t, p, v, time_offset=time_offset)
            dt = j2000 + datetime.timedelta(days=(adjusted_time - 2451545.0))
            date_str = dt.strftime("%Y %m %d %H %M %S") + ".00"
            writer.writerow([date_str,
                             f"{p_ecef[0]:16.6f}",
                             f"{p_ecef[1]:16.6f}",
                             f"{p_ecef[2]:16.6f}"])
    print(f"CSV file '{csv_filename}' with ECEF positions created successfully.")

def teme_to_ecef(time_jd, p_teme, v_teme, time_offset=19):
    """
    Converts a state vector from TEME to ECEF (ITRF) using Skyfield.
    Returns adjusted time, ECEF position, and ECEF velocity.
    """
    # Convert offset from seconds to days:
    offset_days = time_offset / 86400.0
    adjusted_time = time_jd + offset_days
    v_teme_km_day = np.asarray(v_teme) * 86400.0
    p_ecef, v_ecef_km_day = sgp4lib.TEME_to_ITRF(time_jd, np.asarray(p_teme), v_teme_km_day)
    v_ecef = np.asarray(v_ecef_km_day) / 86400.0
    return adjusted_time, p_ecef.tolist(), v_ecef.tolist()

def teme_to_j2000(time_jd, p_teme, v_teme):
    """
    Converts a state vector from TEME to J2000 inertial (GCRS) coordinates using Astropy.
    """
    v_teme_km_day = np.asarray(v_teme) * 86400.0
    p_itrs, v_itrs_km_day = sgp4lib.TEME_to_ITRF(time_jd, np.asarray(p_teme), v_teme_km_day)
    v_itrs = np.asarray(v_itrs_km_day) / 86400.0

    date = datetime.datetime(2000, 1, 1, 12, 0, 0) + datetime.timedelta(days=(time_jd - 2451545.0))
    obstime = Time(date)
    
    itrs = coord.ITRS(x=p_itrs[0]*u.km, y=p_itrs[1]*u.km, z=p_itrs[2]*u.km,
                      v_x=v_itrs[0]*u.km/u.s, v_y=v_itrs[1]*u.km/u.s, v_z=v_itrs[2]*u.km/u.s,
                      obstime=obstime)
    
    gcrs = itrs.transform_to(coord.GCRS(obstime=obstime))
    
    p_j2000 = gcrs.cartesian.xyz.value.tolist()
    v_j2000 = gcrs.velocity.d_xyz.value.tolist()
    
    return p_j2000, v_j2000

def batch_teme_to_j2000(jd_list, p_teme_list, v_teme_list):
    """
    Batch conversion of TEME state vectors to J2000 (GCRS) coordinates.
    """
    p_itrs_list = []
    v_itrs_list = []
    for jd, p_teme, v_teme in zip(jd_list, p_teme_list, v_teme_list):
        v_teme_km_day = np.asarray(v_teme) * 86400.0
        p_itrs, v_itrs_km_day = sgp4lib.TEME_to_ITRF(jd, np.asarray(p_teme), v_teme_km_day)
        p_itrs_list.append(p_itrs)
        v_itrs_list.append(np.asarray(v_itrs_km_day) / 86400.0)
    
    p_itrs_arr = np.array(p_itrs_list)
    v_itrs_arr = np.array(v_itrs_list)
    
    j2000_epoch = datetime.datetime(2000, 1, 1, 12, 0, 0)
    dt_list = [j2000_epoch + datetime.timedelta(days=(jd - 2451545.0)) for jd in jd_list]
    obstimes = Time(dt_list)
    
    rep = coord.CartesianRepresentation(x=p_itrs_arr[:,0]*u.km,
                                        y=p_itrs_arr[:,1]*u.km,
                                        z=p_itrs_arr[:,2]*u.km)
    
    diff = coord.CartesianDifferential(d_x=v_itrs_arr[:,0]*u.km/u.s,
                                       d_y=v_itrs_arr[:,1]*u.km/u.s,
                                       d_z=v_itrs_arr[:,2]*u.km/u.s)
    
    itrs_coords = coord.ITRS(rep.with_differentials(diff), obstime=obstimes)
    
    gcrs_coords = itrs_coords.transform_to(coord.GCRS(obstime=obstimes))
    
    p_j2000 = np.column_stack([gcrs_coords.cartesian.x.value,
                               gcrs_coords.cartesian.y.value,
                               gcrs_coords.cartesian.z.value])
    v_j2000 = np.column_stack([gcrs_coords.velocity.d_x.value,
                               gcrs_coords.velocity.d_y.value,
                               gcrs_coords.velocity.d_z.value])
    
    return p_j2000, v_j2000

def write_j2000_spice_kernel(kernel_filename, times, positions, velocities, tle_filename):
    """
    Converts TEME simulation states to J2000 using batch processing and writes an SPK kernel using spiceypy.
    The satellite NORAD ID is extracted from the TLE filename.
    """
    load_spice_kernels()
    
    base = os.path.basename(tle_filename)
    match = re.match(r"(\d+)", base)
    if match:
        sat_norad = match.group(1)
    else:
        sat_norad = "UNKNOWN"
    
    def jd_to_utc_str(jd):
        j2000 = datetime.datetime(2000, 1, 1, 12, 0, 0)
        dt = j2000 + datetime.timedelta(days=(jd - 2451545.0))
        return dt.strftime("%Y-%m-%dT%H:%M:%S")
    
    j2000_positions_arr, j2000_velocities_arr = batch_teme_to_j2000(times, positions, velocities)
    
    et_times = []
    for jd in times:
        utc_str = jd_to_utc_str(jd)
        et = spice.utc2et(utc_str)
        et_times.append(et)
    et_times = np.array(et_times)
    
    states_matrix = np.hstack((j2000_positions_arr, j2000_velocities_arr)).tolist()
    
    segid = f"SPK_SEGMENT_{sat_norad}"
    handle = spice.spkopn(kernel_filename, f"SPK Kernel for satellite {sat_norad}", 0)
    
    spice.spkw08(handle, int(sat_norad), 399, "J2000", et_times[0], et_times[-1],
                 segid, 7, len(et_times), states_matrix, et_times[0], 30)
    
    spice.spkcls(handle)
    print(f"SPK kernel '{kernel_filename}' created successfully for satellite {sat_norad}.")

def remove_duplicate_epochs(times, positions, velocities, tol=1e-9):
    """
    Removes duplicate epochs from the lists 'times', 'positions', and 'velocities'.
    If the difference between consecutive epochs is less than tol, only the first occurrence is kept.
    
    Returns:
      new_times, new_positions, new_velocities
    """
    if not times:
        return times, positions, velocities

    new_times = [times[0]]
    new_positions = [positions[0]]
    new_velocities = [velocities[0]]
    
    for t, p, v in zip(times[1:], positions[1:], velocities[1:]):
        if abs(t - new_times[-1]) < tol:
            # Duplicate epoch: skip it, keeping the first occurrence
            continue
        else:
            new_times.append(t)
            new_positions.append(p)
            new_velocities.append(v)
    
    return new_times, new_positions, new_velocities


def run_simulation(tle_filename, csv_output=None, output_folder="output"):
    """
    Runs the TLE simulation and SPICE kernel generation.
    If csv_output is not provided, it creates a CSV file with the same base name as the TLE file appended with '_ecef.csv'
    in the specified output folder.
    """
    # Ensure the output folder exists.
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    print("Loading TLE records...")
    tle_list = load_tle_records(tle_filename)
    
    print("Starting simulation of satellite states (TEME)...")
    times, positions, velocities = simulate_tle_states(tle_list, timestep_seconds=1)
    times, positions, velocities = remove_duplicate_epochs(times, positions, velocities)

    base_name = os.path.splitext(os.path.basename(tle_filename))[0]
    kernel_filename = os.path.join(output_folder, base_name + ".bsp")
    print("Writing J2000 SPICE kernel...(be patient!)")
    write_j2000_spice_kernel(kernel_filename, times, positions, velocities, tle_filename)
    
    if not csv_output:
        csv_output = os.path.join(output_folder, base_name + "_ecef.csv")
    
    print("Writing ECEF positions to CSV for comparison...")
    write_positions_csv_ecef(csv_output, times, positions, velocities)
    
    print("Done.")

if __name__ == '__main__':
    # For direct testing, adjust the TLE file path as needed.
    run_simulation(r"C:\Users\User\Desktop\463188_COSMIC22.txt")
