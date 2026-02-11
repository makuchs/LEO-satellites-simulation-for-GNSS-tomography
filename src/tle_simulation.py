"""
TLE Simulation and SPICE Kernel Generation

This module loads TLE records, propagates satellite states using sgp4, converts state vectors
from TEME to J2000 coordinates, writes an SPK kernel using spiceypy, and optionally writes 
simulated positions converted to ECEF coordinates to a CSV file.

It can also generate SPICE kernels from SP3 precise ephemeris files by:
- parsing positions in ITRF/ECEF,
- resampling to a fixed step,
- estimating velocities via finite differences,
- transforming ITRF -> GCRS/J2000 with Astropy,
- writing an SPK (type 8) via SpiceyPy.

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
            continue
        else:
            new_times.append(t)
            new_positions.append(p)
            new_velocities.append(v)

    return new_times, new_positions, new_velocities

def parse_sp3_positions(sp3_file, sat_id, unit_is_km=True):
    """
    Parse an SP3 file and return epochs (datetime list) and ITRF/ECEF positions (Nx3, km).
    Only records for the given sat_id are returned (e.g., 'PG01', 'PL99').
    """
    epochs = []
    pos_ecef = []
    current_epoch = None

    with open(sp3_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith('*'):
                tokens = line[1:].strip().split()
                if len(tokens) < 6:
                    current_epoch = None
                    continue
                year = int(tokens[0]); month = int(tokens[1]); day = int(tokens[2])
                hour = int(tokens[3]); minute = int(tokens[4]); second = float(tokens[5])
                sec_int = int(second)
                micro = int(round((second - sec_int) * 1e6))
                current_epoch = datetime.datetime(year, month, day, hour, minute, sec_int, micro)
                continue

            if current_epoch is None:
                continue

            tokens = line.split()
            if not tokens:
                continue

            if tokens[0] == sat_id:
                x = float(tokens[1]); y = float(tokens[2]); z = float(tokens[3])
                r = np.array([x, y, z], dtype=float)
                if not unit_is_km:
                    r = r / 1000.0
                epochs.append(current_epoch)
                pos_ecef.append(r)

    if len(epochs) < 2:
        raise ValueError(f"Not enough epochs for {sat_id} in {sp3_file} (found {len(epochs)}).")

    return epochs, np.vstack(pos_ecef)

def resample_positions(epochs_dt, pos_km, step_seconds=1):
    """
    Resample positions to a constant time step using linear interpolation.
    Returns (epochs_dt_new, pos_km_new).
    """
    t0 = epochs_dt[0]
    t_sec = np.array([(e - t0).total_seconds() for e in epochs_dt], dtype=float)
    t_new = np.arange(t_sec[0], t_sec[-1] + 1e-9, float(step_seconds), dtype=float)
    pos_new = np.empty((len(t_new), 3), dtype=float)
    for i in range(3):
        pos_new[:, i] = np.interp(t_new, t_sec, pos_km[:, i])
    epochs_new = [t0 + datetime.timedelta(seconds=float(s)) for s in t_new]
    return epochs_new, pos_new

def finite_difference_velocities(pos_km, step_seconds=1):
    """
    Compute velocities (km/s) from positions (km) using central differences.
    """
    dt = float(step_seconds)
    v = np.empty_like(pos_km)
    v[1:-1] = (pos_km[2:] - pos_km[:-2]) / (2.0 * dt)
    v[0] = (pos_km[1] - pos_km[0]) / dt
    v[-1] = (pos_km[-1] - pos_km[-2]) / dt
    return v

def itrf_to_gcrs(epochs_dt, pos_km, vel_km_s):
    """
    Transform ITRF/ITRS positions and velocities to GCRS (J2000-like inertial) using Astropy.
    Returns (pos_j2000_km, vel_j2000_km_s).
    """
    obstimes = Time(epochs_dt)
    rep = coord.CartesianRepresentation(
        x=pos_km[:, 0] * u.km,
        y=pos_km[:, 1] * u.km,
        z=pos_km[:, 2] * u.km
    )
    diff = coord.CartesianDifferential(
        d_x=vel_km_s[:, 0] * u.km / u.s,
        d_y=vel_km_s[:, 1] * u.km / u.s,
        d_z=vel_km_s[:, 2] * u.km / u.s
    )
    itrs = coord.ITRS(rep.with_differentials(diff), obstime=obstimes)
    gcrs = itrs.transform_to(coord.GCRS(obstime=obstimes))
    p = np.column_stack([gcrs.cartesian.x.to_value(u.km),
                         gcrs.cartesian.y.to_value(u.km),
                         gcrs.cartesian.z.to_value(u.km)])
    v = np.column_stack([gcrs.velocity.d_x.to_value(u.km/u.s),
                         gcrs.velocity.d_y.to_value(u.km/u.s),
                         gcrs.velocity.d_z.to_value(u.km/u.s)])
    return p, v

def write_positions_csv_ecef_from_epochs(csv_filename, epochs_dt, pos_ecef_km):
    """
    Write ECEF positions (km) with datetime epochs to a CSV in the same format as the TLE simulation output.
    Columns: Date, X, Y, Z.
    """
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Date", "X", "Y", "Z"])
        for dt, p in zip(epochs_dt, pos_ecef_km):
            date_str = dt.strftime("%Y %m %d %H %M %S") + ".00"
            writer.writerow([date_str,
                             f"{p[0]:16.6f}",
                             f"{p[1]:16.6f}",
                             f"{p[2]:16.6f}"])
    print(f"CSV file '{csv_filename}' with ECEF positions created successfully.")

def write_sp3_spice_kernel(sp3_file, sat_id, naif_id, kernel_filename,
                           step_seconds=1, unit_is_km=True,
                           center_id=399, frame="J2000", segid=None):
    """
    Generate an SPK (type 8) kernel from an SP3 precise ephemeris file for a single satellite.
    The SP3 is assumed to provide positions in ITRF/ECEF; velocities are estimated by finite differences.
    """
    load_spice_kernels()

    epochs, pos_ecef_km = parse_sp3_positions(sp3_file, sat_id=sat_id, unit_is_km=unit_is_km)

    if step_seconds is not None and step_seconds > 0:
        epochs, pos_ecef_km = resample_positions(epochs, pos_ecef_km, step_seconds=step_seconds)

    vel_ecef_km_s = finite_difference_velocities(pos_ecef_km, step_seconds=step_seconds if step_seconds else 1)
    p_j2000, v_j2000 = itrf_to_gcrs(epochs, pos_ecef_km, vel_ecef_km_s)

    et_times = np.array([spice.utc2et(e.strftime("%Y-%m-%dT%H:%M:%S")) for e in epochs], dtype=float)
    states_matrix = np.hstack((p_j2000, v_j2000)).tolist()

    if segid is None:
        segid = f"SP3_SPK_{sat_id}"

    handle = spice.spkopn(kernel_filename, f"SPK from SP3 for {sat_id}", 0)
    spice.spkw08(handle, int(naif_id), int(center_id), frame,
                 float(et_times[0]), float(et_times[-1]),
                 segid, 7, len(et_times), states_matrix,
                 float(et_times[0]), float(step_seconds if step_seconds else 1))
    spice.spkcls(handle)
    print(f"SPK kernel '{kernel_filename}' created successfully from SP3 for sat {sat_id} (NAIF ID {naif_id}).")

def run_simulation(input_file, csv_output=None, output_folder="output",
                   source="auto", sp3_sat_id=None, sp3_naif_id=None,
                   sp3_step_seconds=1, sp3_unit_is_km=True,
                   write_csv=False):
    """
    Run either a TLE-based simulation or an SP3-based kernel build, depending on the input or 'source'.

    Parameters
    ----------
    input_file : str
        Path to a TLE file (.txt) or an SP3 file (.sp3).
    csv_output : str or None
        Optional output CSV path. If None, a default name is created in output_folder.
    output_folder : str
        Folder where outputs (SPK/CSV) will be written.
    source : {'auto','tle','sp3'}
        Select processing mode. 'auto' infers from file extension.
    sp3_sat_id : str or None
        Satellite identifier as it appears in the SP3 (e.g., 'PG01', 'PL99'). Required for SP3 mode.
    sp3_naif_id : int or None
        Integer NAIF ID to assign in the SPK for this satellite. Required for SP3 mode.
    sp3_step_seconds : int
        Resampling step for SP3 positions before kernel generation.
    sp3_unit_is_km : bool
        Whether SP3 XYZ units are kilometers (True) or meters (False).
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if source not in ("auto", "tle", "sp3"):
        raise ValueError("source must be one of: 'auto', 'tle', 'sp3'.")

    ext = os.path.splitext(input_file)[1].lower()
    mode = source
    if source == "auto":
        if ext == ".sp3":
            mode = "sp3"
        else:
            mode = "tle"

    base_name = os.path.splitext(os.path.basename(input_file))[0]
    kernel_filename = os.path.join(output_folder, base_name + ".bsp")
    if write_csv and (not csv_output):
        csv_output = os.path.join(output_folder, base_name + "_ecef.csv")

    if mode == "tle":
        print("Loading TLE records...")
        tle_list = load_tle_records(input_file)

        print("Starting simulation of satellite states (TEME)...")
        times, positions, velocities = simulate_tle_states(tle_list, timestep_seconds=1)
        times, positions, velocities = remove_duplicate_epochs(times, positions, velocities)

        print("Writing J2000 SPICE kernel...(be patient!)")
        write_j2000_spice_kernel(kernel_filename, times, positions, velocities, input_file)

        if write_csv:
            print("Writing ECEF positions to CSV for comparison...")
            write_positions_csv_ecef(csv_output, times, positions, velocities)
        print("Done.")
        return

    if mode == "sp3":
        if sp3_sat_id is None or sp3_naif_id is None:
            raise ValueError("SP3 mode requires sp3_sat_id and sp3_naif_id.")

        print("Parsing SP3 and building SPK...")
        write_sp3_spice_kernel(
            sp3_file=input_file,
            sat_id=sp3_sat_id,
            naif_id=sp3_naif_id,
            kernel_filename=kernel_filename,
            step_seconds=sp3_step_seconds,
            unit_is_km=sp3_unit_is_km
        )

        if write_csv:
            print("Writing ECEF positions to CSV for comparison...")
            epochs, pos_ecef_km = parse_sp3_positions(input_file, sat_id=sp3_sat_id, unit_is_km=sp3_unit_is_km)
            if sp3_step_seconds is not None and sp3_step_seconds > 0:
                epochs, pos_ecef_km = resample_positions(epochs, pos_ecef_km, step_seconds=sp3_step_seconds)
            write_positions_csv_ecef_from_epochs(csv_output, epochs, pos_ecef_km)
        else:
            print("Skipping CSV generation.")
        print("Done.")
        return


if __name__ == '__main__':
    run_simulation(
        r"C:\Users\User\Downloads\GRG0OPSULT_20260341800_02D_05M_ORB.SP3\GRG0OPSULT_20260341800_02D_05M_ORB.SP3",
        source="sp3",
        sp3_sat_id="PC19",
        sp3_naif_id=7000019,
        sp3_step_seconds=300,
        sp3_unit_is_km=True,
        output_folder="output"
    )
