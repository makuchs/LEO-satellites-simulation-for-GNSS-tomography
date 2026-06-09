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
from project_paths import DE432_FILE, EARTH_BPC_FILE, LSK_FILE, PCK_TPC_FILE

def load_spice_kernels():
    """
    Loads the required SPICE kernels from the repository's kernels folder.
    """
    spice.kclear()
    spice.furnsh(str(LSK_FILE))
    spice.furnsh(str(PCK_TPC_FILE))
    spice.furnsh(str(DE432_FILE))
    if EARTH_BPC_FILE.exists():
        spice.furnsh(str(EARTH_BPC_FILE))

def load_tle_records(filename):
    tle_records = []
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        raise ValueError("TLE file is empty.")

    if all((ln.startswith("1 ") or ln.startswith("2 ")) for ln in lines):
        if len(lines) % 2 != 0:
            raise ValueError("TLE file format error: 2-line TLE file has odd number of nonempty lines.")
        base_name = os.path.splitext(os.path.basename(filename))[0]
        rec_idx = 1
        for i in range(0, len(lines), 2):
            l1, l2 = lines[i], lines[i + 1]
            if not l1.startswith("1 ") or not l2.startswith("2 "):
                raise ValueError(f"TLE file format error near lines {i+1}-{i+2}.")
            tle_records.append((f"{base_name}_{rec_idx:05d}", l1, l2))
            rec_idx += 1
        return tle_records

    if len(lines) % 3 != 0:
        raise ValueError("TLE file format error: total nonempty lines not a multiple of 3.")

    for i in range(0, len(lines), 3):
        name, l1, l2 = lines[i], lines[i + 1], lines[i + 2]
        if not l1.startswith("1 ") or not l2.startswith("2 "):
            raise ValueError(f"TLE file format error near lines {i+1}-{i+3}.")
        tle_records.append((name, l1, l2))

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

def write_j2000_spice_kernel(kernel_filename, times, positions, velocities, tle_filename, step_seconds=30):
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

    spice.spkw08(
        handle,
        int(sat_norad),
        399,
        "J2000",
        et_times[0],
        et_times[-1],
        segid,
        7,
        len(et_times),
        states_matrix,
        et_times[0],
        float(step_seconds),
    )

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

def list_sp3_satellite_ids(sp3_file):
    """
    Return sorted unique satellite ID strings from SP3 position records (first column
    of each XYZ line under a valid ``*`` epoch). Uses the same epoch/line rules as
    ``parse_sp3_positions`` so the listed IDs are those that can be kernelized.
    """
    seen = set()
    current_epoch = None

    with open(sp3_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith("*"):
                tokens = line[1:].strip().split()
                if len(tokens) < 6:
                    current_epoch = None
                    continue
                year = int(tokens[0])
                month = int(tokens[1])
                day = int(tokens[2])
                hour = int(tokens[3])
                minute = int(tokens[4])
                second = float(tokens[5])
                sec_int = int(second)
                micro = int(round((second - sec_int) * 1e6))
                current_epoch = datetime.datetime(
                    year, month, day, hour, minute, sec_int, micro
                )
                continue

            if current_epoch is None:
                continue

            tokens = line.split()
            if len(tokens) < 4:
                continue
            try:
                float(tokens[1])
                float(tokens[2])
                float(tokens[3])
            except ValueError:
                continue
            seen.add(tokens[0])

    if not seen:
        raise ValueError(f"No satellite position records found in {sp3_file!r}.")
    return sorted(seen)

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
                           center_id=399, frame="J2000", segid=None,
                           epochs_dt=None, pos_ecef_km=None):
    """
    Generate an SPK (type 8) kernel from an SP3 precise ephemeris file for a single satellite.
    The SP3 is assumed to provide positions in ITRF/ECEF; velocities are estimated by finite differences.

    If epochs_dt and pos_ecef_km are provided, sp3_file is not read (same arrays as from parse_sp3_positions).

    Returns (epochs_dt, pos_ecef_km) after any resampling — the same samples written into the SPK.
    """
    load_spice_kernels()

    if (epochs_dt is None) ^ (pos_ecef_km is None):
        raise ValueError("epochs_dt and pos_ecef_km must both be set or both omitted.")
    if epochs_dt is not None:
        epochs, pos_ecef_km = epochs_dt, pos_ecef_km
    else:
        epochs, pos_ecef_km = parse_sp3_positions(sp3_file, sat_id=sat_id, unit_is_km=unit_is_km)

    if step_seconds is not None and step_seconds > 0:
        epochs, pos_ecef_km = resample_positions(epochs, pos_ecef_km, step_seconds=step_seconds)

    vel_ecef_km_s = finite_difference_velocities(pos_ecef_km, step_seconds=step_seconds if step_seconds else 1)
    p_j2000, v_j2000 = itrf_to_gcrs(epochs, pos_ecef_km, vel_ecef_km_s)

    def _epoch_to_utc_str(e):
        if e.microsecond:
            return e.strftime("%Y-%m-%dT%H:%M:%S.%f")
        return e.strftime("%Y-%m-%dT%H:%M:%S")

    et_raw = np.array([spice.utc2et(_epoch_to_utc_str(e)) for e in epochs], dtype=np.float64)
    n_et = len(et_raw)
    if n_et >= 2:
        # SPKW08 requires a uniform ET grid; per-epoch utc2et floats can drift and trigger COVERAGEGAP.
        et0 = float(et_raw[0])
        step_et = float(et_raw[1] - et_raw[0])
        et_times = et0 + np.arange(n_et, dtype=np.float64) * step_et
    else:
        et_times = et_raw.astype(np.float64)

    states_matrix = np.hstack((p_j2000, v_j2000)).tolist()

    if segid is None:
        segid = f"SP3_SPK_{sat_id}"

    if len(et_times) > 1:
        step_for_spk = float(et_times[1] - et_times[0])
    else:
        step_for_spk = float(step_seconds if step_seconds else 1)

    handle = spice.spkopn(kernel_filename, f"SPK from SP3 for {sat_id}", 0)
    spice.spkw08(handle, int(naif_id), int(center_id), frame,
                 float(et_times[0]), float(et_times[-1]),
                 segid, 7, len(et_times), states_matrix,
                 float(et_times[0]), step_for_spk)
    spice.spkcls(handle)
    print(f"SPK kernel '{kernel_filename}' created successfully from SP3 for sat {sat_id} (NAIF ID {naif_id}).")
    return epochs, pos_ecef_km

def run_simulation(input_file, csv_output=None, output_folder="output",
                   source="auto", sp3_sat_id=None, sp3_naif_id=None,
                   sp3_step_seconds=1, sp3_unit_is_km=True,
                   write_csv=False, tle_step_seconds=60):
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
    if write_csv and csv_output:
        os.makedirs(os.path.dirname(csv_output) or ".", exist_ok=True)

    if mode == "tle":
        print("Loading TLE records...")
        tle_list = load_tle_records(input_file)

        print("Starting simulation of satellite states (TEME)...")
        times, positions, velocities = simulate_tle_states(
            tle_list,
            timestep_seconds=int(tle_step_seconds)
        )
        times, positions, velocities = remove_duplicate_epochs(times, positions, velocities)

        print("Writing J2000 SPICE kernel...(be patient!)")
        write_j2000_spice_kernel(
            kernel_filename,
            times,
            positions,
            velocities,
            input_file,
            step_seconds=int(tle_step_seconds),
        )

        if write_csv:
            print("Writing ECEF positions to CSV for comparison...")
            write_positions_csv_ecef(csv_output, times, positions, velocities)
        print("Done.")
        return kernel_filename

    if mode == "sp3":
        if sp3_sat_id is None or sp3_naif_id is None:
            raise ValueError("SP3 mode requires sp3_sat_id and sp3_naif_id.")

        print("Parsing SP3 and building SPK...")
        epochs, pos_ecef_km = write_sp3_spice_kernel(
            sp3_file=input_file,
            sat_id=sp3_sat_id,
            naif_id=sp3_naif_id,
            kernel_filename=kernel_filename,
            step_seconds=sp3_step_seconds,
            unit_is_km=sp3_unit_is_km
        )

        if write_csv:
            print("Writing ECEF positions to CSV for comparison...")
            write_positions_csv_ecef_from_epochs(csv_output, epochs, pos_ecef_km)
        else:
            print("Skipping CSV generation.")
        print("Done.")
        return kernel_filename


if __name__ == '__main__':
    raise SystemExit("Use main.py to run simulations, or import run_simulation() from Python.")
