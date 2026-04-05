"""
Constellation Simulation Module

This module processes TLE and SP3 satellite files for LEO and GNSS constellations.
It converts satellite positions from TLE/SP3 formats to TEME/ECEF/ENU coordinates,
computes azimuth/elevation from observer points, and calculates DOP values.

Dependencies:
  - os, ast, csv
  - numpy, pandas
  - sgp4, skyfield
  - pathlib
  - datetime
"""

import os
import ast
import csv
import sys
import time
import numpy as np
import pandas as pd
from sgp4.api import Satrec
from skyfield.api import load
from datetime import datetime, timedelta, timezone


# ------------------------------
# 1. TLE processing
# ------------------------------
def process_subconstellation(simulation_start_date, interval_between_epochs, number_of_epochs, tle_folder, output_file):
    """
    Processes TLE files in a folder and generates satellite positions in TEME coordinates.
    
    Args:
        YEAR, MONTH, DAY, HOUR, MINUTE: start epoch
        INTERVAL_IN_SECONDS: time step between epochs
        NUMBER_OF_EPOCHS: number of epochs to compute
        tle_folder: folder containing TLE files
        output_file: path to save resulting CSV
    Returns:
        DataFrame with satellite positions in TEME frame
    """
    start_epoch = datetime.strptime(simulation_start_date, "%Y-%m-%d %H:%M:%S")
    epochs = [start_epoch + timedelta(seconds=interval_between_epochs * i) for i in range(number_of_epochs + 1)]
    ts = load.timescale()
    data = []

    tle_files = [f for f in os.listdir(tle_folder) if f.endswith(".txt")]

    for tle_file_name in tle_files:
        sat_name = tle_file_name.split("_")[1].split(".")[0]
        sat_id = tle_file_name.split("_")[0]

        with open(os.path.join(tle_folder, tle_file_name), "r") as tle_file:
            line0 = tle_file.readline().strip()
            line1 = tle_file.readline().strip()
            line2 = tle_file.readline().strip()

        satellite = Satrec.twoline2rv(line1, line2)

        for epoch in epochs:
            t = ts.utc(epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute, epoch.second)
            jd, fr = t.ut1, t.ut1_fraction
            e, r, v = satellite.sgp4(jd, fr)
            if e == 0:
                data.append({
                    "Constellation": sat_name,
                    "Epoch": epoch,
                    "Satellite_ID": sat_id,
                    "X": r[0],
                    "Y": r[1],
                    "Z": r[2],
                    "vX": v[0],
                    "vY": v[1],
                    "vZ": v[2]
                })
            else:
                print(f"SGP4 error {e} for satellite {sat_id} at epoch {epoch}")

    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    return df

# ------------------------------
# 2. SP3 processing
# ------------------------------
def read_sp3_file(filepath, prn, epoch=None):
    """
    Reads SP3 file for a given satellite PRN and optional epoch filter.
    Returns a DataFrame with X, Y, Z positions.
    """
    data = []
    current_epoch = None
    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('*'):
                parts = line[1:].split()
                current_epoch = f"{int(parts[0]):04d} {int(parts[1]):02d} {int(parts[2]):02d} {int(parts[3]):02d} {int(parts[4]):02d} {int(float(parts[5])):02d}"
            elif line.startswith(f'P{prn}'):
                if (epoch is None) or (current_epoch == epoch):
                    parts = line.split()
                    if len(parts) >= 4:
                        x, y, z = map(float, parts[1:4])
                        data.append([current_epoch, x, y, z])
    if data:
        return pd.DataFrame(data, columns=["Epoch", "X", "Y", "Z"])
    else:
        return None

def process_sp3_folder(simulation_start_date, interval_between_epochs, number_of_epochs, sp3_folder, output_file):
    """
    Processes all SP3 files in a folder, selects satellites and epochs, and saves to CSV.
    """
    start_timer = time.time()

    start_epoch = datetime.strptime(simulation_start_date, "%Y-%m-%d %H:%M:%S")
    epochs = [start_epoch + timedelta(seconds=interval_between_epochs * i) for i in range(number_of_epochs)]

    systems = ["G", "R", "E", "C"]
    prn_range = range(1, 61)

    all_data = pd.DataFrame()
    sp3_files = [f for f in os.listdir(sp3_folder) if f.lower().endswith(".sp3")]

    total_iterations = len(sp3_files) * len(epochs) * len(systems) * len(prn_range)
    current_iteration = 0

    for sp3_file in sp3_files:
        filepath = os.path.join(sp3_folder, sp3_file)

        for epoch in epochs:
            epoch_str = f"{epoch.year:04d} {epoch.month:02d} {epoch.day:02d} {epoch.hour:02d} {epoch.minute:02d} {epoch.second:02d}"

            for system in systems:
                for number in prn_range:

                    current_iteration += 1

                    progress = current_iteration / total_iterations
                    bar_length = 40
                    filled = int(bar_length * progress)
                    bar = "█" * filled + "-" * (bar_length - filled)

                    elapsed = time.time() - start_timer
                    eta = (elapsed / progress - elapsed) if progress > 0 else 0

                    eta_hours = int(eta // 3600)
                    eta_minutes = int((eta % 3600) // 60)
                    eta_seconds = eta % 60
                    
                    eta_formatted = f"{eta_hours:02d}:{eta_minutes:02d}:{eta_seconds:04.1f}"

                    sys.stdout.write(
                        f"\r|{bar}| {progress*100:6.2f}% "
                        f"{current_iteration}/{total_iterations} "
                        f"ETA: {eta_formatted}"
                    )
                    sys.stdout.flush()

                    prn = f"{system}{number:02d}"
                    records = read_sp3_file(filepath, prn, epoch_str)

                    if records is not None and not records.empty:
                        records["Satellite_ID"] = prn
                        cols = records.columns.tolist()
                        cols = [cols[0], "Satellite_ID"] + cols[1:-1]
                        records = records[cols]
                        all_data = pd.concat([all_data, records], ignore_index=True)

    print()

    if not all_data.empty:
        all_data["vX"], all_data["vY"], all_data["vZ"] = 0, 0, 0
        all_data["Epoch"] = pd.to_datetime(
            all_data["Epoch"],
            errors='coerce',
            format="%Y %m %d %H %M %S"
        )
        all_data.to_csv(output_file, index=False)
    else:
        print("No SP3 data to save. Check folder and epochs.")

    end_timer = time.time()
    print(f"\nTotal execution time: {end_timer - start_timer:.2f} seconds\n")

# ------------------------------
# 3. Interpolation
# ------------------------------
def interpolate_sp3(df_file, simulation_start_date, interval_between_epochs, number_of_epochs, output_file):
    """
    Interpolates SP3 satellite positions to match TLE epochs.
    """
    df = pd.read_csv(df_file)
    df["Epoch"] = pd.to_datetime(df["Epoch"])
    start_time = datetime.strptime(simulation_start_date, "%Y-%m-%d %H:%M:%S")
    end_time = start_time + timedelta(seconds=interval_between_epochs * (number_of_epochs - 1))
    df = df[df["Epoch"].between(start_time, end_time)]

    def resample_sat(group):
        sat_id = group["Satellite_ID"].iloc[0]
    
        group = group[(group["Epoch"] >= start_time) & (group["Epoch"] <= end_time)]
    
        if group.empty:
            return pd.DataFrame()
    
        group = group.sort_values("Epoch").drop_duplicates(subset="Epoch").set_index("Epoch").asfreq(f"{interval_between_epochs}s")
    
        numeric_cols = group.select_dtypes(include=[np.number]).columns
        group[numeric_cols] = group[numeric_cols].interpolate()
    
        group["Satellite_ID"] = sat_id
    
        return group.reset_index()

    df_resampled = df.groupby("Satellite_ID", group_keys=False).apply(resample_sat).reset_index(drop=True)
    df_resampled.to_csv(output_file, index=False)
    return df_resampled

# ------------------------------
# 4. TEME to ECEF
# ------------------------------
def julian_date(dt):
    """Computes Julian Date from datetime."""
    year, month = dt.year, dt.month
    day = dt.day + dt.hour / 24 + dt.minute / 1440 + dt.second / 86400
    if month <= 2:
        year -= 1
        month += 12
    A = int(year / 100)
    B = 2 - A + int(A / 4)
    return int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5

def gmst(dt):
    """Computes Greenwich Mean Sidereal Time in radians."""
    jd = julian_date(dt)
    T = (jd - 2451545.0) / 36525
    gmst_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * T**2 - T**3 / 38710000
    return np.deg2rad(gmst_deg % 360)

def teme_to_ecef(r_teme, v_teme, dt_utc):
    """Converts TEME position/velocity to ECEF frame."""
    theta = gmst(dt_utc)
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    r_ecef = rot @ np.array(r_teme)
    omega_earth = 7.2921150e-5
    v_ecef = rot @ v_teme + np.cross([0,0,omega_earth], r_ecef)
    return r_ecef, v_ecef

def convert_tle_to_ecef(tle_file, output_file):
    """Converts TLE TEME coordinates to ECEF and saves CSV."""
    tle_to_teme = pd.read_csv(tle_file)
    with open(output_file, mode='w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Constellation", "Epoch", "Satellite_ID", "X", "Y", "Z", "vX", "vY", "vZ"])
        for (constellation, epoch), group in tle_to_teme.groupby(["Constellation", "Epoch"]):
            dt_utc = pd.to_datetime(epoch).replace(tzinfo=timezone.utc)
            for sat in group.itertuples():
                r_ecef, v_ecef = teme_to_ecef([sat.X, sat.Y, sat.Z], [sat.vX, sat.vY, sat.vZ], dt_utc)
                writer.writerow([constellation, epoch, sat.Satellite_ID, *r_ecef, *v_ecef])

# ------------------------------
# 5. ECEF <-> ENU
# ------------------------------
a = 6378137.0
f = 1/298.257223563
e2 = f*(2-f)

def geodetic_to_ecef(lat, lon, h):
    """Converts geodetic coordinates to ECEF."""
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    N = a / np.sqrt(1 - e2 * np.sin(lat)**2)
    return np.array([
        (N + h) * np.cos(lat) * np.cos(lon),
        (N + h) * np.cos(lat) * np.sin(lon),
        (N*(1-e2)+h) * np.sin(lat)
    ])

def ecef_to_enu(r_ecef, lat, lon, h):
    """Converts ECEF coordinates to local ENU frame."""
    obs_ecef = geodetic_to_ecef(lat, lon, h)
    dx = r_ecef - obs_ecef
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    R = np.array([
        [-np.sin(lon), np.cos(lon), 0],
        [-np.sin(lat)*np.cos(lon), -np.sin(lat)*np.sin(lon), np.cos(lat)],
        [np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]
    ])
    return R @ dx

def enu_files(points_file, ecef_file, output_file, constellation_name_col="Constellation"):
    """
    Computes ENU coordinates from ECEF for a list of observer points.
    
    If constellation_name_col is None, defaults to 'GNSS'.
    """
    points = pd.read_csv(points_file)
    satellites = pd.read_csv(ecef_file)
    with open(output_file, mode='w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["Reference", "Constellation", "Epoch", "Satellite_ID", "E", "N", "U"])
        for point in points.itertuples():
            for sat in satellites.itertuples():
                r_enu = ecef_to_enu([sat.X*1000, sat.Y*1000, sat.Z*1000], point.Latitude, point.Longitude, point.Altitude)
                if constellation_name_col is None:
                    const_name = "GNSS"
                else:
                    const_name = getattr(sat, constellation_name_col, "GNSS")
                
                writer.writerow([point.Name, const_name, sat.Epoch, sat.Satellite_ID, *r_enu])

# ------------------------------
# 6. ENU -> AZ/EL
# ------------------------------
def az_el_from_enu(enu):
    """Converts ENU vector to azimuth and elevation angles."""
    e, n, u = enu
    az = np.degrees(np.arctan2(e, n)) % 360
    el = np.degrees(np.arcsin(u / np.linalg.norm(enu)))
    return az, el

def enu_to_az_el(enu_file, output_file, constellation_name_col="Constellation"):
    """Converts ENU coordinates CSV to azimuth/elevation CSV."""
    df_enu = pd.read_csv(enu_file)
    with open(output_file, mode='w', newline='') as f_out:
        writer = csv.writer(f_out)
        headers = ["Reference", "Constellation", "Epoch", "Satellite_ID", "Azimuth", "Elevation"] if constellation_name_col=="Constellation" else ["Reference", "Epoch", "Satellite_ID", "Azimuth", "Elevation"]
        writer.writerow(headers)
        for sat in df_enu.itertuples():
            az, el = az_el_from_enu([sat.E, sat.N, sat.U])
            row = [sat.Reference, getattr(sat, constellation_name_col, "GNSS"), sat.Epoch, sat.Satellite_ID, az, el] if constellation_name_col=="Constellation" else [sat.Reference, sat.Epoch, sat.Satellite_ID, az, el]
            writer.writerow(row)

# ------------------------------
# 7. DOP calculation
# ------------------------------
def compute_dop(sat_positions_enu):
    """Computes GDOP, PDOP, HDOP, VDOP, TDOP from ENU satellite positions."""
    H = [[e/np.linalg.norm(enu), n/np.linalg.norm(enu), u/np.linalg.norm(enu), 1] for enu in sat_positions_enu for e, n, u in [enu]]
    H = np.array(H)
    Q = np.linalg.inv(H.T @ H)
    return np.sqrt(Q[0,0]+Q[1,1]+Q[2,2]+Q[3,3]), np.sqrt(Q[0,0]+Q[1,1]+Q[2,2]), np.sqrt(Q[0,0]+Q[1,1]), np.sqrt(Q[2,2]), np.sqrt(Q[3,3])

def run_constellation_simulation(observer_position, simulation_start_date, interval_between_epochs, number_of_epochs, dop_results_output, other_results_output):
    """Runs full constellation simulation pipeline and saves outputs to CSV files."""
    # Input folders
    tle_folder = "./examples/constellation_simulation/TLE"
    sp3_folder = "./examples/constellation_simulation/SP3"
    
    # Output folder
    output_folder = "./examples/constellation_simulation/output"
    
    # Output files
    tle_teme_file = f"{other_results_output}/TLE_to_TEME.txt"
    sp3_ecef_file = f"{other_results_output}/SP3_to_ECEF.txt"
    sp3_ecef_int_file = f"{other_results_output}/SP3_to_ECEF_int.txt"
    tle_ecef_file = f"{other_results_output}/TLE_TEME_to_ECEF.txt"
    tle_enu_file = f"{other_results_output}/TLE_ECEF_to_ENU.txt"
    sp3_enu_file = f"{other_results_output}/SP3_ECEF_to_ENU.txt"
    tle_az_el_file = f"{other_results_output}/TLE_ENU_to_AZ_EL.txt"
    sp3_az_el_file = f"{other_results_output}/SP3_ENU_to_AZ_EL.txt"
    dop_file = f"{dop_results_output}/DOP.txt"
    
    if not os.path.exists(dop_results_output): os.makedirs(dop_results_output)
    if not os.path.exists(other_results_output): os.makedirs(other_results_output)
    
    # 1. TLE -> TEME
    print("Processing TLE files...")
    process_subconstellation(simulation_start_date, interval_between_epochs, number_of_epochs,
                             tle_folder, tle_teme_file)
    
    # 2. SP3 -> ECEF
    print("Processing SP3 files...\n")
    process_sp3_folder(simulation_start_date, interval_between_epochs, number_of_epochs,
                       sp3_folder, sp3_ecef_file)
    
    # 3. SP3 interpolation to TLE epochs
    print("Interpolating SP3 to TLE epochs...")
    interpolate_sp3(sp3_ecef_file, simulation_start_date, interval_between_epochs, number_of_epochs,
                    sp3_ecef_int_file)
    
    # 4. TEME (TLE) -> ECEF
    print("Converting TLE TEME -> ECEF...")
    convert_tle_to_ecef(tle_teme_file, tle_ecef_file)
    
    # 5. ECEF (TLE and SP3) -> ENU
    print("Converting ECEF -> ENU...")
    enu_files(observer_position, tle_ecef_file, tle_enu_file, constellation_name_col="Constellation")
    enu_files(observer_position, sp3_ecef_int_file, sp3_enu_file, constellation_name_col=None)
    
    # 6. ENU -> Azimuth/Elevation
    print("Converting ENU -> Azimuth/Elevation...")
    enu_to_az_el(tle_enu_file, tle_az_el_file, constellation_name_col="Constellation")
    enu_to_az_el(sp3_enu_file, sp3_az_el_file, constellation_name_col=None)
    
    # 7. Computing DOP
    print("Computing DOP values...")
    points = pd.read_csv(observer_position)
    tle_enu = pd.read_csv(tle_enu_file)
    sp3_enu = pd.read_csv(sp3_enu_file)
    tle_az_el = pd.read_csv(tle_az_el_file)
    sp3_az_el = pd.read_csv(sp3_az_el_file)
    
    tle_enu_az_el = tle_enu.merge(tle_az_el, how="inner")
    sp3_enu_az_el = sp3_enu.merge(sp3_az_el, how="inner")
    
    sp3_enu_az_el["Constellation"] = "GNSS"
    
    with open(dop_file, mode='w', newline='') as f_out:
        writer = csv.writer(f_out)
        writer.writerow([
            "Reference", "Epoch", "Constellation",
            "LEO Satellites", "GNSS Satellites",
            "GDOP", "PDOP", "HDOP", "VDOP", "TDOP",
            "LEO Satellite IDs", "GNSS Satellite IDs",
            "GPS Satellites", "GLONASS Satellites", "Galileo Satellites", "BeiDou Satellites", "QZSS Satellites", "NavIC Satellites"
        ])
    
    with open(dop_file, mode='a', newline='') as f_out:
        writer = csv.writer(f_out)
        
        for (epoch, constellation, reference), group_leo in tle_enu_az_el.groupby(["Epoch", "Constellation", "Reference"]):
            sats_enu_leo, sats_enu_gnss = [], []
            sats_leo_ids, sats_gnss_ids = [], []

            match = points.loc[points["Name"] == reference, "Elevation"]
            if match.empty:
                continue
            elevation_value = match.iloc[0]
            if isinstance(elevation_value, str):
                elevation_value = ast.literal_eval(elevation_value)

            group_gnss = sp3_enu_az_el[(sp3_enu_az_el["Epoch"] == epoch) & (sp3_enu_az_el["Reference"] == reference)]
            
            # LEO satellites selection
            for sat in group_leo.itertuples():
                for az1, az2, el_thresh in zip(range(0, 315+1, 45), range(45, 360+1, 45), elevation_value):
                    if az1 <= sat.Azimuth < az2 and sat.Elevation > el_thresh:
                        sats_enu_leo.append([sat.E, sat.N, sat.U])
                        sats_leo_ids.append(str(sat.Satellite_ID))
            
            # GNSS satellites selection
            for sat in group_gnss.itertuples():
                for az1, az2, el_thresh in zip(range(0, 315+1, 45), range(45, 360+1, 45), elevation_value):
                    if az1 <= sat.Azimuth < az2 and sat.Elevation > el_thresh:
                        sats_enu_gnss.append([sat.E, sat.N, sat.U])
                        sats_gnss_ids.append(str(sat.Satellite_ID))
            
            num_leo = len(sats_enu_leo)
            num_gnss = len(sats_enu_gnss)
            
            if num_leo + num_gnss >= 4:
                gdop, pdop, hdop, vdop, tdop = compute_dop(sats_enu_leo + sats_enu_gnss)
            else:
                gdop = pdop = hdop = vdop = tdop = np.nan
            
            writer.writerow([
                reference, epoch, constellation,
                num_leo, num_gnss,
                gdop, pdop, hdop, vdop, tdop,
                ', '.join(sats_leo_ids),
                ', '.join(sats_gnss_ids),
                sum(1 for sat_id in sats_gnss_ids if sat_id.startswith('G')),
                sum(1 for sat_id in sats_gnss_ids if sat_id.startswith('R')),
                sum(1 for sat_id in sats_gnss_ids if sat_id.startswith('E')),
                sum(1 for sat_id in sats_gnss_ids if sat_id.startswith('C')),
                sum(1 for sat_id in sats_gnss_ids if sat_id.startswith('J')),
                sum(1 for sat_id in sats_gnss_ids if sat_id.startswith('I'))
            ])
    
    print("Constellation simulation finished. All files saved in:", output_folder)
