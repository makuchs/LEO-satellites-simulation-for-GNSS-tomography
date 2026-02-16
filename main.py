"""
Main entry point for the Satellite Analysis Project.

This file is intentionally simple: edit the CONFIGURATION section and run:

    python main.py

Available modes (set MODE below)
--------------------------------
- "simulation"
    Build a SPICE SPK kernel from a TLE or an SP3 file.
    Optionally export simulated ECEF positions to CSV.

- "occultation"
    Run radio occultation (RO) analysis (requires existing LEO and GNSS kernels).

- "gnssr"
    Run GNSS-R geometry simulation (requires existing LEO and GNSS kernels).

- "compare"
    Compare positions from an SP3 file with positions exported to CSV by "simulation".

- "constellation"
    Creates a satellite constellation from the given TLE and SP3 files and then calculates the DOP coefficients for it.

How to build kernels for GNSS-R / RO
------------------------------------
GNSS-R and occultation need two kernel folders:
- a LEO kernel folder (usually built from a TLE)
- a GNSS kernel folder (usually built from an SP3)

You create them by running MODE="simulation" twice:
1) set SIM_INPUT_FILE = LEO_TLE_FILE, SIM_SOURCE="tle", SIM_OUTPUT_FOLDER = LEO_KERNEL_FOLDER
2) set SIM_INPUT_FILE = GNSS_SP3_FILE, SIM_SOURCE="sp3", SIM_OUTPUT_FOLDER = GNSS_KERNEL_FOLDER
   and fill SIM_SP3_SAT_ID + SIM_SP3_NAIF_ID
"""

import os
import sys

# Make ./src importable when running from the repository root.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

# =================== CONFIGURATION ===================

# Select one: "simulation", "occultation", "gnssr", "compare", "constellation"
MODE = "gnssr"

# Repository-relative folders (recommended for GitHub projects)
OUTPUT_DIR = "output"
LEO_KERNEL_FOLDER = os.path.join("kernels", "leo")
GNSS_KERNEL_FOLDER = os.path.join("kernels", "gnss")

# ---------- Inputs (edit these paths) ----------
LEO_TLE_FILE = os.path.join("data", "leo.tle")          # example path
GNSS_SP3_FILE = os.path.join("data", "gnss.sp3")        # example path
GNSS_SP3_SAT_ID = "PG18"                                # as it appears in SP3 (e.g. PG18, G18)

# ---------- Simulation (MODE="simulation") ----------
# Pick which file to process (TLE or SP3) and where to write the kernel.
SIM_INPUT_FILE = LEO_TLE_FILE                           # change to GNSS_SP3_FILE for SP3 mode
SIM_SOURCE = "tle"                                      # "tle" or "sp3" (or "auto")
SIM_OUTPUT_FOLDER = LEO_KERNEL_FOLDER                   # change to GNSS_KERNEL_FOLDER for GNSS
SIM_WRITE_CSV = False                                   # export ECEF CSV (optional)
SIM_CSV_OUTPUT = None                                   # optional path; None = auto name

# SP3-only parameters (required when SIM_SOURCE="sp3")
SIM_SP3_SAT_ID = GNSS_SP3_SAT_ID
SIM_SP3_NAIF_ID = 7000018                               # integer ID used inside the SPK (choose a unique one)
SIM_SP3_STEP_SECONDS = 1                                # resampling step before kernel generation
SIM_SP3_UNIT_IS_KM = True                               # SP3 XYZ units: True=km, False=m

# ---------- GNSS-R (MODE="gnssr") ----------
GNSSR_CSV_FILE = os.path.join(OUTPUT_DIR, "gnssr_results.csv")
GNSSR_START_DATE = "2023-10-01"                         # YYYY-MM-DD
GNSSR_END_DATE = "2023-10-01"                           # YYYY-MM-DD
GNSSR_START_TIME_UTC = "00:00:00"                       # HH:MM:SS (UTC)
GNSSR_END_TIME_UTC = "23:59:59"                         # HH:MM:SS (UTC)
GNSSR_STEP_SECONDS = 10

# Optional filters (set to None to disable)
GNSSR_BBOX = None                                       # (min_lat, max_lat, min_lon, max_lon) in degrees, e.g. (49, 55, 14, 24)
GNSSR_BBOX_GUIDED_SEARCH = False
GNSSR_BBOX_GRID_DEG = 1.0
GNSSR_INITIAL_LATLON = None                             # (lat, lon) in degrees, e.g. (52.0, 21.0)

GNSSR_SPECULAR_ERR_MAX = 1.0
GNSSR_INC_MAX_DEG = None
GNSSR_BISTATIC_MAX_DEG = None
GNSSR_EXCESS_MAX_KM = None
GNSSR_LOS_TOL_KM = 0.01

# ---------- Occultation (MODE="occultation") ----------
OCC_CSV_FILE = os.path.join(OUTPUT_DIR, "occultation_results.csv")
OCC_START_DATE = "2022-01-03"
OCC_END_DATE = "2022-01-10"

# ---------- Compare (MODE="compare") ----------
CMP_SP3_FILE = GNSS_SP3_FILE
CMP_SIM_ECEF_CSV = os.path.join(OUTPUT_DIR, "simulation_ecef.csv")
CMP_SAT_ID = GNSS_SP3_SAT_ID

#----------- Constellation (MODE="constellation) -----------
OBSERVER = "examples/constellation_simulation/OBSERVER_POINTS/OBSERVER.txt"     # File with the observer's position and horizon obscuration in the format:
                                                                                # Name,Latitude,Longitude,Altitude,Elevation (elevation angles list in degrees, every 45 degrees azimuth)
SIM_START_DATE = "2025-03-31 12:00:00"                                          # YYYY-MM-DD hh:mm:ss
SIM_INTERVAL = 30                                                               # Interval between epochs, in seconds
NUM_OF_EPOCHS = 120                                                             # Number of epochs
DOP_OUTPUT = f"examples/constellation_simulation/output/DOP {SIM_START_DATE.replace(':', '-')} {SIM_INTERVAL} INT {NUM_OF_EPOCHS} EPO"
MISC_OUTPUT = f"examples/constellation_simulation/output/MISC {SIM_START_DATE.replace(':', '-')} {SIM_INTERVAL} INT {NUM_OF_EPOCHS} EPO"

# =====================================================


def _ensure_dirs() -> None:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LEO_KERNEL_FOLDER, exist_ok=True)
    os.makedirs(GNSS_KERNEL_FOLDER, exist_ok=True)


def main() -> None:
    _ensure_dirs()

    if MODE == "simulation":
        import tle_simulation

        print("Running: simulation")
        tle_simulation.run_simulation(
            SIM_INPUT_FILE,
            csv_output=SIM_CSV_OUTPUT,
            output_folder=SIM_OUTPUT_FOLDER,
            source=SIM_SOURCE,
            sp3_sat_id=SIM_SP3_SAT_ID,
            sp3_naif_id=SIM_SP3_NAIF_ID,
            sp3_step_seconds=SIM_SP3_STEP_SECONDS,
            sp3_unit_is_km=SIM_SP3_UNIT_IS_KM,
            write_csv=SIM_WRITE_CSV,
        )

    elif MODE == "occultation":
        import occultation

        print("Running: occultation")
        occultation.run_occultation(
            leo_folder=LEO_KERNEL_FOLDER,
            gnss_folder=GNSS_KERNEL_FOLDER,
            csv_file_path=OCC_CSV_FILE,
            start_date=OCC_START_DATE,
            end_date=OCC_END_DATE,
        )

    elif MODE == "gnssr":
        import gnssr

        print("Running: gnssr")
        gnssr.run_gnssr(
            leo_folder=LEO_KERNEL_FOLDER,
            gnss_folder=GNSS_KERNEL_FOLDER,
            csv_file_path=GNSSR_CSV_FILE,
            start_date=GNSSR_START_DATE,
            end_date=GNSSR_END_DATE,
            start_time_utc=GNSSR_START_TIME_UTC,
            end_time_utc=GNSSR_END_TIME_UTC,
            bbox=GNSSR_BBOX,
            step_seconds=GNSSR_STEP_SECONDS,
            specular_err_max=GNSSR_SPECULAR_ERR_MAX,
            inc_max_deg=GNSSR_INC_MAX_DEG,
            bistatic_max_deg=GNSSR_BISTATIC_MAX_DEG,
            excess_max_km=GNSSR_EXCESS_MAX_KM,
            los_tol_km=GNSSR_LOS_TOL_KM,
            bbox_guided_search=GNSSR_BBOX_GUIDED_SEARCH,
            bbox_grid_deg=GNSSR_BBOX_GRID_DEG,
            initial_latlon=GNSSR_INITIAL_LATLON,
        )

    elif MODE == "compare":
        import compare_positions

        print("Running: compare")
        compare_positions.run_compare(
            sp3_file=CMP_SP3_FILE,
            csv_file=CMP_SIM_ECEF_CSV,
            sat_id=CMP_SAT_ID,
        )

    elif MODE == "constellation":
        import constellation_simulation
        
        print("Running: constellation")
        constellation_simulation.run_constellation_simulation(
            observer_position=OBSERVER,
            simulation_start_date=SIM_START_DATE,
            interval_between_epochs=SIM_INTERVAL,
            number_of_epochs=NUM_OF_EPOCHS,
            dop_results_output=DOP_OUTPUT,
            other_results_output=MISC_OUTPUT)

    else:
        raise ValueError("Invalid MODE. Use: 'simulation', 'occultation', 'gnssr', or 'compare'.")


if __name__ == "__main__":
    main()

