# Satellite Propagation and Analysis

Tools to:
- propagate satellites from **TLE** or **SP3**,
- generate **SPICE SPK (.bsp)** kernels,
- convert state vectors between **TEME / ECEF / J2000 (GCRS)**,
- run **radio occultation (RO)** analysis,
- run basic **GNSS-R** bistatic geometry simulation,
- compare **SP3** positions with simulated **ECEF CSV** output.

---

## Requirements

- Python 3.x
- Packages: `numpy`, `sgp4`, `skyfield`, `astropy`, `spiceypy`, `requests`


---

## Project structure

```
repo/
  kernels/
    naif0012.tls
    pck00011.tpc
    de432s.bsp
    leo/        # generated LEO SPKs (.bsp) go here
    gnss/       # generated GNSS SPKs (.bsp) go here
  data/
    leo.tle     # your input files (example names)
    gnss.sp3
  output/
    *.csv       # results + optional simulated ECEF CSV
  src/
    tle_simulation.py
    occultation.py
    gnssr.py
    compare_positions.py
  misc/
    sp3_tools.py # optional: download + SP3 post-processing helpers
  main.py       # edit CONFIGURATION and run
  README.md
```

Notes:
- Base kernels under `kernels/` (LSK/PCK/planetary SPK) are required by the analysis modules.
- Generated SPKs in `kernels/leo/` and `kernels/gnss/` are not tracked in git in most setups (they can be large).

---

## Running the project

Edit the CONFIGURATION section in `main.py`, set `MODE`, then run:

```bash
python main.py
```

Available modes:

- `simulation` — build an SPK kernel from a **TLE** or **SP3** file (optional ECEF CSV export)
- `occultation` — run RO analysis (requires LEO + GNSS kernel folders)
- `gnssr` — run GNSS-R geometry simulation (requires LEO + GNSS kernel folders)
- `compare` — compare an SP3 file with an ECEF CSV produced by `simulation`
- `constellation` — build a satellite constellation from TLE/SP3 inputs and compute DOP coefficients (PDOP, GDOP, HDOP, VDOP, TDOP) for a given observer / horizon mask.


GNSS-R and occultation need two kernel folders:
- **LEO kernels** (typically built from a TLE)
- **GNSS kernels** (typically built from an SP3)

Do it in two runs of `MODE="simulation"`:

1) **Build LEO kernel from TLE**
- set `MODE = "simulation"`
- set `SIM_INPUT_FILE = LEO_TLE_FILE`
- set `SIM_SOURCE = "tle"`
- set `SIM_OUTPUT_FOLDER = LEO_KERNEL_FOLDER`
- run `python main.py`

2) **Build GNSS kernel from SP3**
- keep `MODE = "simulation"`
- set `SIM_INPUT_FILE = GNSS_SP3_FILE`
- set `SIM_SOURCE = "sp3"`
- set `SIM_OUTPUT_FOLDER = GNSS_KERNEL_FOLDER`
- set SP3-only params: `SIM_SP3_SAT_ID` and `SIM_SP3_NAIF_ID`
- run `python main.py`

Then switch `MODE` to `gnssr` or `occultation` and run again.

---

## Mode details 

### 1) `simulation`

Input:
- TLE: set `SIM_SOURCE="tle"` and point `SIM_INPUT_FILE` to a TLE file.
- SP3: set `SIM_SOURCE="sp3"` and point `SIM_INPUT_FILE` to an SP3 file, plus:
  - `SIM_SP3_SAT_ID` — satellite label as it appears in the SP3 (e.g. `PG18`)
  - `SIM_SP3_NAIF_ID` — integer NAIF ID used inside the generated SPK (pick a unique value)

Output:
- SPK `.bsp` written to `SIM_OUTPUT_FOLDER`
- optional ECEF CSV if `SIM_WRITE_CSV=True`

### 2) `gnssr`

Requires:
- `LEO_KERNEL_FOLDER` with one or more `.bsp` files
- `GNSS_KERNEL_FOLDER` with one or more `.bsp` files

Main outputs:
- results CSV (e.g. `output/gnssr_results.csv`)

Optional filters include date/time range, bbox, and basic geometry thresholds (see `main.py`).

### 3) `occultation`

Requires:
- `LEO_KERNEL_FOLDER` and `GNSS_KERNEL_FOLDER`

Output:
- results CSV (e.g. `output/occultation_results.csv`)

### 4) `compare`

Input:
- an SP3 file
- an ECEF CSV produced by `simulation`
- satellite ID used to match epochs (e.g. `PG18`)

Output:
- a plot of position differences over time

### 5) `constellation`

Input:
- observer file in the format:
  `Name,Latitude,Longitude,Altitude,Elevation`
  where `Elevation` is a list of elevation angles (degrees) sampled every 45° azimuth.
- satellite sources from SP3 and/or TLE (see `examples/constellation_simulation/`)

Output:
- DOP time series (PDOP, GDOP, HDOP, VDOP, TDOP)
- additional per-epoch constellation diagnostics (see output paths in `main.py`)

---

## Input data sources

---

## Misc utilities 

The `misc/` folder contains helper scripts that are not required to run the main pipeline,
but can speed up data preparation (downloading and cleaning SP3/LEO products).

### `misc/sp3_tools.py`

One script that combines:
- downloading GNSS SP3 products from an IGS MGEX FTP mirror (and unzipping),
- downloading LEO orbit tarballs from UCAR (and extracting),
- combining many SP3 files into one time-sorted SP3 (optionally filtering satellites),
- sorting large file dumps into subfolders.

Examples:

```bash
# GNSS SP3 (IGS MGEX): download weeks covering a date range
python misc/sp3_tools.py download-gnss --ac COD --start-date 2021-06-15 --end-date 2021-07-06 --out data/sp3/igs

# LEO orbit products (UCAR): download daily tarballs and extract them
python misc/sp3_tools.py download-leo --system spire --start-date 2023-09-30 --end-date 2023-10-03 --out data/sp3/spire

# Combine many SP3 files into one file (keep only one record ID, e.g. PG18)
python misc/sp3_tools.py combine --in data/sp3/igs --out data/PG18.sp3 --sat PG18

# Sort files into subfolders based on a filename segment (dot-split)
python misc/sp3_tools.py sort-subfolders --in data/sp3/spire --segment-index 2
```

Notes:
- `--sat` expects a SP3 record ID (e.g. `PG18`). If you want velocities too, add `VG18` as well.


### TLE
You can obtain current TLE data from:
- Space-Track 
- CelesTrak 

### SP3
Use precise ephemeris SP3 files from your GNSS data source.

---

## Kernel naming convention

If you generate multiple kernels, it is recommended to keep a consistent filename pattern such as:

```
NORAD_SATNAME.bsp
```

Example: NORAD `1234` and name `lemur2` → `1234_lemur2.bsp`

---

## Contributors

- Thanks to **@marcinwapinski** for adding a satellite constellation simulation and DOP calculation module, including:
  - building constellations from **SP3** and **TLE** inputs,
  - interpolation of satellite positions for missing epochs,
  - coordinate transformations between supported reference frames,
  - DOP metrics: **PDOP, GDOP, HDOP, VDOP, TDOP**.

---

## Funding and acknowledgments

This work was funded by the National Science Centre (NCN) Poland, Grant: **UMO-2020/37/B/ST10/03703**.  
We thank **SPIRE** for delivering RO data.

---

## Acknowledgments

This project makes use of:
- `sgp4` for TLE propagation
- Skyfield for coordinate transformations
- Astropy for time/coordinate utilities
- spiceypy for SPICE kernel operations
