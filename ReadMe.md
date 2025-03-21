# Satellite Propagation and Analysis

This project provides tools to simulate satellite propagation using Two-Line Element (TLE) data, generate SPICE kernel files, convert state vectors between coordinate systems (TEME, ECEF, J2000/GCRS), and perform occultation analysis. It also includes functionality to compare simulated positions with those derived from SP3 files.

## Features

- **TLE Simulation**  
  Propagate satellite states using the [sgp4](https://pypi.org/project/sgp4/) library.

- **Coordinate Conversion**  
  Convert TEME state vectors to:
  - **ECEF** coordinates using [Skyfield](https://rhodesmill.org/skyfield/)  
  - **J2000** coordinates using [Astropy](https://www.astropy.org/)

- **SPICE Kernel Generation**  
  Generate SPK kernel files using [spiceypy](https://pypi.org/project/spiceypy/).

- **Occultation Analysis**  
  Analyze occultation events between satellites (module provided).

- **Position Comparison**  
  Compare simulated positions (from TLE propagation) with positions extracted from SP3 files.

## Project Structure

repo/
├── kernels/
│   ├── naif0012.tls         # SPICE kernel file
│   ├── pck00011.tpc         # SPICE kernel file
│   └── de432s.bsp           # SPICE kernel file
├── examples/
│   ├── SP3/
│   │   └── 44351_COSMIC2-2_2022.sp3   # Example SP3 file
│   └── TLE/
│       └── 46317_LEMUR2SQUAREJAWS.txt  # Example TLE file
├── output/                 # Folder for generated output (SPK kernels, CSV files)
├── src/
│   ├── __init__.py         
│   ├── tle_simulation.py   # Module for TLE propagation, coordinate conversion, and SPICE kernel generation
│   ├── occultation.py      # Module for occultation analysis
│   └── compare_positions.py# Module for comparing SP3 and simulated positions
├── main.py                 # Main entry point – choose mode (simulation, occultation, compare)
└── README.md               # Project documentation

## Dependencies

The project requires Python 3.x and the following packages:

- numpy
- sgp4
- skyfield
- astropy
- spiceypy

## Usage

The project is designed to be run via main.py. This file allows you to select one of three modes:

- **simulation:** Run TLE simulation to generate a SPICE kernel and a CSV file with ECEF positions.
- **occultation:** Run occultation analysis using satellite data from specified folders.
- **compare:** Compare positions from SP3 files with simulated positions from the CSV file.

## Funding and Acknowledgments

This work was funded by the National Science Centre (NCN) Poland, by the Grant: UMO-2020/37/B/ST10/03703.
We thank SPIRE company for delivering RO data.

## Acknowledgments

This project makes use of:

The sgp4 library for TLE propagation.
The Skyfield library for coordinate transformations.
Astropy for time and coordinate management.
spiceypy for SPICE kernel operations.
Contributions and suggestions are welcome!
