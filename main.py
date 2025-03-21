"""
Main entry point for the Satellite Analysis Project.

Modes:
  - "simulation": Runs TLE simulation to generate a SPICE kernel and an ECEF CSV file.
  - "occultation": Runs occultation analysis.
  - "compare": Compares SP3 positions with positions from the simulation CSV file.

Edit the variables in the configuration section below to set file paths, dates, mode, and output folder.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

# =================== Configuration ===================
# Choose mode: "simulation", "occultation", or "compare"
MODE = "occultation"

# For simulation mode:
TLE_FILE = "examples/TLE/46317_LEMUR2SQUAREJAWS.txt"  
OUTPUT_FOLDER = "output"  # All generated files (SPK and CSV) will be stored here.

# For occultation mode:
LEO_FOLDER = "kernels/leo"         # Folder containing LEO kernels (adjust as needed)
GNSS_FOLDER = "kernels/gnss"       # Folder containing GNSS kernels (adjust as needed)
OCC_CSV_FILE = os.path.join(OUTPUT_FOLDER, "occultation_results.csv")
OCC_START_DATE = "2022-01-03"    # YYYY-MM-DD
OCC_END_DATE   = "2022-01-10"    # YYYY-MM-DD

# For compare mode:
SP3_FILE = "examples/SP3/46317_LEMUR2SQUAREJAWS_2022.sp3"
# For compare mode, the simulation CSV file will be used.
SIM_CSV_FILE_FOR_COMPARE = "output/46317_LEMUR2SQUAREJAWS_ecef.csv"  
SAT_ID = "PL99"                # Satellite identifier as it appears in the SP3 file
# =====================================================

def main():
    if MODE == "simulation":
        from src import tle_simulation
        print("Running TLE Simulation...")
        # Call run_simulation with the output folder.
        tle_simulation.run_simulation(TLE_FILE, output_folder=OUTPUT_FOLDER)
    
    elif MODE == "occultation":
        from src import occultation
        print("Running Occultation Analysis...")
        occultation.run_occultation(LEO_FOLDER, GNSS_FOLDER, OCC_CSV_FILE, OCC_START_DATE, OCC_END_DATE)
    
    elif MODE == "compare":
        from src import compare_positions
        print("Running Position Comparison...")
        compare_positions.run_compare(SP3_FILE, SIM_CSV_FILE_FOR_COMPARE, SAT_ID)
    
    else:
        print("Invalid MODE. Please set MODE to 'simulation', 'occultation', or 'compare'.")

if __name__ == '__main__':
    main()
