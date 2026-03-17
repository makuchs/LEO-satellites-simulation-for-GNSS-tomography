"""Shared repository paths used by the analysis modules."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
KERNELS_DIR = REPO_ROOT / "kernels"

LSK_FILE = KERNELS_DIR / "lsk" / "naif0012.tls"
PCK_TPC_FILE = KERNELS_DIR / "pck" / "pck00011.tpc"
DE432_FILE = KERNELS_DIR / "spk" / "de432s.bsp"
EARTH_BPC_FILE = KERNELS_DIR / "pck" / "earth_000101_241106_240813.bpc"
