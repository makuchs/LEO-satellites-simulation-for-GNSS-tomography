"""misc/sp3_tools.py

Simple, repo-friendly helper script (same style as main.py: one MODE + config block).

This file is optional. It is meant to live in:  misc/sp3_tools.py

It can:
  - download GNSS SP3 products from an IGS MGEX FTP mirror (and unzip)
  - download LEO orbir from UCAR COSMIC data (and extract)
  - combine many SP3 files into one, time-sorted SP3 (optionally keep only selected satellites)
  - sort a dump folder into subfolders based on a filename segment

How to use
----------
1) Set MODE below.
2) Adjust the CONFIG section for that mode.
3) Run:
      python misc/sp3_tools.py

"""

from __future__ import annotations

import gzip
import re
import shutil
import tarfile
from collections import defaultdict
from datetime import date, datetime, timedelta
from ftplib import FTP
from pathlib import Path

try:
    import requests
except ImportError as e:
    raise SystemExit("Missing dependency: requests. Install it with: pip install requests") from e


# =============================================================================
# MODE
# =============================================================================
# Choose one:
#   "download_gnss"   - FTP download GNSS SP3 (MGEX) + unzip
#   "download_leo"    - HTTP download LEO orbit tarballs (UCAR) + extract
#   "combine_sp3"     - combine SP3 files into one, time-sorted file
#   "sort_subfolders" - sort files into subfolders
MODE = "download_gnss"


# =============================================================================
# CONFIG
# =============================================================================

# --- Common paths (repo-relative defaults) ---
REPO_ROOT = Path(__file__).resolve().parents[1]  # .../repo/misc/sp3_tools.py -> repo root
DATA_DIR = REPO_ROOT / "data"
OUTPUT_DIR = REPO_ROOT / "output"

# --- download_gnss ---
GNSS_AC = "COD"  # Analysis Center, e.g. COD, GFZ, GRG...
GNSS_START_DATE = "2021-06-15"  # YYYY-MM-DD
GNSS_END_DATE = "2021-07-06"    # YYYY-MM-DD
GNSS_OUT_DIR = DATA_DIR / "sp3" / "igs"

GNSS_FTP_HOST = "igs.ign.fr"
GNSS_FTP_BASE_DIR = "/pub/igs/products/mgex"
# Optional regex. If None, uses the common MGEX naming:
#   <AC>0MGXFIN_YYYYDDD0000_01D_05M_ORB.SP3.gz
GNSS_FILENAME_REGEX = None
GNSS_KEEP_GZ = False

# --- download_leo ---
# Supported: spire, cosmic2, metopb
LEO_SYSTEM = "spire"
LEO_START_DATE = "2023-09-30"  # YYYY-MM-DD
LEO_END_DATE = "2023-10-03"    # YYYY-MM-DD
LEO_OUT_DIR = DATA_DIR / "sp3" / "leo"
LEO_TIMEOUT_S = 60

# --- combine_sp3 ---
COMBINE_IN_DIR = DATA_DIR / "sp3" / "igs"
COMBINE_OUT_FILE = DATA_DIR / "PG18.sp3"
# If you want only one satellite, set:
#   SAT_IDS = ["PG18"]      (position only)
#   SAT_IDS = ["PG18","VG18"]  (position + velocity)
SAT_IDS = ["PG18"]
SAVE_COMMON_ONLY = False  # keep only sat IDs present in every epoch (useful for messy dumps)

# --- sort_subfolders ---
SORT_IN_DIR = DATA_DIR / "sp3" / "leo"
SEGMENT_INDEX = 2  # filename split by '.' and take this segment (0-based)
DRY_RUN = True     # set False to actually move files


# =============================================================================
# Implementation
# =============================================================================

UCAR_BASE = {
    "spire":  ("https://data.cosmic.ucar.edu/gnss-ro/spire/noaa/nrt/level1b", "leoOrb_nrt"),
    "cosmic2":("https://data.cosmic.ucar.edu/gnss-ro/cosmic2/nrt/level1b", "leoOrb_nrt"),
    "metopb": ("https://data.cosmic.ucar.edu/gnss-ro/metopb/postProc/level1b", "leoOrb_postProc"),
}

SAT_ID_RX = re.compile(r"^([PV][A-Z][0-9]{2})")  # e.g. PG18, VG18, PRN-like records


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def parse_ymd(s: str) -> date:
    return datetime.strptime(s, "%Y-%m-%d").date()


def gpsweek_from_date(d: date) -> int:
    gps_epoch = date(1980, 1, 6)
    return ((d - gps_epoch).days) // 7


def safe_extract_tar(tar: tarfile.TarFile, dst: Path) -> None:
    """Prevent path traversal when extracting tar archives."""
    dst = dst.resolve()
    for member in tar.getmembers():
        member_path = (dst / member.name).resolve()
        if not str(member_path).startswith(str(dst)):
            raise RuntimeError(f"Unsafe path in tar: {member.name}")
    tar.extractall(path=str(dst))


def download_gnss_sp3() -> None:
    ac = GNSS_AC.upper()
    out_dir = ensure_dir(Path(GNSS_OUT_DIR))

    start_date = parse_ymd(GNSS_START_DATE)
    end_date = parse_ymd(GNSS_END_DATE)
    gpsweek_start = gpsweek_from_date(start_date)
    gpsweek_stop = gpsweek_from_date(end_date)

    pattern = GNSS_FILENAME_REGEX or rf"{re.escape(ac)}0MGXFIN.*SP3\\.gz$"
    rx = re.compile(pattern)

    print(f"[download_gnss] AC={ac}, weeks {gpsweek_start}..{gpsweek_stop}")
    print(f"[download_gnss] host={GNSS_FTP_HOST}, base={GNSS_FTP_BASE_DIR}")
    print(f"[download_gnss] regex={pattern}")
    print(f"[download_gnss] out={out_dir}")

    ftp = FTP(GNSS_FTP_HOST, timeout=60)
    ftp.login("", "")  # anonymous

    try:
        for gpsweek in range(gpsweek_start, gpsweek_stop + 1):
            week_dir = f"{GNSS_FTP_BASE_DIR}/{gpsweek}/"
            try:
                ftp.cwd(week_dir)
            except Exception as e:
                print(f"  ! cannot cd to {week_dir}: {e}")
                continue

            try:
                file_list = ftp.nlst()
            except Exception as e:
                print(f"  ! cannot list {week_dir}: {e}")
                continue

            matches = [fn for fn in file_list if rx.search(fn)]
            if not matches:
                print(f"  - week {gpsweek}: no matches")
                continue

            print(f"  - week {gpsweek}: {len(matches)} file(s)")
            for filename in matches:
                gz_path = out_dir / filename
                sp3_path = gz_path.with_suffix("")  # strip .gz

                if sp3_path.exists():
                    print(f"    = exists: {sp3_path.name}")
                    continue

                # download .gz
                with open(gz_path, "wb") as f:
                    print(f"    > downloading: {filename}")
                    ftp.retrbinary(f"RETR {filename}", f.write)

                # unzip to .sp3
                try:
                    with gzip.open(gz_path, "rb") as f_in, open(sp3_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    print(f"    + wrote: {sp3_path.name}")
                finally:
                    if (not GNSS_KEEP_GZ) and gz_path.exists():
                        gz_path.unlink()
    finally:
        try:
            ftp.quit()
        except Exception:
            pass


def download_leo_orbits() -> None:
    system = LEO_SYSTEM.lower()
    if system not in UCAR_BASE:
        raise SystemExit(f"Unknown LEO_SYSTEM={system}. Supported: {', '.join(sorted(UCAR_BASE))}")

    base_url, file_prefix = UCAR_BASE[system]
    out_dir = ensure_dir(Path(LEO_OUT_DIR))

    start_date = parse_ymd(LEO_START_DATE)
    end_date = parse_ymd(LEO_END_DATE)

    print(f"[download_leo] system={system}")
    print(f"[download_leo] base={base_url}")
    print(f"[download_leo] range={start_date} .. {end_date}")
    print(f"[download_leo] out={out_dir}")

    d = start_date
    while d <= end_date:
        year = d.year
        doy = f"{d.timetuple().tm_yday:03d}"
        filename = f"{file_prefix}_{year}_{doy}.tar.gz"
        url = f"{base_url}/{year}/{doy}/{filename}"

        day_dir = ensure_dir(out_dir / f"{year}_{doy}")
        tar_path = day_dir / filename
        
        if any(day_dir.iterdir()) and not tar_path.exists():
            print(f"  = {year}-{doy}: seems extracted already, skipping")
            d += timedelta(days=1)
            continue

        print(f"  > downloading: {url}")
        try:
            with requests.get(url, stream=True, timeout=int(LEO_TIMEOUT_S)) as r:
                if r.status_code == 404:
                    print(f"  - not found (404): {filename}")
                    d += timedelta(days=1)
                    continue
                r.raise_for_status()
                with open(tar_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

            print(f"  > extracting: {tar_path.name}")
            with tarfile.open(tar_path, "r:gz") as tar:
                safe_extract_tar(tar, day_dir)

            tar_path.unlink(missing_ok=True)
            print(f"  + extracted to: {day_dir}")
        except Exception as e:
            print(f"  ! error: {e}")

        d += timedelta(days=1)


def iter_sp3_lines(path: Path):
    if path.suffix.lower() == ".gz" or path.name.lower().endswith(".sp3.gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                yield line
    else:
        with open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                yield line


def combine_sp3_files() -> None:
    in_dir = Path(COMBINE_IN_DIR)
    out_file = Path(COMBINE_OUT_FILE)
    ensure_dir(out_file.parent)

    sat_ids = None if SAT_IDS is None else [s.strip() for s in SAT_IDS if s.strip()]
    files = []
    for p in sorted(in_dir.iterdir()):
        if not p.is_file():
            continue
        name = p.name.lower()
        if name.endswith(".sp3") or name.endswith(".sp3.gz") or name.endswith(".gz"):
            files.append(p)

    if not files:
        raise SystemExit(f"No SP3 files found in: {in_dir}")

    print(f"[combine_sp3] in={in_dir}")
    print(f"[combine_sp3] out={out_file}")
    print(f"[combine_sp3] files={len(files)}")
    print(f"[combine_sp3] sat_ids={sat_ids}")
    print(f"[combine_sp3] save_common_only={bool(SAVE_COMMON_ONLY)}")

    # Header: take all leading header/comment lines from the first file
    header_lines = []
    for line in iter_sp3_lines(files[0]):
        if line.startswith("*"):
            break
        header_lines.append(line)

    # epoch -> {sat_id -> record_line}
    data = defaultdict(dict)

    for p in files:
        epoch = None
        for line in iter_sp3_lines(p):
            if line.startswith("*"):
                epoch = line.rstrip("\n")
                continue
            if epoch is None:
                continue

            m = SAT_ID_RX.match(line)
            if not m:
                continue

            sid = m.group(1)
            if sat_ids and sid not in sat_ids:
                continue

            data[epoch][sid] = line.rstrip("\n")

    if not data:
        raise SystemExit("No epochs/records parsed. Check input files and SAT_IDS.")

    # Optionally keep only sat ids that exist in every epoch
    if SAVE_COMMON_ONLY:
        epochs = list(data.keys())
        common = set(data[epochs[0]].keys())
        for ep in epochs[1:]:
            common &= set(data[ep].keys())
        for ep in list(data.keys()):
            data[ep] = {sid: rec for sid, rec in data[ep].items() if sid in common}
        print(f"[combine_sp3] common sat count: {len(common)}")

    epochs_sorted = sorted(data.keys())

    with open(out_file, "wt", encoding="utf-8", newline="\n") as out:
        out.writelines(header_lines)
        for ep in epochs_sorted:
            sats = data[ep]
            if not sats:
                continue
            out.write(ep + "\n")
            for sid in sorted(sats.keys()):
                out.write(sats[sid] + "\n")
        out.write("EOF\n")

    print(f"[combine_sp3] wrote: {out_file}")


def sort_into_subfolders() -> None:
    in_dir = Path(SORT_IN_DIR)
    ensure_dir(in_dir)

    moved = 0
    for p in sorted(in_dir.iterdir()):
        if not p.is_file():
            continue
        parts = p.name.split(".")
        if len(parts) <= int(SEGMENT_INDEX):
            continue

        frag = parts[int(SEGMENT_INDEX)]
        dest_dir = ensure_dir(in_dir / frag)
        dest = dest_dir / p.name

        if DRY_RUN:
            print(f"[dry_run] {p.name} -> {dest_dir.name}/")
            continue

        shutil.move(str(p), str(dest))
        moved += 1

    print(f"[sort_subfolders] moved {moved} file(s) under: {in_dir}")


def main() -> None:
    if MODE == "download_gnss":
        download_gnss_sp3()
        return

    if MODE == "download_leo":
        download_leo_orbits()
        return

    if MODE == "combine_sp3":
        combine_sp3_files()
        return

    if MODE == "sort_subfolders":
        sort_into_subfolders()
        return

    raise SystemExit(f"Unknown MODE: {MODE}")


if __name__ == "__main__":
    main()