"""
Builds per-satellite TLE files for a target NORAD list by filtering local historical TLE bundles
(downloaded manually from Space-Track cloud storage, e.g. yearly ZIP files). No Space-Track API calls.
"""

from __future__ import annotations

import re
import sys
import unicodedata
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

INPUT_TXT: str | None = None
DEFAULT_INPUT_FILENAME = "norad_ids_leo_constellations_gnssr_related.txt"

# Path to a downloaded yearly historical TLE ZIP (or a folder containing ZIP/TXT files)
HISTORY_SOURCE_PATH: str | None = None
DEFAULT_HISTORY_SOURCE = "tle_history_2023.zip"

# Optional strict year filter (set to None to keep all epochs found in the source files)
TARGET_YEAR: int | None = 2023

OUTPUT_SUBFOLDER = "LEO"

NORAD_PAD_WIDTH = 0
REPLACE_SPACES_WITH = " "
UPPERCASE_FILENAMES = True
SKIP_IF_FILE_EXISTS = False


@dataclass
class SatRecord:
    norad_id: int
    name: str


def detect_desktop_dir() -> Path:
    home = Path.home()
    candidates = [
        home / "Desktop",
        home / "Pulpit",
        home / "OneDrive" / "Desktop",
        home / "OneDrive" / "Pulpit",
    ]
    for p in candidates:
        if p.exists() and p.is_dir():
            return p
    return home / "Desktop"


def sanitize_name_for_filename(name: str) -> str:
    s = name.strip()
    s = unicodedata.normalize("NFKD", s)
    s = s.encode("ascii", "ignore").decode("ascii")

    if REPLACE_SPACES_WITH:
        s = s.replace(" ", REPLACE_SPACES_WITH)

    sep_repl = REPLACE_SPACES_WITH if REPLACE_SPACES_WITH else "_"
    s = re.sub(r"[^A-Za-z0-9_-]+", sep_repl, s)

    if REPLACE_SPACES_WITH:
        sep = re.escape(REPLACE_SPACES_WITH)
        s = re.sub(rf"{sep}+", REPLACE_SPACES_WITH, s)

    s = s.strip(" _-")
    if not s:
        s = "NONAME"
    if UPPERCASE_FILENAMES:
        s = s.upper()
    return s


def format_norad_for_filename(norad_id: int) -> str:
    s = str(norad_id)
    if NORAD_PAD_WIDTH and NORAD_PAD_WIDTH > 0:
        return s.zfill(NORAD_PAD_WIDTH)
    return s


def resolve_input_file() -> Path:
    if INPUT_TXT:
        p = Path(INPUT_TXT).expanduser()
        if p.exists():
            return p
        raise FileNotFoundError(f"INPUT_TXT file not found: {p}")

    if len(sys.argv) > 1:
        p = Path(sys.argv[1]).expanduser()
        if p.exists():
            return p
        raise FileNotFoundError(f"Input file from CLI argument not found: {p}")

    p = Path(DEFAULT_INPUT_FILENAME).expanduser()
    if p.exists():
        return p

    script_dir = Path(__file__).resolve().parent
    p2 = script_dir / Path(DEFAULT_INPUT_FILENAME).name
    if p2.exists():
        return p2

    p3 = Path.cwd() / Path(DEFAULT_INPUT_FILENAME).name
    if p3.exists():
        return p3

    raise FileNotFoundError(
        "Input file not found. Pass a CLI argument, set INPUT_TXT, or place the file next to the script."
    )


def resolve_history_source() -> Path:
    candidates = []

    if HISTORY_SOURCE_PATH:
        candidates.append(Path(HISTORY_SOURCE_PATH).expanduser())

    script_dir = Path(__file__).resolve().parent
    cwd = Path.cwd()
    default_name = Path(DEFAULT_HISTORY_SOURCE).name
    candidates.extend(
        [
            Path(DEFAULT_HISTORY_SOURCE).expanduser(),
            script_dir / default_name,
            cwd / default_name,
        ]
    )

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "History source not found. Set HISTORY_SOURCE_PATH or place the ZIP/folder next to the script."
    )


def parse_sat_list(txt_path: Path) -> Tuple[List[SatRecord], Dict[int, List[str]]]:
    unique: Dict[int, SatRecord] = {}
    aliases: Dict[int, List[str]] = {}

    with txt_path.open("r", encoding="utf-8", errors="replace") as f:
        for line_no, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if not re.match(r"^\d+\s+", line):
                continue

            m = re.match(r"^(?P<norad>\d+)\s+(?P<name>.+?)\s*$", line)
            if not m:
                print(f"[WARN] Unparsable line {line_no}: {raw.rstrip()}")
                continue

            norad_id = int(m.group("norad"))
            name = m.group("name").strip()

            aliases.setdefault(norad_id, [])
            if name not in aliases[norad_id]:
                aliases[norad_id].append(name)

            if norad_id not in unique:
                unique[norad_id] = SatRecord(norad_id=norad_id, name=name)

    sats = list(unique.values())
    sats.sort(key=lambda s: s.norad_id)
    return sats, aliases


def alpha5_char_value(ch: str) -> int:
    if ch.isdigit():
        return int(ch)
    letters = "ABCDEFGHJKLMNPQRSTUVWXYZ"
    if ch not in letters:
        raise ValueError(f"Invalid Alpha-5 character: {ch}")
    return 10 + letters.index(ch)


def parse_norad_from_tle_line1(line1: str) -> int:
    if len(line1) < 7 or not line1.startswith("1 "):
        raise ValueError("Invalid TLE line 1")
    token = line1[2:7].strip()
    if not token:
        raise ValueError("Empty NORAD token in TLE line 1")
    if token[0].isdigit():
        return int(token)
    if len(token) != 5:
        raise ValueError(f"Unexpected Alpha-5 token length: {token}")
    return alpha5_char_value(token[0]) * 10000 + int(token[1:])


def parse_tle_epoch(line1: str) -> Tuple[int, float]:
    if len(line1) < 32:
        raise ValueError("TLE line 1 too short for epoch parse")
    yy = int(line1[18:20])
    day = float(line1[20:32])
    full_year = 1900 + yy if yy >= 57 else 2000 + yy
    return full_year, day


def iter_tle_pairs_from_text(text: str) -> Iterable[Tuple[str, str]]:
    lines = [ln.rstrip("\r") for ln in text.splitlines() if ln.strip()]
    i = 0
    while i < len(lines):
        if lines[i].startswith("1 ") and i + 1 < len(lines) and lines[i + 1].startswith("2 "):
            yield lines[i], lines[i + 1]
            i += 2
            continue

        if (
            i + 2 < len(lines)
            and not lines[i].startswith("1 ")
            and lines[i + 1].startswith("1 ")
            and lines[i + 2].startswith("2 ")
        ):
            yield lines[i + 1], lines[i + 2]
            i += 3
            continue

        i += 1


def decode_bytes(data: bytes) -> str:
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            return data.decode(enc)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace")


def iter_source_texts(path_value: str) -> Iterable[Tuple[str, str]]:
    p = Path(path_value).expanduser()
    if not p.exists():
        raise FileNotFoundError(f"History source path not found: {p}")

    if p.is_file():
        if p.suffix.lower() == ".zip":
            with zipfile.ZipFile(p, "r") as zf:
                for info in zf.infolist():
                    if info.is_dir():
                        continue
                    with zf.open(info, "r") as fh:
                        yield f"{p.name}:{info.filename}", decode_bytes(fh.read())
        else:
            yield str(p), p.read_text(encoding="utf-8", errors="replace")
        return

    for child in sorted(p.rglob("*")):
        if not child.is_file():
            continue
        suffix = child.suffix.lower()
        if suffix == ".zip":
            with zipfile.ZipFile(child, "r") as zf:
                for info in zf.infolist():
                    if info.is_dir():
                        continue
                    with zf.open(info, "r") as fh:
                        yield f"{child.name}:{info.filename}", decode_bytes(fh.read())
        elif suffix in {".txt", ".tle"}:
            yield str(child), child.read_text(encoding="utf-8", errors="replace")


def write_sat_tle_file(out_dir: Path, sat: SatRecord, tle_text: str) -> Path:
    norad_str = format_norad_for_filename(sat.norad_id)
    name_safe = sanitize_name_for_filename(sat.name)
    out_fp = out_dir / f"{norad_str}_{name_safe}.tle"
    out_fp.write_text(tle_text if tle_text.endswith("\n") else (tle_text + "\n"), encoding="utf-8", newline="\n")
    return out_fp


def main() -> int:
    try:
        txt_path = resolve_input_file()
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        return 1

    try:
        history_source = resolve_history_source()
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        return 1

    sats, aliases = parse_sat_list(txt_path)
    if not sats:
        print("[ERROR] No valid NORAD+NAME entries found in the input file.")
        return 1

    target_norads = {s.norad_id for s in sats}
    desktop_dir = detect_desktop_dir()
    out_dir = desktop_dir / OUTPUT_SUBFOLDER
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Input NORAD list: {txt_path}")
    print(f"[INFO] History source: {history_source}")
    print(f"[INFO] Satellites (unique by NORAD): {len(sats)}")
    print(f"[INFO] Output folder: {out_dir}")
    print(f"[INFO] TARGET_YEAR filter: {TARGET_YEAR if TARGET_YEAR is not None else 'disabled'}")

    dup_alias_count = sum(1 for names in aliases.values() if len(names) > 1)
    if dup_alias_count:
        print(f"[INFO] {dup_alias_count} NORAD IDs have multiple aliases; first name is used for filenames.")

    grouped_pairs: Dict[int, List[Tuple[Tuple[int, float], str, str]]] = {n: [] for n in target_norads}
    seen_pairs: Dict[int, set[Tuple[str, str]]] = {n: set() for n in target_norads}

    scanned_sources = 0
    scanned_pairs = 0
    matched_pairs = 0
    parse_errors = 0

    try:
        for source_name, text in iter_source_texts(str(history_source)):
            scanned_sources += 1
            print(f"[INFO] Scanning source: {source_name}")

            for l1, l2 in iter_tle_pairs_from_text(text):
                scanned_pairs += 1
                try:
                    norad_id = parse_norad_from_tle_line1(l1)
                    epoch_key = parse_tle_epoch(l1)
                except Exception:
                    parse_errors += 1
                    continue

                if norad_id not in target_norads:
                    continue

                if TARGET_YEAR is not None and epoch_key[0] != TARGET_YEAR:
                    continue

                pair_key = (l1, l2)
                if pair_key in seen_pairs[norad_id]:
                    continue

                seen_pairs[norad_id].add(pair_key)
                grouped_pairs[norad_id].append((epoch_key, l1, l2))
                matched_pairs += 1

    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        return 1
    except zipfile.BadZipFile as e:
        print(f"[ERROR] Invalid ZIP file: {e}")
        return 1

    print(f"[INFO] Scanned sources: {scanned_sources}")
    print(f"[INFO] Scanned TLE pairs: {scanned_pairs}")
    print(f"[INFO] Matched TLE pairs: {matched_pairs}")
    if parse_errors:
        print(f"[INFO] Skipped unparsable pairs: {parse_errors}")

    ok = 0
    empty = 0
    skipped = 0
    fail = 0

    for sat in sats:
        out_fp = out_dir / f"{format_norad_for_filename(sat.norad_id)}_{sanitize_name_for_filename(sat.name)}.tle"
        if SKIP_IF_FILE_EXISTS and out_fp.exists():
            skipped += 1
            print(f"SKIP  {sat.norad_id} {sat.name} (already exists)")
            continue

        rows = grouped_pairs.get(sat.norad_id, [])
        if not rows:
            empty += 1
            print(f"EMPTY {sat.norad_id} {sat.name}")
            continue

        rows.sort(key=lambda x: (x[0][0], x[0][1]))
        tle_text = "\n".join(f"{l1}\n{l2}" for _, l1, l2 in rows) + "\n"

        try:
            out_fp = write_sat_tle_file(out_dir, sat, tle_text)
            ok += 1
            print(f"OK    {sat.norad_id} {sat.name} -> {out_fp.name} ({len(rows)} TLE pairs)")
        except Exception as e:
            fail += 1
            print(f"FAIL  {sat.norad_id} {sat.name}: {e}")

    print("\n=== SUMMARY ===")
    print(f"Saved OK:         {ok}")
    print(f"No TLE (EMPTY):   {empty}")
    print(f"Errors (FAIL):    {fail}")
    print(f"Skipped (SKIP):   {skipped}")
    print(f"Output folder:    {out_dir}")

    return 0 if fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
