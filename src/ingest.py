"""
FAERS Data Ingestion Module

Reads raw $-delimited FAERS ASCII files into pandas DataFrames.
Selects only the columns needed for signal detection analysis
and applies basic type conversions.

FAERS relational model:
    DEMO (1 row per case report)
      └── primaryid ──┬── DRUG (many rows per case, one per drug)
                      ├── REAC (many rows per case, one per reaction)
                      └── OUTC (many rows per case, one per outcome)
"""

import pandas as pd
from pathlib import Path


# ── Column selections ────────────────────────────────────────
# We only load the columns relevant to signal detection.
# This saves memory and makes the data easier to work with.

DEMO_COLS = [
    "primaryid",       # Unique report ID (join key)
    "caseid",          # Case ID (shared across duplicate reports)
    "caseversion",     # Version number (higher = more recent update)
    "i_f_code",        # I=Initial, F=Follow-up
    "event_dt",        # Date the adverse event occurred (YYYYMMDD)
    "fda_dt",          # Date FDA received the report (YYYYMMDD)
    "age",             # Patient age
    "age_cod",         # Age unit (YR, MON, DY, etc.)
    "sex",             # Patient sex (M, F, UNK)
    "reporter_country",# Country of the reporter
    "occr_country",    # Country where event occurred
]

DRUG_COLS = [
    "primaryid",       # Join key
    "caseid",
    "drug_seq",        # Sequence number of drug within the case
    "role_cod",        # PS=Primary Suspect, SS=Secondary, C=Concomitant, I=Interacting
    "drugname",        # Drug name as reported (may be brand, generic, or misspelled)
    "prod_ai",         # Active ingredient (more standardized)
]

REAC_COLS = [
    "primaryid",       # Join key
    "caseid",
    "pt",              # MedDRA Preferred Term for the adverse reaction
]

OUTC_COLS = [
    "primaryid",       # Join key
    "caseid",
    "outc_cod",        # Outcome code: DE=Death, HO=Hospitalization, LT=Life-threatening,
                       #   DS=Disability, CA=Congenital anomaly, RI=Required intervention, OT=Other
]


def _find_table_file(extract_dir: Path, table_prefix: str) -> Path:
    """Find the data file for a specific FAERS table within the extraction dir.

    FAERS ZIPs have a nested 'ASCII/' subdirectory containing the .txt files.
    File names follow the pattern: {TABLE_PREFIX}{YY}Q{Q}.txt
    (e.g., DEMO24Q1.txt, DRUG24Q1.txt)

    Args:
        extract_dir: Path to the extracted FAERS quarter directory.
        table_prefix: Table name prefix (e.g., "DEMO", "DRUG").

    Returns:
        Path to the matching .txt file.

    Raises:
        FileNotFoundError: If no matching file is found.
    """
    # Search recursively for the file
    pattern = f"{table_prefix}*.txt"
    matches = list(extract_dir.rglob(pattern)) + list(extract_dir.rglob(pattern.upper()))

    # Filter to only .txt files in ASCII or Deleted subdirectories
    txt_files = [m for m in matches if m.suffix.lower() == ".txt"]

    if not txt_files:
        raise FileNotFoundError(
            f"No {table_prefix}*.txt file found in {extract_dir}"
        )

    # If multiple matches, prefer the one in the ASCII subdirectory
    for f in txt_files:
        if "ascii" in str(f.parent).lower():
            return f

    return txt_files[0]


def load_table(extract_dir: Path, table_prefix: str, columns: list[str]) -> pd.DataFrame:
    """Load a single FAERS table from a quarterly extract.

    Args:
        extract_dir: Path to the extracted FAERS quarter directory.
        table_prefix: Table name prefix (e.g., "DEMO", "DRUG").
        columns: List of column names to keep.

    Returns:
        DataFrame containing only the selected columns.
    """
    filepath = _find_table_file(extract_dir, table_prefix)

    print(f"  Loading {filepath.name} ...", end="", flush=True)

    # Read the $-delimited file
    # - encoding='latin-1' handles international characters
    # - on_bad_lines='warn' skips malformed rows instead of crashing
    # - low_memory=False prevents mixed-type warnings on large files
    df = pd.read_csv(
        filepath,
        sep="$",
        encoding="latin-1",
        low_memory=False,
        on_bad_lines="warn",
    )

    # Standardize column names to lowercase
    df.columns = df.columns.str.strip().str.lower()

    # Select only the columns we need (ignore missing optional columns)
    available_cols = [c for c in columns if c in df.columns]
    missing_cols = [c for c in columns if c not in df.columns]

    if missing_cols:
        print(f" (missing columns: {missing_cols})", end="")

    df = df[available_cols]

    print(f"  {len(df):,} rows")
    return df


def load_quarter(extract_dir: Path) -> dict[str, pd.DataFrame]:
    """Load all relevant FAERS tables for a single quarter.

    Args:
        extract_dir: Path to the extracted FAERS quarter directory
                     (e.g., data/raw/faers_ascii_2024q1/)

    Returns:
        Dictionary with keys 'demo', 'drug', 'reac', 'outc' mapping
        to their respective DataFrames.
    """
    print(f"\nLoading FAERS data from: {extract_dir.name}")
    print("-" * 50)

    tables = {
        "demo": load_table(extract_dir, "DEMO", DEMO_COLS),
        "drug": load_table(extract_dir, "DRUG", DRUG_COLS),
        "reac": load_table(extract_dir, "REAC", REAC_COLS),
        "outc": load_table(extract_dir, "OUTC", OUTC_COLS),
    }

    print("-" * 50)
    print(f"  Total memory usage: "
          f"{sum(df.memory_usage(deep=True).sum() for df in tables.values()) / 1e6:.1f} MB")

    return tables


def load_delete_list(extract_dir: Path) -> set[int]:
    """Load the list of deleted/retracted case IDs for a quarter.

    These are cases that the FDA has removed from the database.
    We should exclude them from our analysis.

    Args:
        extract_dir: Path to the extracted FAERS quarter directory.

    Returns:
        Set of caseid integers that should be excluded.
    """
    try:
        filepath = _find_table_file(extract_dir, "DELETE")
    except FileNotFoundError:
        return set()

    delete_ids = set()
    with open(filepath, "r", encoding="latin-1") as f:
        for line in f:
            line = line.strip()
            if line and line.isdigit():
                delete_ids.add(int(line))

    print(f"  Loaded {len(delete_ids):,} deleted case IDs from {filepath.name}")
    return delete_ids


# ── CLI entry point ──────────────────────────────────────────
if __name__ == "__main__":
    from pathlib import Path

    extract_dir = Path("data/raw/faers_ascii_2024q1")

    # Load all tables
    tables = load_quarter(extract_dir)

    # Load delete list
    deletes = load_delete_list(extract_dir)

    # Print summary statistics for each table
    print("\n" + "=" * 60)
    print("DATA SUMMARY")
    print("=" * 60)

    for name, df in tables.items():
        print(f"\n{name.upper()} table:")
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {list(df.columns)}")
        print(f"  First 3 rows:")
        print(df.head(3).to_string(index=False))

    print(f"\nDELETE list: {len(deletes):,} case IDs to exclude")
    print(f"\nUnique cases in DEMO: {tables['demo']['primaryid'].nunique():,}")
    print(f"Unique drugs in DRUG: {tables['drug']['drugname'].nunique():,}")
    print(f"Unique reactions in REAC: {tables['reac']['pt'].nunique():,}")
