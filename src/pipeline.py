"""
FAERS Multi-Quarter Pipeline

Downloads, loads, cleans, and computes signal detection across
multiple quarters of FAERS data. Handles cross-quarter deduplication
and produces quarterly signal snapshots for temporal analysis.
"""

import pandas as pd
from pathlib import Path
from download import download_and_extract
from ingest import load_quarter, load_delete_list
from clean import clean_quarter
from detection import detect_signals


def generate_quarter_list(start_year: int, start_q: int,
                          end_year: int, end_q: int) -> list[tuple[int, int]]:
    """Generate a list of (year, quarter) tuples for a date range."""
    quarters = []
    year, q = start_year, start_q
    while (year, q) <= (end_year, end_q):
        quarters.append((year, q))
        q += 1
        if q > 4:
            q = 1
            year += 1
    return quarters


def download_all_quarters(quarters: list[tuple[int, int]],
                          data_dir: Path = None) -> list[Path]:
    """Download and extract all specified quarters."""
    extract_dirs = []
    for year, q in quarters:
        try:
            extract_dir = download_and_extract(year, q, data_dir)
            extract_dirs.append(extract_dir)
        except Exception as e:
            print(f"  ERROR downloading {year} Q{q}: {e}")
            extract_dirs.append(None)
    return extract_dirs


def process_single_quarter(extract_dir: Path) -> pd.DataFrame | None:
    """Process a single quarter: load, clean, detect signals.

    Returns the signals DataFrame, or None if processing fails.
    """
    try:
        # Load
        tables = load_quarter(extract_dir)
        delete_ids = load_delete_list(extract_dir)

        # Clean
        tables = clean_quarter(tables, delete_ids)

        # Detect signals
        total_cases = tables["demo"]["primaryid"].nunique()
        signals = detect_signals(tables["drug_ps"], tables["reac"], total_cases)

        # Add quarter metadata
        dir_name = extract_dir.name  # e.g., "faers_ascii_2024q1"
        signals["quarter"] = dir_name.replace("faers_ascii_", "").upper()

        return signals

    except Exception as e:
        print(f"  ERROR processing {extract_dir.name}: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_multi_quarter_pipeline(
    start_year: int = 2018,
    start_q: int = 1,
    end_year: int = 2024,
    end_q: int = 1,
    data_dir: Path = None,
    output_dir: Path = None,
) -> pd.DataFrame:
    """Run the complete pipeline across multiple quarters.

    For each quarter:
      1. Download and extract the data
      2. Load, clean, and compute signals
      3. Save quarterly results

    Then combine all quarterly results for temporal analysis.

    Args:
        start_year, start_q: First quarter to process.
        end_year, end_q: Last quarter to process.
        data_dir: Directory for raw data. Defaults to data/raw/.
        output_dir: Directory for results. Defaults to data/processed/.

    Returns:
        Combined DataFrame of all quarterly signals.
    """
    if data_dir is None:
        data_dir = Path("data/raw")
    if output_dir is None:
        output_dir = Path("data/processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    quarters = generate_quarter_list(start_year, start_q, end_year, end_q)
    print(f"\nMulti-Quarter Pipeline")
    print(f"=" * 60)
    print(f"Processing {len(quarters)} quarters: "
          f"{start_year}Q{start_q} through {end_year}Q{end_q}")
    print(f"=" * 60)

    all_signals = []

    for i, (year, q) in enumerate(quarters):
        print(f"\n{'#' * 60}")
        print(f"# QUARTER {i+1}/{len(quarters)}: {year} Q{q}")
        print(f"{'#' * 60}")

        # Download
        try:
            extract_dir = download_and_extract(year, q, data_dir)
        except Exception as e:
            print(f"  ERROR downloading {year} Q{q}: {e}")
            continue

        # Check if we already have results for this quarter
        quarter_label = f"{year}q{q}"
        quarter_file = output_dir / f"signals_{quarter_label}.csv"

        if quarter_file.exists():
            print(f"\n  Already processed - loading cached results from {quarter_file.name}")
            signals = pd.read_csv(quarter_file)
        else:
            # Process
            signals = process_single_quarter(extract_dir)
            if signals is None:
                continue

            # Save quarterly results
            signals.to_csv(quarter_file, index=False)
            print(f"\n  Saved: {quarter_file.name} ({len(signals):,} pairs)")

        all_signals.append(signals)

    # Combine all quarters
    if all_signals:
        combined = pd.concat(all_signals, ignore_index=True)
        combined_path = output_dir / "signals_all_quarters.csv"
        combined.to_csv(combined_path, index=False)
        print(f"\n{'=' * 60}")
        print(f"PIPELINE COMPLETE")
        print(f"  Total quarters processed: {len(all_signals)}")
        print(f"  Total drug-reaction-quarter records: {len(combined):,}")
        print(f"  Saved combined results: {combined_path.name}")
        print(f"{'=' * 60}")
        return combined
    else:
        print("No quarters were successfully processed.")
        return pd.DataFrame()


# -- CLI entry point --
if __name__ == "__main__":
    # Process 2018-2024 (covers Zantac withdrawal 2020,
    # Valsartan recall 2018, Pentosan warning 2020)
    combined = run_multi_quarter_pipeline(
        start_year=2018, start_q=1,
        end_year=2024, end_q=1,
    )
