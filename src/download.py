"""
FAERS Data Download Module

Downloads quarterly ASCII data extracts from the FDA AEMS/FAERS portal.
Each quarterly ZIP contains $-delimited text files for:
  DEMO (demographics), DRUG (drugs), REAC (reactions),
  OUTC (outcomes), INDI (indications), THER (therapy), RPSR (report sources)
"""

import os
import zipfile
import requests
from pathlib import Path
from tqdm import tqdm


# Base URL for FAERS ASCII quarterly data files
FAERS_BASE_URL = "https://fis.fda.gov/content/Exports"

# Default directory for downloaded data
DEFAULT_DATA_DIR = Path(__file__).parent.parent / "data" / "raw"


def build_download_url(year: int, quarter: int) -> str:
    """Build the download URL for a specific FAERS quarterly file.

    Args:
        year: Four-digit year (e.g. 2024)
        quarter: Quarter number (1, 2, 3, or 4)

    Returns:
        Full URL to the ZIP file on the FDA server.
    """
    # FDA uses inconsistent casing in filenames (Q4 vs q1),
    # but lowercase works for all recent files
    filename = f"faers_ascii_{year}q{quarter}.zip"
    return f"{FAERS_BASE_URL}/{filename}"


def download_quarter(year: int, quarter: int, data_dir: Path = None) -> Path:
    """Download a single FAERS quarterly ZIP file.

    Args:
        year: Four-digit year
        quarter: Quarter number (1-4)
        data_dir: Directory to save the file. Defaults to data/raw/.

    Returns:
        Path to the downloaded ZIP file.
    """
    if data_dir is None:
        data_dir = DEFAULT_DATA_DIR

    data_dir.mkdir(parents=True, exist_ok=True)

    url = build_download_url(year, quarter)
    zip_filename = f"faers_ascii_{year}q{quarter}.zip"
    zip_path = data_dir / zip_filename

    # Skip if already downloaded
    if zip_path.exists():
        print(f"  Already downloaded: {zip_path.name}")
        return zip_path

    print(f"  Downloading {url} ...")
    response = requests.get(url, stream=True)
    response.raise_for_status()

    # Get total file size for progress bar
    total_size = int(response.headers.get("content-length", 0))

    with open(zip_path, "wb") as f:
        with tqdm(total=total_size, unit="B", unit_scale=True, desc=zip_filename) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))

    print(f"  Saved: {zip_path}")
    return zip_path


def extract_quarter(zip_path: Path) -> Path:
    """Extract a FAERS quarterly ZIP file.

    Args:
        zip_path: Path to the ZIP file.

    Returns:
        Path to the extraction directory.
    """
    extract_dir = zip_path.parent / zip_path.stem
    if extract_dir.exists():
        print(f"  Already extracted: {extract_dir.name}/")
        return extract_dir

    print(f"  Extracting {zip_path.name} ...")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)

    print(f"  Extracted to: {extract_dir}")
    return extract_dir


def list_extracted_files(extract_dir: Path) -> list[Path]:
    """List all data files in an extracted FAERS directory.

    Returns:
        Sorted list of file paths.
    """
    # FAERS extracts often have a nested directory inside the ZIP
    all_files = sorted(extract_dir.rglob("*.txt")) + sorted(extract_dir.rglob("*.TXT"))
    return all_files


def download_and_extract(year: int, quarter: int, data_dir: Path = None) -> Path:
    """Download and extract a FAERS quarterly file.

    This is the main convenience function.

    Args:
        year: Four-digit year
        quarter: Quarter number (1-4)
        data_dir: Directory for downloads. Defaults to data/raw/.

    Returns:
        Path to the extraction directory.
    """
    print(f"\n{'='*60}")
    print(f"FAERS {year} Q{quarter}")
    print(f"{'='*60}")

    zip_path = download_quarter(year, quarter, data_dir)
    extract_dir = extract_quarter(zip_path)

    files = list_extracted_files(extract_dir)
    print(f"\n  Files found ({len(files)}):")
    for f in files:
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"    {f.name:30s}  {size_mb:7.1f} MB")

    return extract_dir


# ── CLI entry point ──────────────────────────────────────────
if __name__ == "__main__":
    # Download ONE quarter to start — we'll inspect it before scaling up
    print("FAERS Data Acquisition")
    print("=" * 60)
    print("Downloading a single quarter to inspect data structure.")
    print()

    extract_dir = download_and_extract(year=2024, quarter=1)

    print()
    print("Done! Next step: inspect the raw data files.")
