# Drug Safety Signal Detection from FDA Adverse Event Data

> **Research question**: How early and how reliably do statistical disproportionality
> signals (PRR/ROR) emerge in FDA FAERS data for drugs that were later subject to
> regulatory action — and what are the primary sources of noise and false signals?

## Status

🚧 **Under development**

## Dataset

FDA Adverse Event Monitoring System (AEMS/FAERS) quarterly data extracts.

- Source: https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
- Format: `$`-delimited ASCII text files in ZIP archives
- Coverage: 2012–2026 (quarterly releases)

## Setup

```bash
# Create virtual environment
python -m venv .venv

# Activate (Windows)
.venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Project Structure

```
faers-signal-detection/
├── data/
│   ├── raw/          # Downloaded FAERS ZIP/extracted files (gitignored)
│   └── processed/    # Cleaned, deduplicated data (gitignored)
├── src/              # Python source modules
│   ├── download.py   # FAERS data acquisition
│   ├── ingest.py     # Parsing and loading raw files
│   ├── clean.py      # Deduplication and normalization
│   └── signal.py     # PRR/ROR computation
├── notebooks/        # Jupyter analysis notebooks
├── figures/          # Generated visualizations
├── requirements.txt
└── README.md
```

## Methodology

- **Disproportionality analysis**: Proportional Reporting Ratio (PRR) and Reporting
  Odds Ratio (ROR) with 95% confidence intervals
- **Temporal analysis**: Signal emergence tracked across quarterly time windows
- **Validation**: Tested against drugs with known regulatory actions

## Limitations

- FAERS is a spontaneous reporting system; a report does not prove causation.
- Duplicate reports, drug-name inconsistencies, and reporting biases are known
  limitations that are documented and addressed in this analysis.
