# Drug Safety Signal Detection from FDA Adverse Event Data

> **Research question**: How early and how reliably do statistical disproportionality
> signals (PRR/ROR) emerge in FDA FAERS data for drugs that were later subject to
> regulatory action — and what are the primary sources of noise and false signals?

## Key Findings

1. **Signal detection works.** PRR/ROR consistently identified known drug safety risks across 24 quarters of real FDA data (7.7M+ drug-reaction observations).

2. **Signals preceded regulatory action.** For all 5 case-study drugs, statistical signals were detectable in FAERS data *before* the FDA took action:
   - Pentosan (Elmiron) + maculopathy: **28+ months** before label warning
   - Ranitidine (Zantac) + cancer: **25+ months** before market withdrawal
   - Metformin + NDMA: **27+ months** before FDA investigation
   - Fluoroquinolones + aortic dissection: **10+ months** before safety communication
   - Valsartan + contamination: **5+ months** before recall

3. **Major sources of noise identified:**
   - **Low case counts**: 40% of Evans signals have fewer than 5 cases — extreme PRR but unreliable
   - **Stimulated reporting**: Ranitidine reports surged **21.5x** *after* the FDA withdrawal, inflating post-action signals
   - **Notoriety bias**: Widely-prescribed drugs generate more signals due to higher reporting volume

## Dataset

- **Source**: [FDA AEMS/FAERS Quarterly Data Extracts](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)
- **Coverage**: 2018 Q1 through 2024 Q1 (25 quarters)
- **Scale**: ~400,000 case reports per quarter, 7.7M+ drug-reaction-quarter observations
- **Tables used**: DEMO (demographics), DRUG (drug information), REAC (reactions)

## Methodology

- **Disproportionality analysis**: Proportional Reporting Ratio (PRR) and Reporting Odds Ratio (ROR) with 95% confidence intervals
- **Signal classification**: Evans criteria (PRR ≥ 2, χ² ≥ 4, N ≥ 3)
- **Temporal analysis**: Signal strength tracked across quarterly time windows
- **Validation**: Tested against 5 drugs with documented FDA regulatory actions
- **Noise investigation**: Case count reliability, stimulated reporting, notoriety bias

### Reference

Evans SJW, Waller PC, Davis S. "Use of proportional reporting ratios (PRRs) for signal generation from spontaneous adverse drug reaction reports." *Pharmacoepidemiology and Drug Safety*, 2001.

## Setup

```bash
# Create virtual environment
python -m venv .venv

# Activate (Windows)
.venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Reproduction

```bash
# Step 1: Download FAERS data (25 quarters, ~1.5 GB total)
python src/pipeline.py

# Step 2: Run temporal analysis on case-study drugs
python src/temporal.py

# Step 3: Run noise/false signal analysis
python src/noise_analysis.py
```

Results are cached: once a quarter is processed, re-running skips it.

## Project Structure

```
faers-signal-detection/
├── data/
│   ├── raw/              # Downloaded FAERS ZIP/extracted files (gitignored)
│   └── processed/        # Cleaned signal detection results (gitignored)
├── src/
│   ├── download.py       # FAERS data acquisition from FDA servers
│   ├── ingest.py         # Parsing $-delimited files into DataFrames
│   ├── clean.py          # Deduplication, normalization, PS filtering
│   ├── detection.py      # PRR/ROR computation and Evans classification
│   ├── pipeline.py       # Multi-quarter orchestration with caching
│   ├── temporal.py       # Case study timelines and signal emergence
│   └── noise_analysis.py # False signal investigation
├── figures/              # Generated visualizations
├── requirements.txt
└── README.md
```

## Limitations

- **FAERS is a spontaneous reporting system**; a report does not prove causation.
- **Under-reporting** is inherent; absence of a signal does not mean absence of risk.
- **Cross-quarter deduplication** was applied within quarters only (same caseid, keep latest caseversion). True cross-quarter deduplication would require cumulative processing.
- **Drug name normalization** uses the `prod_ai` (active ingredient) field, which covers 98% of records. The remaining 2% fall back to cleaned `drugname`.
- **Stimulated reporting** means post-action signals are inflated and should not be interpreted as new evidence of risk.

## Technologies

Python · pandas · NumPy · SciPy · matplotlib · seaborn · requests
