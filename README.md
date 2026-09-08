# Drug Safety Signal Detection from FDA Adverse Event Data

> **Research question**: How early and how reliably do statistical disproportionality
> signals (PRR/ROR) appear in FDA FAERS data for drugs that were later subject to
> regulatory action — and what are the primary sources of estimate instability
> and potential confounding?

## Key Findings

1. **Disproportionality signals are retrospectively consistent with known safety events.** PRR/ROR signals for clinically relevant reactions were present in FAERS data for all 5 case-study drugs. Evans signals averaged ~6.4% of drug-reaction pairs per quarter.

2. **Signals were retrospectively observable before regulatory action dates.** For case-study drugs with well-defined relevant reactions, Evans signals for those reactions were present in FAERS quarterly data in quarters preceding the regulatory action:
   - Pentosan (Elmiron) + retinal reactions: present **19 months** before label warning
   - Fluoroquinolones + aortic reactions: present **10 months** before safety communication
   - Ranitidine (Zantac) + cancer reactions: present **8 months** before market withdrawal
   - Valsartan + contamination reactions: present **2 months** before recall
   - Metformin: peak signal was for lactic acidosis (a known pre-existing risk, not NDMA-related) — **weakest case study**

   > These are retrospective observations. Case-study drugs were selected because
   > their regulatory actions are already known. This does not demonstrate
   > prospective predictive capability.

3. **Sources of estimate instability identified:**
   - **Low case counts**: 40% of Evans signals are based on 3–4 case reports. These may reflect genuine safety concerns, but their PRR estimates have wide confidence intervals and should be interpreted cautiously.
   - **Reporting-volume surges**: Ranitidine-related Evans signals increased ~3.3x in quarters following the FDA withdrawal request, consistent with stimulated reporting.
   - **Reporting-volume effects**: Drugs with broad adverse effect profiles (immunosuppressants, antipsychotics) generate more signals. Without prescription-volume data, we cannot separate genuine multi-reaction risk from reporting-volume artifacts.

4. **Negative-control comparison**: Widely-prescribed drugs without major regulatory actions (levothyroxine, omeprazole, amlodipine) also produce persistent Evans signals in every quarter. This illustrates that Evans signals alone do not indicate a drug safety problem requiring regulatory action.

## Dataset

- **Source**: [FDA FAERS Quarterly Data Extracts](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html)
- **Coverage**: 2018 Q1 through 2024 Q1 (25 quarters downloaded; 24 quarters analyzed)
- **Scale**: ~400,000 case reports per quarter, 7.7M+ drug-reaction-quarter observations
- **Tables used**: DEMO (demographics), DRUG (drug information), REAC (reactions)
- **Temporal unit**: FAERS quarterly data file (not individual event date)

## Methodology

- **Disproportionality analysis**: Proportional Reporting Ratio (PRR) and Reporting Odds Ratio (ROR) with 95% confidence intervals via log-normal approximation
- **Signal classification**: Evans criteria (PRR ≥ 2, χ² ≥ 4, a ≥ 3), where a = number of case reports for the drug-reaction pair
- **Temporal analysis**: Signal presence tracked across quarterly time windows for clinically relevant reactions
- **Retrospective comparison**: Compared against 5 drugs with documented FDA regulatory actions
- **Negative controls**: 3 widely-prescribed drugs without major regulatory actions (2018–2024)
- **Limitations analysis**: Case count stability, reporting-volume effects

### Reference

Evans SJW, Waller PC, Davis S. "Use of proportional reporting ratios (PRRs) for signal generation from spontaneous adverse drug reaction reports." *Pharmacoepidemiology and Drug Safety*, 2001; 10:483–486.

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
# Step 1: Download and process FAERS data (25 quarters, ~1.5 GB total)
python src/pipeline.py

# Step 2: Run temporal analysis on case-study drugs
python src/temporal.py

# Step 3: Run signal reliability and limitations analysis
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
│   ├── temporal.py       # Retrospective case study timelines
│   └── noise_analysis.py # Signal reliability and limitations analysis
├── figures/              # Generated visualizations
├── requirements.txt
└── README.md
```

## Limitations

- **FAERS is a spontaneous reporting system**; a report does not prove causation and the database lacks a true exposure denominator.
- **Under-reporting** is inherent; absence of a signal does not mean absence of risk.
- **Cross-quarter deduplication**: Each quarter's data was processed independently. Cases updated across quarters (same `caseid`, higher `caseversion`) may appear in multiple quarters, which could moderately affect temporal signal persistence estimates.
- **Drug name normalization** uses the `prod_ai` (active ingredient) field, which covers 98% of records. The remaining 2% fall back to cleaned `drugname`.
- **Retrospective case-study selection**: Case-study drugs were selected because their regulatory actions are already known. The analysis does not demonstrate prospective prediction.
- **No exposure denominator**: Without prescription-volume data, high signal counts for widely-used drugs cannot be distinguished from genuine multi-reaction risk profiles.

## Future Work

- Cross-quarter deduplication using a cumulative caseid registry
- MedDRA hierarchy grouping for related reactions
- Prescription-volume normalization using external exposure data
- Formal negative-control analysis with more drugs and false-positive-rate estimation

## Technologies

Python · pandas · NumPy · SciPy · matplotlib · seaborn · requests
