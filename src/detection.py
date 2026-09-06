"""
FAERS Signal Detection Module

Computes disproportionality measures for drug-reaction pairs:
  - PRR (Proportional Reporting Ratio)
  - ROR (Reporting Odds Ratio)
  - Chi-squared statistic
  - 95% confidence intervals

These are the standard methods used by pharmacovigilance agencies
(FDA, EMA, WHO-UMC) to detect safety signals in spontaneous
reporting databases.

Reference:
  Evans SJW, Waller PC, Davis S. "Use of proportional reporting ratios
  (PRRs) for signal generation from spontaneous adverse drug reaction
  reports." Pharmacoepidemiology and Drug Safety, 2001.
"""

import numpy as np
import pandas as pd
from scipy import stats


def build_drug_reaction_counts(
    drug_ps: pd.DataFrame,
    reac: pd.DataFrame,
) -> pd.DataFrame:
    """Count the number of cases for each drug-reaction pair.

    Joins the Primary Suspect drug table with the reaction table on primaryid,
    then counts distinct cases for each (drug, reaction) combination.

    Args:
        drug_ps: DRUG table filtered to Primary Suspect only,
                 with 'drug_normalized' column.
        reac: REAC table with 'pt' (Preferred Term) column.

    Returns:
        DataFrame with columns: drug, reaction, pair_count
    """
    print("  Building drug-reaction pair counts...")

    # Join drugs and reactions on primaryid
    # Each case can have multiple drugs and multiple reactions,
    # so this creates all drug-reaction combinations within each case
    merged = drug_ps[["primaryid", "drug_normalized"]].merge(
        reac[["primaryid", "pt"]],
        on="primaryid",
        how="inner",
    )

    # Count unique cases (not rows) for each drug-reaction pair
    pair_counts = (
        merged
        .drop_duplicates(subset=["primaryid", "drug_normalized", "pt"])
        .groupby(["drug_normalized", "pt"])
        .size()
        .reset_index(name="pair_count")
    )

    pair_counts.columns = ["drug", "reaction", "pair_count"]

    print(f"    {len(pair_counts):,} unique drug-reaction pairs found")
    print(f"    Pair count range: {pair_counts['pair_count'].min()} - "
          f"{pair_counts['pair_count'].max()}")

    return pair_counts


def build_marginal_counts(pair_counts: pd.DataFrame, total_cases: int) -> pd.DataFrame:
    """Add marginal counts needed for the 2x2 contingency table.

    For each drug-reaction pair (a), we need:
      a = cases with BOTH this drug AND this reaction
      b = cases with this drug but NOT this reaction
      c = cases with this reaction but NOT this drug
      d = cases with NEITHER this drug NOR this reaction

    These are derived from:
      drug_total = total cases reporting this drug (a + b)
      reaction_total = total cases reporting this reaction (a + c)
      N = total number of cases

    Args:
        pair_counts: DataFrame with drug, reaction, pair_count columns.
        total_cases: Total number of unique cases in the dataset.

    Returns:
        DataFrame with a, b, c, d columns added.
    """
    print("  Computing marginal counts...")

    # Total cases per drug (a + b)
    drug_totals = (
        pair_counts
        .groupby("drug")["pair_count"]
        .sum()
        .reset_index(name="drug_total")
    )

    # Total cases per reaction (a + c)
    reaction_totals = (
        pair_counts
        .groupby("reaction")["pair_count"]
        .sum()
        .reset_index(name="reaction_total")
    )

    # Merge marginals back
    df = pair_counts.merge(drug_totals, on="drug", how="left")
    df = df.merge(reaction_totals, on="reaction", how="left")
    df["N"] = total_cases

    # Build the 2x2 table
    # a = pair_count (cases with both drug and reaction)
    # b = drug_total - a (cases with drug but not this reaction)
    # c = reaction_total - a (cases with reaction but not this drug)
    # d = N - a - b - c (cases with neither)
    df["a"] = df["pair_count"]
    df["b"] = df["drug_total"] - df["a"]
    df["c"] = df["reaction_total"] - df["a"]
    df["d"] = df["N"] - df["a"] - df["b"] - df["c"]

    # Safety check: d should never be negative
    if (df["d"] < 0).any():
        n_negative = (df["d"] < 0).sum()
        print(f"    WARNING: {n_negative} pairs have negative d values (data issue)")
        df.loc[df["d"] < 0, "d"] = 0

    print(f"    Total cases (N): {total_cases:,}")
    print(f"    Pairs with a >= 3: {(df['a'] >= 3).sum():,}")

    return df


def compute_prr(df: pd.DataFrame) -> pd.DataFrame:
    """Compute the Proportional Reporting Ratio (PRR) with 95% CI.

    PRR = (a / (a+b)) / (c / (c+d))

    Interpretation:
      PRR = 1  -> drug-reaction reported at the expected rate
      PRR > 1  -> drug-reaction reported MORE than expected (potential signal)
      PRR < 1  -> drug-reaction reported LESS than expected

    95% CI computed using the log-normal approximation:
      ln(PRR) +/- 1.96 * sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))

    Args:
        df: DataFrame with a, b, c, d columns.

    Returns:
        DataFrame with PRR, PRR_lower, PRR_upper columns added.
    """
    a, b, c, d = df["a"], df["b"], df["c"], df["d"]

    # Avoid division by zero
    # Add 0.5 continuity correction to cells with 0
    a_safe = a.clip(lower=0.5)
    b_safe = b.clip(lower=0.5)
    c_safe = c.clip(lower=0.5)
    d_safe = d.clip(lower=0.5)

    # PRR = (a/(a+b)) / (c/(c+d))
    prop_drug = a_safe / (a_safe + b_safe)
    prop_other = c_safe / (c_safe + d_safe)
    df["PRR"] = prop_drug / prop_other

    # 95% CI using log-normal approximation
    ln_prr = np.log(df["PRR"])
    se_ln_prr = np.sqrt(1/a_safe - 1/(a_safe + b_safe) + 1/c_safe - 1/(c_safe + d_safe))

    df["PRR_lower"] = np.exp(ln_prr - 1.96 * se_ln_prr)
    df["PRR_upper"] = np.exp(ln_prr + 1.96 * se_ln_prr)

    return df


def compute_ror(df: pd.DataFrame) -> pd.DataFrame:
    """Compute the Reporting Odds Ratio (ROR) with 95% CI.

    ROR = (a * d) / (b * c)

    Interpretation is similar to PRR but uses odds instead of proportions.
    ROR is often preferred because it has better statistical properties
    for rare events.

    95% CI using the log-normal approximation:
      ln(ROR) +/- 1.96 * sqrt(1/a + 1/b + 1/c + 1/d)

    Args:
        df: DataFrame with a, b, c, d columns.

    Returns:
        DataFrame with ROR, ROR_lower, ROR_upper columns added.
    """
    a, b, c, d = df["a"], df["b"], df["c"], df["d"]

    # Continuity correction
    a_safe = a.clip(lower=0.5)
    b_safe = b.clip(lower=0.5)
    c_safe = c.clip(lower=0.5)
    d_safe = d.clip(lower=0.5)

    # ROR = (a*d) / (b*c)
    df["ROR"] = (a_safe * d_safe) / (b_safe * c_safe)

    # 95% CI
    ln_ror = np.log(df["ROR"])
    se_ln_ror = np.sqrt(1/a_safe + 1/b_safe + 1/c_safe + 1/d_safe)

    df["ROR_lower"] = np.exp(ln_ror - 1.96 * se_ln_ror)
    df["ROR_upper"] = np.exp(ln_ror + 1.96 * se_ln_ror)

    return df


def compute_chi_squared(df: pd.DataFrame) -> pd.DataFrame:
    """Compute Yates-corrected chi-squared statistic for each pair.

    Chi-squared with Yates correction:
      chi2 = N * (|a*d - b*c| - N/2)^2 / ((a+b)(c+d)(a+c)(b+d))

    A chi-squared >= 4 corresponds roughly to p < 0.05 (1 df).

    Args:
        df: DataFrame with a, b, c, d, N columns.

    Returns:
        DataFrame with chi2 column added.
    """
    a, b, c, d, N = df["a"], df["b"], df["c"], df["d"], df["N"]

    numerator = N * (np.abs(a * d - b * c) - N / 2) ** 2
    denominator = (a + b) * (c + d) * (a + c) * (b + d)

    # Avoid division by zero
    denominator = denominator.replace(0, np.nan)

    df["chi2"] = numerator / denominator

    return df


def classify_signals(df: pd.DataFrame) -> pd.DataFrame:
    """Classify drug-reaction pairs as signals using Evans criteria.

    A pair is flagged as a statistical signal if ALL of:
      1. PRR >= 2  (at least twice the expected rate)
      2. chi2 >= 4  (statistically significant at ~p < 0.05)
      3. a >= 3     (at least 3 case reports of the pair)

    This is intentionally conservative to reduce false positives.

    Also flags based on the lower 95% CI of PRR being > 1
    (a stricter criterion used by some agencies).

    Args:
        df: DataFrame with PRR, chi2, a columns.

    Returns:
        DataFrame with 'signal_evans' and 'signal_ci' boolean columns.
    """
    # Evans criteria (standard)
    df["signal_evans"] = (
        (df["PRR"] >= 2) &
        (df["chi2"] >= 4) &
        (df["a"] >= 3)
    )

    # CI-based criterion (stricter): lower bound of 95% CI > 1
    df["signal_ci"] = (
        (df["PRR_lower"] > 1) &
        (df["a"] >= 3)
    )

    evans_count = df["signal_evans"].sum()
    ci_count = df["signal_ci"].sum()
    total = len(df)

    print(f"  Signal classification:")
    print(f"    Evans criteria (PRR>=2, chi2>=4, a>=3): "
          f"{evans_count:,} signals ({evans_count/total*100:.1f}%)")
    print(f"    CI-based (PRR_lower>1, a>=3): "
          f"{ci_count:,} signals ({ci_count/total*100:.1f}%)")

    return df


def detect_signals(
    drug_ps: pd.DataFrame,
    reac: pd.DataFrame,
    total_cases: int,
) -> pd.DataFrame:
    """Run the complete signal detection pipeline.

    This is the main entry point that chains together all steps:
    pair counting -> marginals -> PRR -> ROR -> chi2 -> classification.

    Args:
        drug_ps: DRUG table (Primary Suspect only, with drug_normalized).
        reac: REAC table.
        total_cases: Total number of unique cases.

    Returns:
        DataFrame of all drug-reaction pairs with signal scores.
    """
    print("\n" + "=" * 60)
    print("SIGNAL DETECTION")
    print("=" * 60)

    # Step 1: Count drug-reaction pairs
    pair_counts = build_drug_reaction_counts(drug_ps, reac)

    # Step 2: Compute marginal counts for 2x2 table
    signals = build_marginal_counts(pair_counts, total_cases)

    # Step 3: Compute PRR
    print("\n  Computing PRR...")
    signals = compute_prr(signals)

    # Step 4: Compute ROR
    print("  Computing ROR...")
    signals = compute_ror(signals)

    # Step 5: Compute chi-squared
    print("  Computing chi-squared...")
    signals = compute_chi_squared(signals)

    # Step 6: Classify signals
    print()
    signals = classify_signals(signals)

    # Sort by PRR descending for readability
    signals = signals.sort_values("PRR", ascending=False).reset_index(drop=True)

    print(f"\n  Top 10 signals by PRR (a >= 3):")
    top = signals[signals["a"] >= 3].head(10)
    for _, row in top.iterrows():
        print(f"    {row['drug']:30s} + {row['reaction']:30s}  "
              f"PRR={row['PRR']:8.1f}  ROR={row['ROR']:8.1f}  "
              f"a={int(row['a']):4d}  chi2={row['chi2']:8.1f}")

    print("=" * 60)
    return signals


# -- CLI entry point --
if __name__ == "__main__":
    from pathlib import Path
    from ingest import load_quarter, load_delete_list
    from clean import clean_quarter

    extract_dir = Path("data/raw/faers_ascii_2024q1")

    # Load
    tables = load_quarter(extract_dir)
    delete_ids = load_delete_list(extract_dir)

    # Clean
    tables = clean_quarter(tables, delete_ids)

    # Detect signals
    total_cases = tables["demo"]["primaryid"].nunique()
    signals = detect_signals(tables["drug_ps"], tables["reac"], total_cases)

    # Save results
    output_path = Path("data/processed/signals_2024q1.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    signals.to_csv(output_path, index=False)
    print(f"\nResults saved to: {output_path}")
    print(f"Total pairs: {len(signals):,}")
    print(f"Evans signals: {signals['signal_evans'].sum():,}")
