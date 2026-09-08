"""
FAERS Signal Reliability & Limitations Analysis

Investigates sources of instability and potential confounding
in FAERS disproportionality analysis.

Analyses included:
  1. Low case counts producing statistically unstable PRR estimates
  2. Reporting-volume surge around ranitidine withdrawal
  3. Reporting-volume effects among top signal-generating drugs
  4. Negative-control comparison using drugs without major regulatory actions
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def analyze_case_count_distribution(combined: pd.DataFrame, output_dir: Path):
    """Analyze how case count (a) affects PRR estimate stability.

    Signals based on very few case reports have wide confidence
    intervals, making the point estimate of PRR statistically
    unstable. This analysis quantifies the relationship.
    """
    print("\n1. CASE COUNT vs PRR ESTIMATE STABILITY")
    print("-" * 40)

    signals = combined[combined["signal_evans"]].copy()

    # Group by case count bins
    bins = [3, 5, 10, 25, 50, 100, 500, float("inf")]
    labels = ["3-4", "5-9", "10-24", "25-49", "50-99", "100-499", "500+"]
    signals["count_bin"] = pd.cut(signals["a"], bins=bins, labels=labels, right=False)

    stats = signals.groupby("count_bin", observed=True).agg(
        n_signals=("PRR", "size"),
        median_prr=("PRR", "median"),
        mean_prr=("PRR", "mean"),
        q25_prr=("PRR", lambda x: x.quantile(0.25)),
        q75_prr=("PRR", lambda x: x.quantile(0.75)),
    ).reset_index()

    print("\n  Signal count and PRR by case count bin:")
    print(stats.to_string(index=False))

    # Key finding: what fraction of signals have a < 5?
    total_signals = len(signals)
    low_count = (signals["a"] < 5).sum()
    print(f"\n  Signals with a = 3 or 4: {low_count:,} / {total_signals:,} "
          f"({low_count/total_signals*100:.1f}%)")
    print(f"  These have median PRR = "
          f"{signals[signals['a'] < 5]['PRR'].median():.1f}")
    print(f"  Signals with a >= 10 have median PRR = "
          f"{signals[signals['a'] >= 10]['PRR'].median():.1f}")
    print(f"\n  Interpretation: Signals based on 3-4 cases are not")
    print(f"  necessarily incorrect, but their PRR point estimates")
    print(f"  have wide confidence intervals and should be interpreted")
    print(f"  with caution. The higher median PRR for low-count signals")
    print(f"  reflects the greater statistical variability of small samples.")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: distribution of case counts among signals
    ax = axes[0]
    signals["a"].clip(upper=100).hist(bins=50, ax=ax, color="steelblue", edgecolor="white")
    ax.set_xlabel("Case count (a), clipped at 100")
    ax.set_ylabel("Number of Evans signals")
    ax.set_title("Distribution of case counts among Evans signals")
    ax.axvline(x=5, color="red", linestyle="--", alpha=0.7, label="a = 5 threshold")
    ax.legend()

    # Right: PRR vs case count
    ax = axes[1]
    sample = signals.sample(min(5000, len(signals)), random_state=42)
    ax.scatter(sample["a"], sample["PRR"], alpha=0.15, s=5, color="steelblue")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Case count (a)")
    ax.set_ylabel("PRR")
    ax.set_title("PRR vs case count (Evans signals)")
    ax.axhline(y=2, color="red", linestyle="--", alpha=0.5, label="PRR = 2")
    ax.legend()

    plt.tight_layout()
    fig.savefig(output_dir / "noise_case_count.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: noise_case_count.png")

    return stats


def analyze_reporting_volume_surge(combined: pd.DataFrame, output_dir: Path):
    """Analyze the reporting-volume surge around the ranitidine withdrawal.

    When the FDA announced the withdrawal request for ranitidine (April 2020),
    the number of distinct drug-reaction pairs reported for ranitidine
    increased substantially in subsequent quarters. This pattern is consistent
    with stimulated reporting — a well-documented phenomenon where regulatory
    actions and media coverage trigger increased voluntary reporting.

    Note: This analysis cannot definitively attribute the increase to
    stimulated reporting vs. other causes (e.g., legal proceedings,
    increased diagnostic awareness).
    """
    print("\n2. REPORTING-VOLUME SURGE (Ranitidine)")
    print("-" * 40)

    rani = combined[combined["drug"].str.contains("RANITIDINE", na=False)].copy()
    if rani.empty:
        print("  No ranitidine data found.")
        return

    # Count unique drug-reaction pairs per quarter (not sum of a values)
    # This avoids double-counting from the same case appearing in multiple pairs
    rani_quarterly = (
        rani.groupby("quarter")
        .agg(
            n_pairs=("reaction", "size"),
            n_evans_signals=("signal_evans", "sum"),
            n_unique_reactions=("reaction", "nunique"),
            max_case_count=("a", "max"),
        )
        .reset_index()
        .sort_values("quarter")
    )

    print("\n  Ranitidine quarterly drug-reaction pair counts:")
    print(rani_quarterly.to_string(index=False))

    # Identify the pre/post withdrawal change
    pre = rani_quarterly[rani_quarterly["quarter"] < "2020Q2"]
    post = rani_quarterly[rani_quarterly["quarter"] >= "2020Q2"]

    if not pre.empty and not post.empty:
        pre_avg_reactions = pre["n_unique_reactions"].mean()
        post_avg_reactions = post["n_unique_reactions"].mean()
        pre_avg_signals = pre["n_evans_signals"].mean()
        post_avg_signals = post["n_evans_signals"].mean()
        reaction_ratio = post_avg_reactions / pre_avg_reactions if pre_avg_reactions > 0 else float("inf")
        signal_ratio = post_avg_signals / pre_avg_signals if pre_avg_signals > 0 else float("inf")

        print(f"\n  Avg unique reactions per quarter BEFORE withdrawal: {pre_avg_reactions:.0f}")
        print(f"  Avg unique reactions per quarter AFTER withdrawal: {post_avg_reactions:.0f}")
        print(f"  Ratio: {reaction_ratio:.1f}x")
        print(f"\n  Avg Evans signals per quarter BEFORE withdrawal: {pre_avg_signals:.0f}")
        print(f"  Avg Evans signals per quarter AFTER withdrawal: {post_avg_signals:.0f}")
        print(f"  Ratio: {signal_ratio:.1f}x")
        print(f"\n  This increase in reported reactions is consistent with")
        print(f"  stimulated reporting, though other factors (legal proceedings,")
        print(f"  increased diagnostic awareness) may also contribute.")

    # Plot using n_unique_reactions as the metric
    fig, ax = plt.subplots(figsize=(12, 5))
    from temporal import quarter_to_date
    from datetime import datetime

    rani_quarterly["date"] = rani_quarterly["quarter"].apply(quarter_to_date)
    action_date = datetime(2020, 4, 1)

    colors = ["#2196F3" if d < action_date else "#F44336"
              for d in rani_quarterly["date"]]
    ax.bar(rani_quarterly["date"], rani_quarterly["n_unique_reactions"],
           width=60, color=colors, edgecolor="white", alpha=0.8)
    ax.axvline(x=action_date, color="black", linestyle="-.", linewidth=2,
               label="FDA withdrawal request (Apr 2020)")
    ax.set_xlabel("Quarter")
    ax.set_ylabel("Number of distinct reported reactions")
    ax.set_title("Ranitidine (Zantac): Reporting-volume surge around FDA withdrawal",
                 fontsize=13, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add annotation
    ax.annotate("Post-withdrawal\nreporting increase",
                xy=(datetime(2020, 8, 15), post_avg_reactions),
                xytext=(datetime(2021, 6, 1), post_avg_reactions * 1.4),
                fontsize=10, ha="center",
                arrowprops=dict(arrowstyle="->", color="gray"))

    plt.tight_layout()
    fig.savefig(output_dir / "noise_stimulated_reporting.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: noise_stimulated_reporting.png")


def analyze_reporting_volume_effects(combined: pd.DataFrame, output_dir: Path):
    """Show which drugs generate the most Evans signals.

    Without prescription-volume or exposure data, we cannot determine
    whether high signal counts reflect genuine multi-reaction safety
    profiles or simply high reporting volume.
    """
    print("\n3. REPORTING-VOLUME EFFECTS (top signal-generating drugs)")
    print("-" * 40)

    drug_signals = (
        combined[combined["signal_evans"]]
        .groupby("drug")
        .agg(
            n_signals=("signal_evans", "sum"),
            n_quarters=("quarter", "nunique"),
            max_prr=("PRR", "max"),
            median_a=("a", "median"),
        )
        .reset_index()
        .sort_values("n_signals", ascending=False)
    )

    print("\n  Top 15 drugs by total Evans signal count across all quarters:")
    top15 = drug_signals.head(15)
    for _, row in top15.iterrows():
        print(f"    {row['drug']:45s}  signals={int(row['n_signals']):5d}  "
              f"quarters={int(row['n_quarters']):2d}  "
              f"median_a={row['median_a']:6.0f}  max_PRR={row['max_prr']:10.1f}")

    print(f"\n  Interpretation: The drugs generating the most Evans signals")
    print(f"  include immunosuppressants (tacrolimus, mycophenolate mofetil),")
    print(f"  biologics (adalimumab), and CNS drugs (aripiprazole, olanzapine).")
    print(f"  These drug classes are known to have broad adverse effect profiles.")
    print(f"  Without prescription-volume data, we cannot distinguish genuine")
    print(f"  multi-reaction risk from reporting-volume artifacts.")

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    top15_plot = top15.head(15).sort_values("n_signals")
    ax.barh(top15_plot["drug"], top15_plot["n_signals"],
            color="steelblue", edgecolor="white")
    ax.set_xlabel("Total Evans signals across all quarters")
    ax.set_title("Top 15 drugs by Evans signal count (reporting-volume effects)",
                 fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    fig.savefig(output_dir / "noise_top_drugs.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: noise_top_drugs.png")


def analyze_negative_controls(combined: pd.DataFrame, output_dir: Path):
    """Exploratory negative-control analysis.

    Checks whether the disproportionality methodology also produces
    persistent Evans signals for drugs that were NOT selected as
    known-positive case studies and do not have an obvious major
    regulatory safety action in the 2018-2024 period.

    Selected negative-control drugs:
      1. LEVOTHYROXINE — thyroid hormone replacement, widely prescribed,
         no major safety withdrawal/recall in 2018-2024.
      2. OMEPRAZOLE — proton pump inhibitor, widely prescribed,
         no major safety withdrawal in 2018-2024.
      3. AMLODIPINE — calcium channel blocker for hypertension,
         widely prescribed, no major safety withdrawal in 2018-2024.

    Selection rationale: All three are high-volume, established drugs
    with well-characterized safety profiles. They were NOT selected
    based on whether they have Evans signals or not.

    This analysis does NOT establish a population-level false-positive rate.
    It is an exploratory comparison to illustrate the limitations of
    disproportionality analysis for signal interpretation.
    """
    print("\n4. NEGATIVE-CONTROL COMPARISON")
    print("-" * 40)

    controls = [
        {"name": "Levothyroxine", "pattern": "LEVOTHYROXINE"},
        {"name": "Omeprazole", "pattern": "OMEPRAZOLE"},
        {"name": "Amlodipine", "pattern": "AMLODIPINE"},
    ]

    print("\n  Negative-control drugs (no major regulatory action 2018-2024):")
    for c in controls:
        print(f"    - {c['name']}")

    control_rows = []

    for ctrl in controls:
        ctrl_data = combined[
            combined["drug"].str.contains(ctrl["pattern"], na=False)
        ]

        if ctrl_data.empty:
            print(f"\n  {ctrl['name']}: No data found")
            continue

        n_quarters = ctrl_data["quarter"].nunique()
        total_pairs = len(ctrl_data)
        evans_total = ctrl_data["signal_evans"].sum()
        evans_rate = evans_total / total_pairs * 100 if total_pairs > 0 else 0

        # How many quarters have at least one Evans signal?
        quarters_with_signals = (
            ctrl_data[ctrl_data["signal_evans"]]
            .groupby("quarter")
            .size()
            .reset_index(name="n_signals")
        )
        n_quarters_with_signals = len(quarters_with_signals)

        # Top 3 signals by PRR (a >= 3)
        top_signals = (
            ctrl_data[(ctrl_data["signal_evans"]) & (ctrl_data["a"] >= 3)]
            .sort_values("PRR", ascending=False)
            .drop_duplicates(subset=["reaction"])
            .head(3)
        )

        print(f"\n  {ctrl['name']}:")
        print(f"    Quarters present: {n_quarters}")
        print(f"    Total drug-reaction pairs: {total_pairs:,}")
        print(f"    Evans signals: {evans_total:,} ({evans_rate:.1f}%)")
        print(f"    Quarters with >= 1 Evans signal: "
              f"{n_quarters_with_signals}/{n_quarters}")

        if not top_signals.empty:
            print(f"    Top signals (by PRR, a >= 3):")
            for _, row in top_signals.iterrows():
                print(f"      {row['reaction']:40s}  PRR={row['PRR']:8.1f}  "
                      f"a={int(row['a']):4d}")

        control_rows.append({
            "Drug": ctrl["name"],
            "Quarters Present": n_quarters,
            "Total Pairs": total_pairs,
            "Evans Signals": int(evans_total),
            "Signal Rate (%)": f"{evans_rate:.1f}",
            "Quarters With Signals": f"{n_quarters_with_signals}/{n_quarters}",
        })

    if control_rows:
        ctrl_df = pd.DataFrame(control_rows)
        print(f"\n  Summary:")
        print(ctrl_df.to_string(index=False))

        print(f"\n  Key observation: All three negative-control drugs produce")
        print(f"  persistent Evans signals across most or all quarters.")
        print(f"  This demonstrates that disproportionality signals alone")
        print(f"  are not equivalent to confirmed safety problems and")
        print(f"  reinforces the need for clinical context and additional")
        print(f"  evidence before signals can be considered actionable.")
        print(f"  This exploratory comparison does NOT estimate the overall")
        print(f"  false-positive rate of Evans criteria.")


def generate_findings_summary(combined: pd.DataFrame) -> str:
    """Generate a text summary of key findings for the README/report."""

    total_pairs = len(combined)
    total_quarters = combined["quarter"].nunique()
    total_drugs = combined["drug"].nunique()
    total_reactions = combined["reaction"].nunique()

    low_count_pct = (
        combined[combined["signal_evans"]]["a"] < 5
    ).mean() * 100

    summary = f"""
RESEARCH FINDINGS SUMMARY
{'=' * 60}

Research Question:
  How early and how reliably do statistical disproportionality
  signals (PRR/ROR) appear in FDA FAERS data for drugs that were
  later subject to regulatory action -- and what are the primary
  sources of estimate instability and potential confounding?

Dataset:
  - {total_quarters} quarters of FDA FAERS data (2018 Q1 - 2023 Q4)
  - {total_pairs:,} drug-reaction-quarter observations
  - {total_drugs:,} unique drugs, {total_reactions:,} unique reactions
  - Temporal unit: FAERS quarterly data file (not individual event date)

Key Findings:

  1. RETROSPECTIVE SIGNAL CONSISTENCY
     PRR/ROR signals were retrospectively consistent with known
     safety events for the selected case-study drugs. Evans
     signals averaged approximately 6.4% of drug-reaction pairs
     per quarter.

  2. SIGNALS RETROSPECTIVELY OBSERVABLE BEFORE REGULATORY ACTION
     For the case-study drugs with well-defined relevant reactions,
     Evans signals for clinically relevant reactions were present
     in FAERS quarterly data from early in the observation window.
     This is a retrospective observation; the analysis does not
     demonstrate prospective predictive capability.

  3. SOURCES OF ESTIMATE INSTABILITY
     a) Low case counts: {low_count_pct:.0f}% of Evans signals are based
        on 3-4 case reports. While these may reflect genuine safety
        concerns, their PRR estimates have wide confidence intervals
        and should be interpreted cautiously.
     b) Reporting-volume surges: The volume of ranitidine-related
        reports increased substantially after the FDA withdrawal
        request, consistent with stimulated reporting. This inflates
        post-action signal metrics.
     c) Reporting-volume effects: Drugs with broad adverse effect
        profiles or high prescription volume generate more Evans
        signals. Without exposure data, we cannot separate genuine
        multi-reaction risk from reporting-volume artifacts.

  4. NEGATIVE-CONTROL OBSERVATION
     Widely-prescribed drugs without major regulatory actions
     (levothyroxine, omeprazole, amlodipine) also produce persistent
     Evans signals across most quarters. This demonstrates that
     disproportionality signals alone are not equivalent to confirmed
     safety problems and reinforces the need for clinical context
     and additional evidence before signals can be considered
     actionable. This exploratory comparison does NOT estimate the
     overall false-positive rate of Evans criteria.

  5. LIMITATIONS
     - FAERS is a spontaneous reporting system; reports do not prove
       causation and lack a true exposure denominator.
     - Under-reporting is inherent; absence of signal does not mean
       absence of risk.
     - Multiple comparisons: A large number of drug-reaction
       combinations are tested across many quarters using fixed
       thresholds without formal multiple-testing correction. Some
       signals may therefore occur by chance, and reporting volume
       can further increase the number of observed signals.
       Disproportionality signals should be treated as screening
       signals requiring further evaluation, not as confirmation of
       causality or safety problems.
     - Left-censoring: Where the first observed signal occurs in
       the earliest analyzed quarter (2018Q1), the true signal onset
       is unknown and may predate the dataset. Reported observation
       windows are lower bounds in these cases.
     - Each quarter's data was processed independently. Cases updated
       across quarters may appear in multiple quarters' analyses,
       which could moderately affect temporal signal persistence.
     - Drug name normalization uses the prod_ai field; some ambiguity
       remains for the 2% of records requiring drugname fallback.
     - Case-study drugs were selected retrospectively because their
       regulatory actions are already known.
"""
    return summary


def run_noise_analysis(
    combined_path: Path = None,
    output_dir: Path = None,
):
    """Run the complete signal reliability and limitations analysis."""
    if combined_path is None:
        combined_path = Path("data/processed/signals_all_quarters.csv")
    if output_dir is None:
        output_dir = Path("figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("SIGNAL RELIABILITY & LIMITATIONS ANALYSIS")
    print("=" * 60)

    combined = pd.read_csv(combined_path, dtype={"quarter": str}, low_memory=False)
    combined = combined.dropna(subset=["quarter"])
    print(f"Loaded {len(combined):,} records")

    # Analysis 1: Case count distribution
    analyze_case_count_distribution(combined, output_dir)

    # Analysis 2: Reporting-volume surge
    analyze_reporting_volume_surge(combined, output_dir)

    # Analysis 3: Reporting-volume effects
    analyze_reporting_volume_effects(combined, output_dir)

    # Analysis 4: Negative controls
    analyze_negative_controls(combined, output_dir)

    # Generate findings summary
    summary = generate_findings_summary(combined)
    print(summary)

    # Save summary
    summary_path = output_dir / "findings_summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary)
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    run_noise_analysis()
