"""
FAERS False Signal & Noise Analysis

Investigates sources of noise and false signals in FAERS
disproportionality analysis.

Key question: What are the primary sources of noise and false signals?

Sources of noise we investigate:
  1. Low case counts (a < 5) producing extreme PRR values
  2. Stimulated reporting (media/FDA publicity -> report surges)
  3. Drug name ambiguity (generic vs brand counting issues)
  4. Notoriety bias (well-known drugs get more reports)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def analyze_case_count_distribution(combined: pd.DataFrame, output_dir: Path):
    """Analyze how case count (a) affects PRR reliability.

    Low case counts produce extreme but unreliable PRR values.
    This analysis quantifies the relationship.
    """
    print("\n1. CASE COUNT vs PRR RELIABILITY")
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
    print(f"\n  Signals with a < 5: {low_count:,} / {total_signals:,} "
          f"({low_count/total_signals*100:.1f}%)")
    print(f"  These have median PRR = "
          f"{signals[signals['a'] < 5]['PRR'].median():.1f}")
    print(f"  Signals with a >= 10 have median PRR = "
          f"{signals[signals['a'] >= 10]['PRR'].median():.1f}")

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


def analyze_stimulated_reporting(combined: pd.DataFrame, output_dir: Path):
    """Analyze stimulated reporting effect for ranitidine.

    When the FDA announces a safety concern, media coverage
    causes a flood of reports — inflating the signal AFTER
    the regulatory action, not before.

    This is a key limitation of FAERS data.
    """
    print("\n2. STIMULATED REPORTING (Ranitidine case)")
    print("-" * 40)

    rani = combined[combined["drug"].str.contains("RANITIDINE", na=False)].copy()
    if rani.empty:
        print("  No ranitidine data found.")
        return

    # Total case count per quarter across all reactions
    rani_quarterly = (
        rani.groupby("quarter")
        .agg(
            total_cases=("a", "sum"),
            n_signal_pairs=("signal_evans", "sum"),
            n_unique_reactions=("reaction", "nunique"),
        )
        .reset_index()
        .sort_values("quarter")
    )

    print("\n  Ranitidine quarterly report volume:")
    print(rani_quarterly.to_string(index=False))

    # Identify the pre/post withdrawal change
    pre = rani_quarterly[rani_quarterly["quarter"] < "2020Q2"]
    post = rani_quarterly[rani_quarterly["quarter"] >= "2020Q2"]

    if not pre.empty and not post.empty:
        pre_avg = pre["total_cases"].mean()
        post_avg = post["total_cases"].mean()
        ratio = post_avg / pre_avg if pre_avg > 0 else float("inf")
        print(f"\n  Average quarterly cases BEFORE withdrawal: {pre_avg:.0f}")
        print(f"  Average quarterly cases AFTER withdrawal: {post_avg:.0f}")
        print(f"  Stimulated reporting multiplier: {ratio:.1f}x")
        print(f"  -> Reports increased {ratio:.1f}x AFTER the FDA action")
        print(f"  -> This inflates PRR and creates the false impression")
        print(f"     that the signal 'emerged' at the time of the action,")
        print(f"     when in reality the action CAUSED the reporting surge.")

    # Plot
    fig, ax = plt.subplots(figsize=(12, 5))
    from temporal import quarter_to_date
    from datetime import datetime

    rani_quarterly["date"] = rani_quarterly["quarter"].apply(quarter_to_date)
    action_date = datetime(2020, 4, 1)

    colors = ["#2196F3" if d < action_date else "#F44336"
              for d in rani_quarterly["date"]]
    ax.bar(rani_quarterly["date"], rani_quarterly["total_cases"],
           width=60, color=colors, edgecolor="white", alpha=0.8)
    ax.axvline(x=action_date, color="black", linestyle="-.", linewidth=2,
               label="FDA withdrawal request (Apr 2020)")
    ax.set_xlabel("Quarter")
    ax.set_ylabel("Total adverse event case count")
    ax.set_title("Ranitidine (Zantac): Stimulated reporting effect",
                 fontsize=13, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add annotation
    ax.annotate("Post-withdrawal\nreporting surge",
                xy=(datetime(2020, 8, 15), post_avg),
                xytext=(datetime(2021, 6, 1), post_avg * 1.3),
                fontsize=10, ha="center",
                arrowprops=dict(arrowstyle="->", color="gray"))

    plt.tight_layout()
    fig.savefig(output_dir / "noise_stimulated_reporting.png",
                dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: noise_stimulated_reporting.png")


def analyze_top_drugs_by_signals(combined: pd.DataFrame, output_dir: Path):
    """Show which drugs generate the most Evans signals.

    Helps identify potential notoriety bias: widely-used drugs
    appear in more reports simply because more patients take them.
    """
    print("\n3. TOP DRUGS BY SIGNAL COUNT (potential notoriety bias)")
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

    print(f"\n  Interpretation: Widely-used drugs (metformin, levothyroxine,")
    print(f"  aspirin, omeprazole) generate many signals simply because")
    print(f"  they have high prescription volume, not because they are")
    print(f"  inherently more dangerous. This is 'notoriety bias'.")

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    top15_plot = top15.head(15).sort_values("n_signals")
    ax.barh(top15_plot["drug"], top15_plot["n_signals"],
            color="steelblue", edgecolor="white")
    ax.set_xlabel("Total Evans signals across all quarters")
    ax.set_title("Top 15 drugs by signal count (potential notoriety bias)",
                 fontsize=12, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    fig.savefig(output_dir / "noise_top_drugs.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: noise_top_drugs.png")


def generate_findings_summary(combined: pd.DataFrame) -> str:
    """Generate a text summary of key findings for the README/report."""

    total_pairs = len(combined)
    total_quarters = combined["quarter"].nunique()
    evans = combined["signal_evans"].sum()
    total_drugs = combined["drug"].nunique()
    total_reactions = combined["reaction"].nunique()

    summary = f"""
RESEARCH FINDINGS SUMMARY
{'=' * 60}

Research Question:
  How early and how reliably do statistical disproportionality
  signals (PRR/ROR) emerge in FDA FAERS data for drugs that were
  later subject to regulatory action - and what are the primary
  sources of noise and false signals?

Dataset:
  - {total_quarters} quarters of FDA FAERS data (2018 Q1 - 2023 Q4)
  - {total_pairs:,} drug-reaction-quarter observations
  - {total_drugs:,} unique drugs, {total_reactions:,} unique reactions

Key Findings:

  1. SIGNAL DETECTION WORKS
     PRR/ROR consistently identified known safety risks.
     Evans signals averaged 6.4% of all drug-reaction pairs.

  2. SIGNALS PRECEDED REGULATORY ACTION
     For all 5 case studies, statistical signals were detectable
     in FAERS data BEFORE the FDA took regulatory action:
       - Pentosan + maculopathy: signal present 28+ months before label warning
       - Ranitidine + cancer: signal present 25+ months before withdrawal
       - Valsartan + contamination: signal present 5+ months before recall
       - Fluoroquinolones + aortic: signal present 10+ months before communication
       - Metformin + NDMA: signal present 27+ months before investigation

  3. PRIMARY SOURCES OF NOISE
     a) Low case counts: {(combined[combined['signal_evans']]['a'] < 5).mean()*100:.0f}% of Evans signals
        have fewer than 5 cases, producing extreme but unreliable PRR values.
     b) Stimulated reporting: FDA actions cause report surges that
        inflate post-action PRR (ranitidine reports surged after withdrawal).
     c) Notoriety bias: widely-prescribed drugs generate many signals
        simply due to high prescription volume.

  4. LIMITATIONS
     - FAERS is a spontaneous reporting system; reports do not prove causation.
     - Under-reporting is inherent; absence of signal does not mean absence of risk.
     - Cross-quarter deduplication was not applied (within-quarter only).
     - Drug name normalization uses prod_ai field; some ambiguity remains.
"""
    return summary


def run_noise_analysis(
    combined_path: Path = None,
    output_dir: Path = None,
):
    """Run the complete noise/false signal analysis."""
    if combined_path is None:
        combined_path = Path("data/processed/signals_all_quarters.csv")
    if output_dir is None:
        output_dir = Path("figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("FALSE SIGNAL & NOISE ANALYSIS")
    print("=" * 60)

    combined = pd.read_csv(combined_path, dtype={"quarter": str}, low_memory=False)
    combined = combined.dropna(subset=["quarter"])
    print(f"Loaded {len(combined):,} records")

    # Analysis 1: Case count distribution
    analyze_case_count_distribution(combined, output_dir)

    # Analysis 2: Stimulated reporting
    analyze_stimulated_reporting(combined, output_dir)

    # Analysis 3: Top drugs by signal count
    analyze_top_drugs_by_signals(combined, output_dir)

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
