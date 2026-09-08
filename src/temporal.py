"""
FAERS Temporal Analysis & Retrospective Case Study Comparison

Tracks signal strength (PRR/ROR) over time for specific drugs
and compares against known regulatory actions.

This module investigates the core research question:
  "How early and how reliably do disproportionality signals appear
  in FAERS data for drugs that were later subject to regulatory action?"

Note: This is a retrospective analysis. Case-study drugs were selected
because their regulatory actions are already known. The analysis does
not demonstrate prospective predictive capability.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for saving figures
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
from pathlib import Path
from datetime import datetime


# ── Case Study Drugs ─────────────────────────────────────────
# Each entry: (drug_normalized pattern, reaction pattern,
#              regulatory action date, action description)
#
# Drug names must match the 'drug_normalized' column (uppercase active ingredient).
# Reaction names must match the 'reaction' column (MedDRA Preferred Term).

CASE_STUDIES = [
    {
        "name": "Ranitidine (Zantac) - NDMA contamination",
        "drug_pattern": "RANITIDINE",
        "reactions": ["Gastrointestinal carcinoma", "Bladder cancer",
                      "Gastric cancer", "Hepatocellular carcinoma",
                      "Colorectal cancer", "Oesophageal carcinoma",
                      "Neoplasm malignant"],
        "action_date": "2020-04-01",
        "action_desc": "FDA requests market withdrawal",
    },
    {
        "name": "Valsartan - NDMA contamination",
        "drug_pattern": "VALSARTAN",
        "reactions": ["Neoplasm malignant", "Product contamination",
                      "Product quality issue", "Hepatocellular carcinoma",
                      "Bladder cancer", "Gastric cancer",
                      "Colorectal cancer"],
        "action_date": "2018-07-13",
        "action_desc": "FDA announces voluntary recall",
    },
    {
        "name": "Pentosan polysulfate (Elmiron) - Maculopathy",
        "drug_pattern": "PENTOSAN POLYSULFATE SODIUM",
        "reactions": ["Macular degeneration", "Maculopathy",
                      "Pigmentary maculopathy", "Retinal pigmentation",
                      "Retinal pigment epitheliopathy",
                      "Visual impairment", "Vision blurred"],
        "action_date": "2020-06-16",
        "action_desc": "FDA adds label warning for retinal damage",
    },
    {
        "name": "Fluoroquinolones - Aortic dissection",
        "drug_pattern": "LEVOFLOXACIN",
        "reactions": ["Aortic dissection", "Aortic aneurysm",
                      "Aortic rupture", "Tendon rupture",
                      "Peripheral neuropathy"],
        "action_date": "2018-12-20",
        "action_desc": "FDA safety communication on aortic risks",
    },
    {
        "name": "Metformin - NDMA impurity",
        "drug_pattern": "METFORMIN",
        "reactions": ["Product contamination", "Product quality issue",
                      "Neoplasm malignant", "Bladder cancer",
                      "Lactic acidosis"],
        "action_date": "2020-05-28",
        "action_desc": "FDA investigating NDMA in metformin",
    },
]


def quarter_to_date(quarter_str: str) -> datetime:
    """Convert quarter string like '2024Q1' to a datetime.

    Maps to the middle of the quarter for plotting.
    """
    q_str = quarter_str.upper().strip()
    year = int(q_str[:4])
    q = int(q_str[-1])
    # Middle of each quarter
    month = {1: 2, 2: 5, 3: 8, 4: 11}[q]
    return datetime(year, month, 15)


def extract_drug_signals(
    combined: pd.DataFrame,
    drug_pattern: str,
    reactions: list[str] | None = None,
    top_n: int = 5,
) -> pd.DataFrame:
    """Extract quarterly signal data for a specific drug.

    Args:
        combined: Combined signals DataFrame from all quarters.
        drug_pattern: Drug name to match (exact match on drug_normalized).
        reactions: Optional list of specific reactions to track.
                   If None, finds the top N reactions by signal strength.
        top_n: Number of top reactions to track if reactions is None.

    Returns:
        DataFrame with quarterly PRR/ROR for the drug's top reactions.
    """
    # Filter to this drug
    drug_data = combined[combined["drug"] == drug_pattern].copy()

    if drug_data.empty:
        # Try partial match
        drug_data = combined[combined["drug"].str.contains(drug_pattern, na=False)].copy()

    if drug_data.empty:
        print(f"  WARNING: No data found for drug pattern '{drug_pattern}'")
        return pd.DataFrame()

    # If no specific reactions given, find the top N by average PRR
    if reactions is None:
        top_reactions = (
            drug_data[drug_data["a"] >= 3]
            .groupby("reaction")["PRR"]
            .mean()
            .nlargest(top_n)
            .index.tolist()
        )
    else:
        # Match case-insensitively
        reactions_upper = [r.lower() for r in reactions]
        available = drug_data["reaction"].str.lower().unique()
        top_reactions = [
            drug_data[drug_data["reaction"].str.lower() == r]["reaction"].iloc[0]
            for r in reactions_upper
            if r in available
        ]
        if not top_reactions:
            # Fallback to top N
            top_reactions = (
                drug_data[drug_data["a"] >= 3]
                .groupby("reaction")["PRR"]
                .mean()
                .nlargest(top_n)
                .index.tolist()
            )

    # Filter to top reactions
    result = drug_data[drug_data["reaction"].isin(top_reactions)].copy()

    # Add date column
    result["date"] = result["quarter"].apply(quarter_to_date)

    return result.sort_values("date")


def plot_signal_timeline(
    drug_signals: pd.DataFrame,
    case_study: dict,
    output_dir: Path,
) -> Path | None:
    """Plot PRR over time for a case-study drug.

    Shows:
      - PRR trajectory for each tracked reaction
      - Evans signal threshold (PRR=2)
      - Regulatory action date (vertical line)
      - 95% confidence intervals

    Args:
        drug_signals: Quarterly signal data for this drug.
        case_study: Case study dictionary with metadata.
        output_dir: Directory to save the figure.

    Returns:
        Path to the saved figure, or None if no data.
    """
    if drug_signals.empty:
        return None

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)

    reactions = drug_signals["reaction"].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, max(len(reactions), 3)))

    action_date = datetime.strptime(case_study["action_date"], "%Y-%m-%d")

    # ── Top panel: PRR ──
    ax = axes[0]
    for i, reaction in enumerate(reactions):
        rdata = drug_signals[drug_signals["reaction"] == reaction].sort_values("date")
        if len(rdata) < 2:
            continue

        ax.plot(rdata["date"], rdata["PRR"], "o-",
                color=colors[i], label=reaction, markersize=4, linewidth=1.5)

        # Confidence interval band
        ax.fill_between(rdata["date"], rdata["PRR_lower"], rdata["PRR_upper"],
                        alpha=0.15, color=colors[i])

    ax.axhline(y=2, color="red", linestyle="--", alpha=0.5, label="Evans threshold (PRR=2)")
    ax.axvline(x=action_date, color="black", linestyle="-.", alpha=0.7,
               label=f"Regulatory action ({case_study['action_date'][:7]})")

    ax.set_ylabel("PRR (log scale)")
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.1)
    ax.legend(loc="upper left", fontsize=7, ncol=2)
    ax.set_title(case_study["name"], fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3)

    # ── Bottom panel: Case count (a) ──
    ax2 = axes[1]
    for i, reaction in enumerate(reactions):
        rdata = drug_signals[drug_signals["reaction"] == reaction].sort_values("date")
        if len(rdata) < 2:
            continue
        ax2.bar(rdata["date"], rdata["a"], width=60, alpha=0.6,
                color=colors[i], label=reaction)

    ax2.axvline(x=action_date, color="black", linestyle="-.", alpha=0.7)
    ax2.set_ylabel("Case count (a)")
    ax2.set_xlabel("Quarter")
    ax2.legend(loc="upper left", fontsize=7, ncol=2)
    ax2.grid(True, alpha=0.3)

    # Format x-axis
    ax2.xaxis.set_major_locator(mdates.YearLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax2.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))

    plt.tight_layout()

    # Save
    safe_name = case_study["name"].split(" - ")[0].lower().replace(" ", "_").replace("(", "").replace(")", "")
    fig_path = output_dir / f"timeline_{safe_name}.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"    Saved: {fig_path.name}")
    return fig_path


def generate_summary_table(
    combined: pd.DataFrame,
    case_studies: list[dict],
) -> pd.DataFrame:
    """Generate a summary table showing signal emergence for each case study.

    For each drug, reports:
      - First quarter where Evans signal was detected for a RELEVANT reaction
      - Regulatory action date
      - Observation window (quarters between first observed signal and action)
      - Peak PRR value among relevant reactions

    IMPORTANT: "First signal quarter" is restricted to the clinically relevant
    reactions defined for each case study, NOT all reactions for the drug.
    The temporal unit is the FAERS file quarter (which quarter's data file
    the case appears in), not the individual event date.
    """
    rows = []

    for cs in case_studies:
        drug_data = combined[
            combined["drug"].str.contains(cs["drug_pattern"], na=False)
        ]

        if drug_data.empty:
            rows.append({
                "Drug": cs["name"].split(" - ")[0],
                "Action": cs["action_desc"],
                "Action Date": cs["action_date"][:7],
                "First Signal Quarter": "No data",
                "Observation Window": "N/A",
                "Peak PRR": "N/A",
                "Peak Reaction": "N/A",
            })
            continue

        # Restrict to case-study-specific reactions if defined
        if cs.get("reactions"):
            reactions_lower = [r.lower() for r in cs["reactions"]]
            drug_data = drug_data[
                drug_data["reaction"].str.lower().isin(reactions_lower)
            ]

        if drug_data.empty:
            rows.append({
                "Drug": cs["name"].split(" - ")[0],
                "Action": cs["action_desc"],
                "Action Date": cs["action_date"][:7],
                "First Signal Quarter": "No relevant reaction data",
                "Observation Window": "N/A",
                "Peak PRR": "N/A",
                "Peak Reaction": "N/A",
            })
            continue

        # Find Evans signals (a >= 3) among RELEVANT reactions only
        signals = drug_data[(drug_data["signal_evans"]) & (drug_data["a"] >= 3)]

        if signals.empty:
            first_q = "No Evans signal for relevant reactions"
            lead = "N/A"
        else:
            signals = signals.copy()
            signals["date"] = signals["quarter"].apply(quarter_to_date)
            first_q = signals.sort_values("date")["quarter"].iloc[0]
            first_date = quarter_to_date(first_q)
            action_date = datetime.strptime(cs["action_date"], "%Y-%m-%d")
            lead_months = (action_date - first_date).days / 30.44
            if lead_months > 0:
                lead = f"{lead_months:.0f} months before action"
            else:
                lead = f"{abs(lead_months):.0f} months after action"

        # Peak PRR among relevant reactions (with a >= 3)
        peaks = drug_data[drug_data["a"] >= 3]
        if not peaks.empty:
            peak_idx = peaks["PRR"].idxmax()
            peak_prr = peaks.loc[peak_idx, "PRR"]
            peak_reaction = peaks.loc[peak_idx, "reaction"]
        else:
            peak_prr = "N/A"
            peak_reaction = "N/A"

        rows.append({
            "Drug": cs["name"].split(" - ")[0],
            "Action": cs["action_desc"],
            "Action Date": cs["action_date"][:7],
            "First Signal Quarter": first_q,
            "Observation Window": lead,
            "Peak PRR": f"{peak_prr:,.1f}" if isinstance(peak_prr, (int, float)) else peak_prr,
            "Peak Reaction": peak_reaction,
        })

    return pd.DataFrame(rows)


def run_temporal_analysis(
    combined_path: Path = None,
    output_dir: Path = None,
):
    """Run the complete temporal analysis and case study validation.

    Args:
        combined_path: Path to the combined signals CSV.
        output_dir: Directory for figures and results.
    """
    if combined_path is None:
        combined_path = Path("data/processed/signals_all_quarters.csv")
    if output_dir is None:
        output_dir = Path("figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("TEMPORAL ANALYSIS & RETROSPECTIVE CASE STUDY COMPARISON")
    print("=" * 60)

    # Load combined data
    print(f"\nLoading combined signals from {combined_path.name}...")
    combined = pd.read_csv(combined_path, dtype={"quarter": str}, low_memory=False)
    # Drop rows with missing quarter (shouldn't happen, but defensive)
    combined = combined.dropna(subset=["quarter"])
    combined["quarter"] = combined["quarter"].astype(str).str.strip().str.upper()
    print(f"  {len(combined):,} drug-reaction-quarter records")
    quarters_sorted = sorted(combined["quarter"].unique())
    print(f"  Quarters: {quarters_sorted}")

    # ── Case Study Analysis ──
    print("\n" + "-" * 60)
    print("CASE STUDY SIGNAL TIMELINES")
    print("-" * 60)

    for cs in CASE_STUDIES:
        print(f"\n  {cs['name']}")
        print(f"  Drug pattern: {cs['drug_pattern']}")

        drug_signals = extract_drug_signals(
            combined, cs["drug_pattern"], cs["reactions"]
        )

        if drug_signals.empty:
            print(f"  -> No data found")
            continue

        n_quarters = drug_signals["quarter"].nunique()
        n_reactions = drug_signals["reaction"].nunique()
        print(f"  -> {n_quarters} quarters, {n_reactions} reactions tracked")

        plot_signal_timeline(drug_signals, cs, output_dir)

    # ── Summary Table ──
    print("\n" + "-" * 60)
    print("SIGNAL EMERGENCE SUMMARY")
    print("-" * 60)

    summary = generate_summary_table(combined, CASE_STUDIES)
    print()
    print(summary.to_string(index=False))

    # Save summary
    summary_path = output_dir / "case_study_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"\n  Saved: {summary_path.name}")

    # ── Overview Statistics ──
    print("\n" + "-" * 60)
    print("OVERALL SIGNAL STATISTICS ACROSS QUARTERS")
    print("-" * 60)

    quarter_stats = (
        combined
        .groupby("quarter")
        .agg(
            total_pairs=("drug", "size"),
            evans_signals=("signal_evans", "sum"),
            median_prr=("PRR", "median"),
            pairs_a_gte_3=("a", lambda x: (x >= 3).sum()),
        )
        .reset_index()
    )
    quarter_stats["signal_rate"] = (
        quarter_stats["evans_signals"] / quarter_stats["total_pairs"] * 100
    )
    quarter_stats = quarter_stats.sort_values("quarter")
    print()
    print(quarter_stats.to_string(index=False))

    print("\n" + "=" * 60)
    print("TEMPORAL ANALYSIS COMPLETE")
    print(f"  Figures saved to: {output_dir}/")
    print("=" * 60)


# -- CLI entry point --
if __name__ == "__main__":
    run_temporal_analysis()
