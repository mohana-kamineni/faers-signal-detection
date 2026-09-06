"""
FAERS Data Cleaning Module

Handles the three biggest data quality problems in FAERS:
1. Duplicate reports (same caseid, multiple primaryids)
2. Deleted/retracted cases
3. Drug name inconsistencies

After cleaning, each caseid appears exactly once (the most recent version),
deleted cases are excluded, and drug names are normalized.
"""

import pandas as pd
from pathlib import Path


def deduplicate_demo(demo: pd.DataFrame) -> pd.DataFrame:
    """Remove duplicate case reports from the DEMO table.

    FAERS contains duplicate reports: the same adverse event (same caseid)
    reported multiple times with different primaryids. Each report may have
    a different caseversion.

    Strategy: For each caseid, keep only the row with the HIGHEST caseversion
    (the most recent update). This is the standard deduplication method used
    in pharmacovigilance research.

    Args:
        demo: Raw DEMO DataFrame.

    Returns:
        Deduplicated DEMO DataFrame (one row per caseid).
    """
    before = len(demo)

    # Sort by caseid and caseversion, then keep the last (highest version)
    demo_sorted = demo.sort_values(["caseid", "caseversion"], ascending=True)
    demo_dedup = demo_sorted.drop_duplicates(subset=["caseid"], keep="last")

    after = len(demo_dedup)
    removed = before - after
    print(f"  Deduplication: {before:,} -> {after:,} reports "
          f"({removed:,} duplicates removed, {removed/before*100:.1f}%)")

    return demo_dedup.reset_index(drop=True)


def remove_deleted_cases(
    tables: dict[str, pd.DataFrame],
    delete_ids: set[int],
) -> dict[str, pd.DataFrame]:
    """Remove deleted/retracted cases from all tables.

    The FDA publishes a DELETE file listing caseids that have been
    retracted from the database. We exclude these from analysis.

    Args:
        tables: Dictionary of DataFrames (keys: 'demo', 'drug', 'reac', 'outc').
        delete_ids: Set of caseid integers to exclude.

    Returns:
        Updated tables dictionary with deleted cases removed.
    """
    if not delete_ids:
        print("  No deleted cases to remove.")
        return tables

    cleaned = {}
    for name, df in tables.items():
        before = len(df)
        # Convert caseid to numeric for comparison
        mask = ~df["caseid"].astype(str).isin({str(d) for d in delete_ids})
        cleaned[name] = df[mask].reset_index(drop=True)
        after = len(cleaned[name])
        removed = before - after
        if removed > 0:
            print(f"  Removed {removed:,} rows from {name.upper()} (deleted cases)")

    return cleaned


def normalize_drug_names(drug: pd.DataFrame) -> pd.DataFrame:
    """Normalize drug names to a single standardized identifier.

    Problem: The 'drugname' field contains brand names, generic names,
    abbreviations, and typos. The same drug can appear as:
      "LIPITOR", "Lipitor", "ATORVASTATIN CALCIUM", "atorvastatin", etc.

    Strategy:
      1. Use 'prod_ai' (active ingredient) as primary identifier — it's
         more standardized than 'drugname'.
      2. Where 'prod_ai' is missing, fall back to uppercase-stripped 'drugname'.
      3. Store the normalized name in a new column 'drug_normalized'.

    Args:
        drug: DRUG DataFrame.

    Returns:
        DRUG DataFrame with added 'drug_normalized' column.
    """
    # Clean prod_ai: uppercase, strip whitespace
    drug["prod_ai_clean"] = (
        drug["prod_ai"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    # Clean drugname: uppercase, strip whitespace
    drug["drugname_clean"] = (
        drug["drugname"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    # Use prod_ai_clean when available, otherwise drugname_clean
    drug["drug_normalized"] = drug["prod_ai_clean"].where(
        drug["prod_ai_clean"] != "",
        drug["drugname_clean"],
    )

    # Report stats
    total = len(drug)
    has_prod_ai = (drug["prod_ai_clean"] != "").sum()
    used_drugname = total - has_prod_ai
    unique_before = drug["drugname"].nunique()
    unique_after = drug["drug_normalized"].nunique()

    print(f"  Drug normalization:")
    print(f"    {has_prod_ai:,} rows used prod_ai ({has_prod_ai/total*100:.1f}%)")
    print(f"    {used_drugname:,} rows fell back to drugname ({used_drugname/total*100:.1f}%)")
    print(f"    Unique drug names: {unique_before:,} -> {unique_after:,} after normalization")

    # Drop intermediate columns
    drug = drug.drop(columns=["prod_ai_clean", "drugname_clean"])

    return drug


def filter_suspect_drugs(drug: pd.DataFrame) -> pd.DataFrame:
    """Filter to only Primary Suspect (PS) drugs.

    In FAERS, each drug in a case has a role:
      PS = Primary Suspect (most likely caused the adverse event)
      SS = Secondary Suspect
      C  = Concomitant (taken at the same time but not suspected)
      I  = Interacting

    For disproportionality analysis, we typically analyze only PS drugs.
    This is because including concomitant drugs would dilute the signal.

    Args:
        drug: DRUG DataFrame.

    Returns:
        Filtered DataFrame containing only Primary Suspect drugs.
    """
    before = len(drug)
    drug_ps = drug[drug["role_cod"] == "PS"].reset_index(drop=True)
    after = len(drug_ps)
    print(f"  Filtered to Primary Suspect drugs: {before:,} -> {after:,} rows")

    # Show role distribution for context
    role_counts = drug["role_cod"].value_counts()
    print(f"    Role distribution: {dict(role_counts)}")

    return drug_ps


def clean_quarter(
    tables: dict[str, pd.DataFrame],
    delete_ids: set[int],
) -> dict[str, pd.DataFrame]:
    """Apply all cleaning steps to a quarter's data.

    This is the main function that chains together all cleaning steps.

    Args:
        tables: Dictionary from load_quarter() with keys 'demo', 'drug', 'reac', 'outc'.
        delete_ids: Set of deleted caseid integers.

    Returns:
        Cleaned tables dictionary.
    """
    print("\n" + "=" * 60)
    print("DATA CLEANING")
    print("=" * 60)

    # Step 1: Deduplicate DEMO
    print("\n1. Deduplicating reports...")
    tables["demo"] = deduplicate_demo(tables["demo"])

    # Get the set of valid primaryids after deduplication
    valid_primaryids = set(tables["demo"]["primaryid"])
    print(f"   Valid primaryids after dedup: {len(valid_primaryids):,}")

    # Step 2: Remove deleted cases
    print("\n2. Removing deleted cases...")
    tables = remove_deleted_cases(tables, delete_ids)

    # Step 3: Filter DRUG, REAC, OUTC to only valid primaryids
    print("\n3. Aligning tables to deduplicated DEMO...")
    for name in ["drug", "reac", "outc"]:
        before = len(tables[name])
        tables[name] = tables[name][
            tables[name]["primaryid"].isin(valid_primaryids)
        ].reset_index(drop=True)
        after = len(tables[name])
        print(f"   {name.upper()}: {before:,} -> {after:,} rows")

    # Step 4: Normalize drug names
    print("\n4. Normalizing drug names...")
    tables["drug"] = normalize_drug_names(tables["drug"])

    # Step 5: Filter to Primary Suspect drugs
    print("\n5. Filtering to suspect drugs...")
    tables["drug_ps"] = filter_suspect_drugs(tables["drug"])

    # Final summary
    print("\n" + "=" * 60)
    print("CLEANING COMPLETE")
    print(f"  Unique cases: {tables['demo']['caseid'].nunique():,}")
    print(f"  Drug-reaction pairs possible: "
          f"{len(tables['drug_ps']):,} drugs × {len(tables['reac']):,} reactions")
    print(f"  Unique normalized drugs (PS only): "
          f"{tables['drug_ps']['drug_normalized'].nunique():,}")
    print(f"  Unique reactions: {tables['reac']['pt'].nunique():,}")
    print("=" * 60)

    return tables


# ── CLI entry point ──────────────────────────────────────────
if __name__ == "__main__":
    from ingest import load_quarter, load_delete_list

    extract_dir = Path("data/raw/faers_ascii_2024q1")

    # Load raw data
    tables = load_quarter(extract_dir)
    delete_ids = load_delete_list(extract_dir)

    # Clean
    tables = clean_quarter(tables, delete_ids)
