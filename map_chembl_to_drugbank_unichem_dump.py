"""
Offline ChEMBL → DrugBank mapping using UniChem dump (src1src2.txt).

Inputs:
  - data/processed/nodes_drug_chembl_phase4.csv
        (must have column: chembl_id)
  - src1src2.txt
        with columns: chembl_id    drugbank_id

Outputs:
  - data/processed/chembl_to_drugbank.csv
  - data/processed/chembl_unmapped.csv
"""

import pandas as pd
from pathlib import Path

# ---- Path Setup ----
BASE = Path(".")
PHASE4 = BASE / "data" / "processed" / "nodes_drug_chembl_phase4.csv"
UNICHEM_MAP = BASE / "src1src2.txt"
OUT_MAP  = BASE / "data" / "processed" / "chembl_to_drugbank.csv"
OUT_MISS = BASE / "data" / "processed" / "chembl_unmapped.csv"


def main():

    # -------------------------
    # 1. Load Phase-4 ChEMBL drugs
    # -------------------------
    print(f"Loading phase-4 ChEMBL drugs: {PHASE4}")
    df_drugs = pd.read_csv(PHASE4, dtype=str)
    df_drugs["chembl_id"] = df_drugs["chembl_id"].astype(str).str.strip()
    df_drugs = df_drugs.drop_duplicates(subset=["chembl_id"])
    print(f"  → {len(df_drugs):,} ChEMBL IDs")

    # -------------------------
    # 2. Load UniChem dump mapping
    # -------------------------
    print(f"Loading UniChem mapping: {UNICHEM_MAP}")
    df_map = pd.read_csv(UNICHEM_MAP, sep="\t", dtype=str)
    
    # Ensure expected columns exist
    if not {"chembl_id", "drugbank_id"}.issubset(df_map.columns):
        raise ValueError(f"ERROR: src1src2.txt must contain columns "
                         f"'chembl_id' and 'drugbank_id'. Found {df_map.columns}")

    df_map["chembl_id"] = df_map["chembl_id"].astype(str).str.strip()
    df_map["drugbank_id"] = df_map["drugbank_id"].astype(str).str.strip()

    print(f"  → {len(df_map):,} rows in UniChem mapping")

    # -------------------------
    # 3. Inner join: restrict mapping to our actual drug list
    # -------------------------
    merged = df_drugs.merge(df_map, on="chembl_id", how="inner")
    merged = merged.drop_duplicates(subset=["chembl_id", "drugbank_id"])

    print(f"  → {len(merged):,} mapped ChEMBL→DrugBank pairs")

    # -------------------------
    # 4. Identify unmapped ChEMBL IDs
    # -------------------------
    mapped_ids = set(merged["chembl_id"])
    all_ids = set(df_drugs["chembl_id"])
    missing_ids = sorted(all_ids - mapped_ids)

    # -------------------------
    # 5. Save results
    # -------------------------
    OUT_MAP.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_MAP, index=False)
    pd.DataFrame({"chembl_id": missing_ids}).to_csv(OUT_MISS, index=False)

    print("\n Mapping complete!")
    print(f"Total ChEMBL IDs        : {len(all_ids):,}")
    print(f"Mapped to DrugBank      : {len(mapped_ids):,}")
    print(f"Unmapped                : {len(missing_ids):,}")
    print(f"Mapping CSV saved to    : {OUT_MAP}")
    print(f"Unmapped list saved to  : {OUT_MISS}")

    # preview
    print("\nSample of mapped rows:")
    print(merged.head())


if __name__ == "__main__":
    main()
