"""
Map SIDER STITCH IDs (CID#########) to ChEMBL IDs via PubChem CID and UniChem dump.

Pipeline:
  1. edges_drug_sideeffect_stitch.csv
       - stitch_id (e.g., 'CID000010917') -> pubchem_cid (int, e.g., 10917)
  2. UniChem src1src22.txt  (ChEMBL [1] -> PubChem [22])
       - chembl_id, pubchem_cid
  3. Invert UniChem mapping to get PubChem [22] -> ChEMBL [1]
       - pubchem_cid, chembl_id
  4. Join on pubchem_cid to create:
       sider_stitch_to_chembl_full.csv
         stitch_id, pubchem_cid, chembl_id, umls_cui
       sider_stitch_to_chembl_phase4.csv (optional filter)
"""

import re
import pandas as pd
from pathlib import Path

# ---- Paths  ----
BASE = Path(".")
PROCESSED = BASE / "data" / "processed"
UNICHEM_DIR = BASE / "unichem"

SIDER_EDGES = PROCESSED / "edges_drug_sideeffect_stitch.csv"
# This is the file WE DOWNLOADED FROM https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/: src1src22.txt  (ChEMBL -> PubChem)
UNICHEM_CHEMBL_TO_PUBCHEM = UNICHEM_DIR / "src1src22.txt"
CHEMBL_PHASE4 = PROCESSED / "nodes_drug_chembl_phase4.csv"  # optional filter

OUT_MAP_FULL = PROCESSED / "sider_stitch_to_chembl_full.csv"
OUT_MAP_PHASE4 = PROCESSED / "sider_stitch_to_chembl_phase4.csv"


def extract_pubchem_cid(stitch_id: str):
    """
    Convert STITCH-style SIDER ID (e.g., 'CID000010917') -> integer PubChem CID (10917).
    Returns None if the pattern doesn't match.
    """
    if not isinstance(stitch_id, str):
        return None
    sid = stitch_id.strip()
    # Pattern: 'CID' followed by digits, we strip leading zeros from the number part
    m = re.match(r"CID0*([1-9]\d*|0)$", sid)
    if not m:
        return None
    return int(m.group(1))


def main():
    # ---- 1) Load SIDER edges & derive PubChem CID ----
    print(f"Loading SIDER STITCH edges from: {SIDER_EDGES}")
    edges = pd.read_csv(SIDER_EDGES, dtype=str)
    edges["stitch_id"] = edges["stitch_id"].astype(str).str.strip()

    edges["pubchem_cid"] = edges["stitch_id"].apply(extract_pubchem_cid)
    before = len(edges)
    edges = edges.dropna(subset=["pubchem_cid"])
    edges["pubchem_cid"] = edges["pubchem_cid"].astype(int)

    print(f"  → Rows in SIDER edges: {before:,}")
    print(f"  → Rows with valid PubChem CID: {len(edges):,}")
    print(f"  → Unique PubChem CIDs in SIDER: {edges['pubchem_cid'].nunique():,}")

    # ---- 2) Load UniChem ChEMBL→PubChem and invert it ----
    print(f"\nLoading UniChem ChEMBL→PubChem mapping from: {UNICHEM_CHEMBL_TO_PUBCHEM}")
    if not UNICHEM_CHEMBL_TO_PUBCHEM.exists():
        raise FileNotFoundError(f"UniChem mapping file not found: {UNICHEM_CHEMBL_TO_PUBCHEM}")

    # The file looks like:
    # From src:'1'\tTo src:'22'
    # CHEMBL456354\t24949403
    # So we skip the first header line and treat the rest as two columns.
    uc_raw = pd.read_csv(
        UNICHEM_CHEMBL_TO_PUBCHEM,
        sep="\t",
        skiprows=1,
        header=None,
        names=["chembl_id", "pubchem_cid"],
        dtype=str,
    )

    print(f"  → Raw UniChem rows loaded: {len(uc_raw):,}")

    # Clean & convert
    uc_raw["pubchem_cid"] = pd.to_numeric(uc_raw["pubchem_cid"], errors="coerce")
    uc_raw = uc_raw.dropna(subset=["pubchem_cid"])
    uc_raw["pubchem_cid"] = uc_raw["pubchem_cid"].astype(int)
    uc_raw["chembl_id"] = uc_raw["chembl_id"].astype(str).str.strip()

    # Invert mapping: now each row is PubChem CID → ChEMBL
    uc = uc_raw[["pubchem_cid", "chembl_id"]].drop_duplicates()

    print(f"  → Inverted UniChem mapping rows: {len(uc):,}")
    print(f"  → Unique PubChem CIDs in UniChem: {uc['pubchem_cid'].nunique():,}")
    print(f"  → Unique ChEMBL IDs in UniChem: {uc['chembl_id'].nunique():,}")

    # ---- 3) Join SIDER PubChem → ChEMBL ----
    print("\nMerging SIDER PubChem CIDs with UniChem mapping…")
    merged = edges.merge(uc, on="pubchem_cid", how="inner")

    print(f"  → Mapped SIDER rows (any ChEMBL): {len(merged):,}")
    print(f"  → Unique SIDER STITCH IDs mapped: {merged['stitch_id'].nunique():,}")
    print(f"  → Unique mapped ChEMBL IDs: {merged['chembl_id'].nunique():,}")

    # Keep useful columns
    cols = ["stitch_id", "pubchem_cid", "chembl_id"]
    if "umls_cui" in merged.columns:
        cols.append("umls_cui")

    map_full = (
        merged[cols]
        .drop_duplicates()
        .sort_values(["stitch_id", "chembl_id"] + (["umls_cui"] if "umls_cui" in cols else []))
    )

    OUT_MAP_FULL.parent.mkdir(parents=True, exist_ok=True)
    map_full.to_csv(OUT_MAP_FULL, index=False)
    print(f"\n Full SIDER→PubChem→ChEMBL mapping saved to: {OUT_MAP_FULL}")

    # ---- 4) Filter to phase-4 ChEMBL drugs ----
    if CHEMBL_PHASE4.exists():
        print(f"\nFiltering mappings to ChEMBL phase-4 drugs using: {CHEMBL_PHASE4}")
        p4 = pd.read_csv(CHEMBL_PHASE4, dtype=str)
        p4["chembl_id"] = p4["chembl_id"].astype(str).str.strip()
        
        p4_ids = set(p4["chembl_id"].unique())
        print(f"  → Phase-4 ChEMBL IDs: {len(p4_ids):,}")

        map_phase4 = map_full[map_full["chembl_id"].isin(p4_ids)].copy()
        print(f"  → Rows mapped to phase-4 drugs: {len(map_phase4):,}")
        print(f"  → Unique SIDER STITCH IDs with phase-4 mapping: {map_phase4['stitch_id'].nunique():,}")
        print(f"  → Unique phase-4 ChEMBL IDs observed: {map_phase4['chembl_id'].nunique():,}")

        map_phase4.to_csv(OUT_MAP_PHASE4, index=False)
        print(f" Phase-4 filtered mapping saved to: {OUT_MAP_PHASE4}")
    else:
        print(f"\n Phase-4 file not found at {CHEMBL_PHASE4}, skipping phase-4 filter.")


if __name__ == "__main__":
    main()
