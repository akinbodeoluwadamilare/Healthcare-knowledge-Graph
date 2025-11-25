# extract_sider.py
# Extracts drug–side effect associations from SIDER TSVs
# Compatible with SIDER 4.1+ (EMBL)

import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR
RAW = PROJECT_ROOT / "data" / "raw"
OUT = PROJECT_ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

# Column names differ slightly across releases; these are safe defaults
all_se_cols = ["stitch_id","umls_cui","side_effect_name","method","meddra_type"]
freq_cols   = ["stitch_id","umls_cui","frequency_lower","frequency_upper","meddra_type","sample_size"]

# 1) Read SIDER full associations
all_se = pd.read_csv(RAW/"meddra_all_se.tsv", sep="\t", header=None, names=all_se_cols, low_memory=False)
# Convert IDs to strings to ensure consistent merge keys
all_se["stitch_id"] = all_se["stitch_id"].astype(str)
all_se["umls_cui"]  = all_se["umls_cui"].astype(str)

edges = all_se[["stitch_id","umls_cui"]].drop_duplicates()

# 2) Read frequency file (optional)
freq = pd.read_csv(RAW/"meddra_freq.tsv", sep="\t", header=None, names=freq_cols, low_memory=False)
freq["stitch_id"] = freq["stitch_id"].astype(str)
freq["umls_cui"]  = freq["umls_cui"].astype(str)

# Keep only relevant columns and merge on string keys
freq = freq[["stitch_id","umls_cui","frequency_lower","frequency_upper"]].drop_duplicates()

edges_freq = edges.merge(freq, on=["stitch_id","umls_cui"], how="left")

# 3) Side effect nodes (UMLS CUI, name)
se_nodes = all_se[["umls_cui","side_effect_name"]].drop_duplicates().rename(
    columns={"side_effect_name":"name"}
)

# 4) Save results
se_nodes.to_csv(OUT/"nodes_side_effect_sider.csv", index=False)
edges_freq.to_csv(OUT/"edges_drug_sideeffect_stitch.csv", index=False)

print("Wrote:")
print(" -", (OUT/"nodes_side_effect_sider.csv").resolve())
print(" -", (OUT/"edges_drug_sideeffect_stitch.csv").resolve())
print(f"Nodes: {len(se_nodes):,} | Edges: {len(edges_freq):,}")
