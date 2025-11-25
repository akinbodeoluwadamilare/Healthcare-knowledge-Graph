# One-pass load of hetionet-v1.0.json.bz2 → Neo4j-ready CSVs
# - Normalizes Gene–Disease (associates) to Gene→Disease regardless of input order
# - Extracts evidence from multiple keys (source, sources, pmids), with compact JSON fallback

import bz2, json, pandas as pd
from pathlib import Path

# -------------------------
# Paths (auto-discover file)
# -------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROC_DIR = PROJECT_ROOT / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)

json_bz2 = RAW_DIR / "hetionet-v1.0.json.bz2"
if not json_bz2.exists():
    # try to find it anywhere under the project
    candidates = list(PROJECT_ROOT.rglob("hetionet-v1.0.json.bz2"))
    if candidates:
        json_bz2 = candidates[0]
    else:
        raise FileNotFoundError(
            f"Could not find {RAW_DIR/'hetionet-v1.0.json.bz2'}. "
            f"Place the file in ./data/raw/ or update this script's path."
        )

print(f"Using file: {json_bz2}")

# -------------------------
# Helpers
# -------------------------
def node_key(kind, identifier):
    return f"{kind}::{str(identifier)}"

def evidence_string(d):
    """Make a readable evidence string from Hetionet edge 'data'."""
    if not d:
        return ""
    if isinstance(d, dict):
        # Most concise, human-readable options first
        if "source" in d and isinstance(d["source"], str):
            return d["source"]
        if "sources" in d:
            s = d["sources"]
            if isinstance(s, (list, tuple)):
                return ";".join(sorted(map(str, s)))
            return str(s)
        if "pmids" in d:
            try:
                return f"pmids:{len(d['pmids'])}"
            except Exception:
                return "pmids"
        # Fallback: compact JSON (keeps provenance without bloating CSV)
        return json.dumps(d, ensure_ascii=False, separators=(",", ":"))
    return str(d)

# -------------------------
# Load JSON (one-pass)
# -------------------------
print("Loading JSON … (this can take a moment)")
with bz2.open(json_bz2, "rt", encoding="utf-8") as f:
    data = json.load(f)

nodes = data["nodes"]
edges = data["edges"]

# -------------------------
# Nodes (4 kinds)
# -------------------------
KEEP_KINDS = {"Compound", "Gene", "Disease", "Side Effect"}

compound_rows, gene_rows, disease_rows, se_rows = [], [], [], []

for n in nodes:
    k = n["kind"]
    if k not in KEEP_KINDS:
        continue
    ident = n["identifier"]
    name = n.get("name", "")

    if k == "Compound":
        # DrugBank IDs
        compound_rows.append({"drugbank_id": str(ident), "name": name})
    elif k == "Gene":
        # Entrez Gene IDs
        gene_rows.append({"entrez_id": str(ident), "symbol": name})
    elif k == "Disease":
        # DOID IDs
        disease_rows.append({"doid": str(ident), "name": name})
    elif k == "Side Effect":
        # UMLS CUI
        se_rows.append({"umls_cui": str(ident), "name": name})

pd.DataFrame(compound_rows).drop_duplicates().to_csv(PROC_DIR/"nodes_compound.csv", index=False)
pd.DataFrame(gene_rows).drop_duplicates().to_csv(PROC_DIR/"nodes_gene.csv", index=False)
pd.DataFrame(disease_rows).drop_duplicates().to_csv(PROC_DIR/"nodes_disease.csv", index=False)
pd.DataFrame(se_rows).drop_duplicates().to_csv(PROC_DIR/"nodes_sideeffect.csv", index=False)

# -------------------------
# Edges (4 relations)
# -------------------------
E_DRUG_TARGETS, E_DRUG_TREATS, E_GENE_DISEASE, E_DRUG_SE = [], [], [], []

for e in edges:
    ekind = (e.get("kind") or "").lower().strip()
    direc = (e.get("direction") or "both").lower().strip()

    (src_kind, src_id) = e["source_id"]
    (tgt_kind, tgt_id) = e["target_id"]
    src_kind, tgt_kind = str(src_kind), str(tgt_kind)

    ev = evidence_string(e.get("data"))

    # ---- Compound–binds–Gene (CbG). Prefer Compound→Gene; normalize if flipped when undirected ('both')
    if ekind == "binds" and direc in {"forward", "both"}:
        if src_kind == "Compound" and tgt_kind == "Gene":
            E_DRUG_TARGETS.append({"drugbank_id": str(src_id), "entrez_id": str(tgt_id), "evidence": ev})
        elif direc == "both" and src_kind == "Gene" and tgt_kind == "Compound":
            E_DRUG_TARGETS.append({"drugbank_id": str(tgt_id), "entrez_id": str(src_id), "evidence": ev})

    # ---- Compound–treats–Disease (CtD). Prefer Compound→Disease; normalize if flipped and 'both'
    elif ekind == "treats" and direc in {"forward", "both"}:
        if src_kind == "Compound" and tgt_kind == "Disease":
            E_DRUG_TREATS.append({"drugbank_id": str(src_id), "doid": str(tgt_id), "evidence": ev})
        elif direc == "both" and src_kind == "Disease" and tgt_kind == "Compound":
            E_DRUG_TREATS.append({"drugbank_id": str(tgt_id), "doid": str(src_id), "evidence": ev})

    # ---- Gene–associates–Disease (GaD). UNDIRECTED → normalize to Gene→Disease
    elif ekind == "associates" and direc in {"both", "forward"}:
        if src_kind == "Gene" and tgt_kind == "Disease":
            E_GENE_DISEASE.append({"entrez_id": str(src_id), "doid": str(tgt_id), "evidence": ev})
        elif src_kind == "Disease" and tgt_kind == "Gene":
            E_GENE_DISEASE.append({"entrez_id": str(tgt_id), "doid": str(src_id), "evidence": ev})

    # ---- Compound–causes–Side Effect (CcSE). Prefer Compound→Side Effect; normalize if flipped and 'both'
    elif ekind == "causes" and direc in {"forward", "both"}:
        if src_kind == "Compound" and tgt_kind == "Side Effect":
            E_DRUG_SE.append({"drugbank_id": str(src_id), "umls_cui": str(tgt_id), "evidence": ev})
        elif direc == "both" and src_kind == "Side Effect" and tgt_kind == "Compound":
            E_DRUG_SE.append({"drugbank_id": str(tgt_id), "umls_cui": str(src_id), "evidence": ev})

# Write edge CSVs
pd.DataFrame(E_DRUG_TARGETS).drop_duplicates().to_csv(PROC_DIR/"edges_drug_targets.csv", index=False)
pd.DataFrame(E_DRUG_TREATS).drop_duplicates().to_csv(PROC_DIR/"edges_drug_treats.csv", index=False)
pd.DataFrame(E_GENE_DISEASE).drop_duplicates().to_csv(PROC_DIR/"edges_gene_disease.csv", index=False)
pd.DataFrame(E_DRUG_SE).drop_duplicates().to_csv(PROC_DIR/"edges_drug_sideeffect.csv", index=False)

# -------------------------
# Quick sanity log
# -------------------------
print("Done. Files written to:", PROC_DIR.resolve())
for fn in [
    "nodes_compound.csv","nodes_gene.csv","nodes_disease.csv","nodes_sideeffect.csv",
    "edges_drug_targets.csv","edges_drug_treats.csv","edges_gene_disease.csv","edges_drug_sideeffect.csv"
]:
    p = PROC_DIR / fn
    try:
        n = sum(1 for _ in open(p, encoding="utf-8")) - 1  # minus header
        print(f"{fn:30s} rows: {n:,}")
    except Exception as ex:
        print(f"{fn:30s} (could not count): {ex}")
