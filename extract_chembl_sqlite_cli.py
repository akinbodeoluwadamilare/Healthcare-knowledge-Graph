"""
ChEMBL 36 → Neo4j-ready CSVs with phase filter and name fill.

Outputs (suffix varies by phase):
  nodes_drug_chembl_<suffix>.csv
  nodes_target_uniprot_<suffix>.csv
  edges_drug_targets_chembl_<suffix>.csv

Phase options:
  --phase all   : no filter (everything)
  --phase 4     : approved drugs (max_phase = 4)
  --phase 1-4   : any clinical phase (1..4)

Names:
  name = COALESCE(pref_name, synonym, chembl_id)
"""

import argparse
import sqlite3
import pandas as pd
from pathlib import Path

def get_cols(con, table):
    return set(pd.read_sql(f"PRAGMA table_info({table});", con)["name"].tolist())

def table_exists(con, table):
    q = "SELECT name FROM sqlite_master WHERE type='table' AND name=?"
    return pd.read_sql(q, con, params=[table]).shape[0] == 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="Path to chembl_36.db")
    ap.add_argument("--outdir", default="data/processed", help="Output folder")
    ap.add_argument("--phase", choices=["all","4","1-4"], default="4",
                    help="Phase filter (default: 4)")
    ap.add_argument("--fill-synonyms", action="store_true", default=True,
                    help="Fill missing names from molecule_synonyms (default: on)")
    ap.add_argument("--no-fill-synonyms", dest="fill_synonyms", action="store_false")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    suffix = {"all":"all", "4":"phase4", "1-4":"phase1-4"}[args.phase]

    con = sqlite3.connect(args.db)

    # ---------- schema introspection ----------
    td_cols = get_cols(con, "target_dictionary")
    dm_cols = get_cols(con, "drug_mechanism")
    md_cols = get_cols(con, "molecule_dictionary")

    # dm -> td (prefer numeric tid)
    if "tid" in dm_cols and "tid" in td_cols:
        dm_td_join = "td.tid = dm.tid"
    elif "target_chembl_id" in dm_cols and ("chembl_id" in td_cols or "target_chembl_id" in td_cols):
        td_id = "chembl_id" if "chembl_id" in td_cols else "target_chembl_id"
        dm_td_join = f"td.{td_id} = dm.target_chembl_id"
    else:
        raise RuntimeError("Cannot determine join between drug_mechanism and target_dictionary")

    # dm -> md (get CHEMBL ID from molregno)
    if "molregno" in dm_cols and "molregno" in md_cols:
        dm_md_join = "md.molregno = dm.molregno"
    elif "molecule_chembl_id" in dm_cols and "chembl_id" in md_cols:
        dm_md_join = "md.chembl_id = dm.molecule_chembl_id"
    else:
        raise RuntimeError("Cannot determine join between drug_mechanism and molecule_dictionary")

    # ---------- synonym helper table ----------
    syn_clause = ""
    name_expr = "COALESCE(md.pref_name, md.chembl_id) AS name"  # ensure not null
    if args.fill_synonyms and table_exists(con, "molecule_synonyms"):
        # one deterministic synonym per molregno
        q_syn = """
        SELECT ms.molregno, MIN(ms.synonyms) AS syn
        FROM molecule_synonyms ms
        GROUP BY ms.molregno
        """
        pd.read_sql(q_syn, con).to_sql("_tmp_syn", con, if_exists="replace", index=False)
        syn_clause = "LEFT JOIN _tmp_syn s ON s.molregno = md.molregno"
        name_expr = "COALESCE(md.pref_name, s.syn, md.chembl_id) AS name"

    # ---------- phase WHERE fragment ----------
    # We apply the filter via molecule_dictionary (md)
    where_phase = "1=1"
    if args.phase == "4":
        where_phase = "md.max_phase = 4"
    elif args.phase == "1-4":
        # keep any clinical phase (>=1), cap at 4 just in case
        where_phase = "md.max_phase >= 1 AND md.max_phase <= 4"

    # ---------- 1) drug nodes ----------
    q_drugs = f"""
    SELECT md.chembl_id AS chembl_id,
           {name_expr}
    FROM molecule_dictionary md
    {syn_clause}
    WHERE md.chembl_id IS NOT NULL
      AND {where_phase}
    """
    drugs = pd.read_sql(q_drugs, con).drop_duplicates()
    drugs.to_csv(outdir / f"nodes_drug_chembl_{suffix}.csv", index=False)

    # ---------- 2) drug→target edges ----------
    q_edges = f"""
    SELECT DISTINCT
      md.chembl_id            AS chembl_id,
      cs.accession            AS uniprot,
      dm.mechanism_of_action  AS mechanism,
      dm.action_type          AS action_type
    FROM drug_mechanism dm
    JOIN molecule_dictionary md
      ON {dm_md_join}
    JOIN target_dictionary td
      ON {dm_td_join}
    JOIN target_components tc
      ON tc.tid = td.tid
    JOIN component_sequences cs
      ON cs.component_id = tc.component_id
    WHERE md.chembl_id IS NOT NULL
      AND cs.accession IS NOT NULL
      AND {where_phase}
    """
    edges = pd.read_sql(q_edges, con).drop_duplicates()
    edges.to_csv(outdir / f"edges_drug_targets_chembl_{suffix}.csv", index=False)

    # ---------- 3) target nodes (only those referenced by edges) ----------
    # Build from edges to keep file small & consistent with phase
    targets = (
        edges[["uniprot"]]
        .drop_duplicates()
        .merge(
            pd.read_sql(f"""
                SELECT DISTINCT cs.accession AS uniprot,
                                td.pref_name AS target_name
                FROM target_dictionary td
                JOIN target_components tc ON tc.tid = td.tid
                JOIN component_sequences cs ON cs.component_id = tc.component_id
            """, con),
            on="uniprot", how="left"
        )
        .fillna({"target_name": ""})
    )
    targets.to_csv(outdir / f"nodes_target_uniprot_{suffix}.csv", index=False)

    con.close()

    print(" Done.")
    print("Rows:")
    print(" - drugs   :", len(drugs))
    print(" - edges   :", len(edges))
    print(" - targets :", len(targets))
    print("Wrote:")
    print(" -", (outdir / f"nodes_drug_chembl_{suffix}.csv").resolve())
    print(" -", (outdir / f"edges_drug_targets_chembl_{suffix}.csv").resolve())
    print(" -", (outdir / f"nodes_target_uniprot_{suffix}.csv").resolve())

if __name__ == "__main__":
    main()
