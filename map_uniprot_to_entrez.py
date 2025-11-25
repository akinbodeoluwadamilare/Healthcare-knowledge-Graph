"""
Create UniProt → Entrez Gene mapping using UniProt ID mapping dump.

Inputs:
  - idmapping.dat.gz   (from UniProt FTP)
  - data/processed/nodes_target_uniprot.csv   (from your ChEMBL extraction)

Outputs:
  - data/processed/uniprot_to_entrez.csv
"""

import gzip
import csv
import pandas as pd
from pathlib import Path

BASE = Path(".")
PROCESSED = BASE / "data" / "processed"

UNIPROT_INPUT = BASE / "idmapping.dat.gz"
TARGETS = PROCESSED / "nodes_target_uniprot_phase4.csv"
OUT = PROCESSED / "uniprot_to_entrez.csv"


def main():
    print(f"Loading ChEMBL protein targets from: {TARGETS}")
    df = pd.read_csv(TARGETS, dtype=str)
    df["uniprot"] = df["uniprot"].str.strip()
    target_uniprots = set(df["uniprot"].dropna().unique())
    print(f"  → {len(target_uniprots):,} UniProt IDs needed")

    print(f"\nStreaming UniProt mapping file: {UNIPROT_INPUT}")
    mapping = []

    with gzip.open(UNIPROT_INPUT, "rt", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) != 3:
                continue

            uniprot, mapping_type, value = row

            # Only keep GeneID (Entrez)
            if mapping_type == "GeneID" and uniprot in target_uniprots:
                mapping.append((uniprot, value))

    print(f"  → Found {len(mapping):,} UniProt→Entrez mappings")

    # Reformat into DataFrame
    df_map = pd.DataFrame(mapping, columns=["uniprot", "entrez_id"])
    df_map.drop_duplicates(inplace=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    df_map.to_csv(OUT, index=False)

    print(f"\n Mapping saved to: {OUT}")
    print(df_map.head())


if __name__ == "__main__":
    main()
