// Drugs
CREATE CONSTRAINT drug_chembl_unique IF NOT EXISTS
FOR (d:Drug)
REQUIRE d.chembl_id IS UNIQUE;

// Targets
CREATE CONSTRAINT target_uniprot_unique IF NOT EXISTS
FOR (t:Target)
REQUIRE t.uniprot_id IS UNIQUE;

// Genes
CREATE CONSTRAINT gene_entrez_unique IF NOT EXISTS
FOR (g:Gene)
REQUIRE g.entrez_id IS UNIQUE;

// Diseases
CREATE CONSTRAINT disease_doid_unique IF NOT EXISTS
FOR (d:Disease)
REQUIRE d.doid IS UNIQUE;

// Side effects
CREATE CONSTRAINT sideeffect_umls_unique IF NOT EXISTS
FOR (s:SideEffect)
REQUIRE s.umls_cui IS UNIQUE;

// Compounds (Hetionet)
CREATE CONSTRAINT compound_drugbank_unique IF NOT EXISTS
FOR (c:Compound)
REQUIRE c.drugbank_id IS UNIQUE;

// Stitch bridge nodes
CREATE CONSTRAINT stitch_id_unique IF NOT EXISTS
FOR (s:Stitch)
REQUIRE s.stitch_id IS UNIQUE;
