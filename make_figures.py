import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path("figures")
OUT_DIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# Fig 1 – Top 20 Degree Genes (Horizontal Bar Plot)
def plot_top_genes_degree():
    df = pd.read_csv("results/top_genes_degree.csv")
    # Sort ascending so highest degree is at top of horizontal bar chart
    df = df.sort_values("degree")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.barh(df["symbol"], df["degree"])
    ax.set_xlabel("Degree (number of associations)")
    ax.set_ylabel("Gene")
    ax.set_title("Top 20 Genes by Degree in Disease–Gene Network", fontsize=16, pad=20)
    plt.tight_layout(pad=2.5)
    fig.savefig(OUT_DIR / "fig1_top_genes_degree.png", dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / "fig1_top_genes_degree.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)

# Fig 2 – Top 20 Degree Diseases
def plot_top_diseases_degree():
    df = pd.read_csv("results/top_diseases_degree.csv")
    df = df.sort_values("degree")

    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.barh(df["disease_name"], df["degree"])
    ax.set_xlabel("Degree (number of associated genes)")
    ax.set_ylabel("Disease")
    ax.set_title("Top 20 Diseases by Degree in Disease–Gene Network", fontsize=16, pad=20)
    # Make labels less cramped
    plt.tight_layout(pad=2.5)
    fig.savefig(OUT_DIR / "fig2_top_diseases_degree.png", dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / "fig2_top_diseases_degree.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)

# Fig 3 – Histogram of Disease Community Sizes
def plot_disease_community_sizes():
    df = pd.read_csv("results/disease_communities.csv")
    sizes = df["disease_count"]

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.hist(sizes, bins=20)
    ax.set_xlabel("Community size (number of diseases)")
    ax.set_ylabel("Frequency")
    ax.set_title("Distribution of Disease Community Sizes (Louvain)", fontsize=16, pad=20)
    plt.tight_layout(pad=2.5)
    fig.savefig(OUT_DIR / "fig3_disease_community_sizes.png", dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / "fig3_disease_community_sizes.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)

# Fig 4 – Top 20 Drugs by Number of Side Effects
def plot_top_drugs_sideeffects():
    df = pd.read_csv("results/top_drugs_sideeffects.csv")
    df = df.sort_values("sideeffect_count")

    # If many names are missing, fall back to chembl_id
    labels = df["drug_name"].fillna(df["chembl_id"])

    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.barh(labels, df["sideeffect_count"])
    ax.set_xlabel("Number of side-effect edges")
    ax.set_ylabel("Drug")
    ax.set_title("Top 20 Drugs by Reported Side Effects", fontsize=16, pad=20)
    plt.tight_layout(pad=2.5)
    fig.savefig(OUT_DIR / "fig4_top_drugs_sideeffects.png", dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / "fig4_top_drugs_sideeffects.pdf", dpi=300, bbox_inches='tight')
    plt.close(fig)

# Main (Runs everything)
if __name__ == "__main__":
    plot_top_genes_degree()
    plot_top_diseases_degree()
    plot_disease_community_sizes()
    plot_top_drugs_sideeffects()
    print("Saved figures to:", OUT_DIR.resolve())

