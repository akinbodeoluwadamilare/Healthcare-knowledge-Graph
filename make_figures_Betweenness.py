import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("results/top_genes_betweenness.csv")

# Sort so the lowest is at the bottom of the barh plot
df = df.sort_values("betweenness")

plt.figure(figsize=(10, 6))
plt.barh(df["symbol"], df["betweenness"])

plt.xlabel("Betweenness centrality")
plt.ylabel("Gene")
plt.title("Top 20 Genes by Betweenness in Disease–Gene Network")

plt.tight_layout()
plt.savefig("figures/fig5_top_genes_betweenness.png", dpi=300)
plt.close()
