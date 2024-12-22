import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the CSV files
gene_counts = pd.read_csv("gene_counts_summary.csv", index_col=0)
top_genes = pd.read_csv("top_marker_genes.csv")

# Bar Chart: Number of Marker Genes per Cluster
plt.figure(figsize=(10, 6))
gene_counts.plot(kind="bar", legend=False, color="skyblue", edgecolor="black")
plt.title("Number of Marker Genes per Cluster", fontsize=16)
plt.ylabel("Number of Marker Genes", fontsize=12)
plt.xlabel("Cluster", fontsize=12)
plt.xticks(rotation=45, fontsize=10)
plt.tight_layout()
plt.savefig("marker_genes_bar_chart.png")
plt.show()

# Table: Highlight Top Marker Genes for a Few Clusters
# Select a few clusters for demonstration
selected_clusters = ["Astro-adult", "Immature-Pyr", "Endothelial"]
highlight_genes = top_genes[top_genes["Cluster"].isin(selected_clusters)]
print("Highlighted Top Marker Genes for Selected Clusters:")
print(highlight_genes)

# Scatter Plot: COSG Scores vs Log Fold Changes for a Specific Cluster
# Choose a specific cluster to plot
cluster_to_plot = "Astro-adult"
astro_genes = top_genes[top_genes["Cluster"] == cluster_to_plot]

plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=astro_genes,
    x="Log Fold Change",
    y="COSG Score",
    hue="Gene Name",
    palette="viridis",
    s=100,
)
plt.title(f"COSG Scores vs Log Fold Changes: {cluster_to_plot}", fontsize=16)
plt.xlabel("Log Fold Change", fontsize=12)
plt.ylabel("COSG Score", fontsize=12)
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=10)
plt.tight_layout()
plt.savefig(f"{cluster_to_plot}_scatter_plot.png")
plt.show()
