import pandas as pd

# Load the marker genes file
file_path = "marker_genes.csv"
data = pd.read_csv(file_path)

# Display the first few rows of the data
print("First few rows of the data:")
print(data.head())

# Check the structure of the data
print("\nData structure:")
print(data.info())

# Group the data by clusters and count the number of genes per cluster
print("\nNumber of marker genes per cluster:")
gene_counts = data["Cluster"].value_counts()
print(gene_counts)

# Find the top marker genes for each cluster based on COSG Scores
print("\nTop marker genes per cluster:")
top_genes = data.sort_values(by="COSG Score", ascending=False).groupby("Cluster").head(3)  # Top 3 genes per cluster
print(top_genes)

# Save insights into a summary file
gene_counts.to_csv("gene_counts_summary.csv", index=True, header=["Number of Marker Genes"])
top_genes.to_csv("top_marker_genes.csv", index=False)

print("\nAnalysis complete. Summaries saved as 'gene_counts_summary.csv' and 'top_marker_genes.csv'.")
