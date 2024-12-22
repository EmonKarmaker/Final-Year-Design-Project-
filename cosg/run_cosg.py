import pandas as pd
import scanpy as sc
import cosg
import os

# File paths
expression_file = r"D:\Gene_marker_identification2\GSE104323_10X_expression_data_V2.tab.gz"
metadata_file = r"D:\Gene_marker_identification2\GSE104323_metadata_barcodes_24185cells.txt.gz"

# Check if files exist
if not os.path.exists(expression_file):
    raise FileNotFoundError(f"Expression data file not found: {expression_file}")
if not os.path.exists(metadata_file):
    raise FileNotFoundError(f"Metadata file not found: {metadata_file}")
print("Files found. Proceeding...")

# Load data
expression_data = pd.read_csv(expression_file, sep="\t", compression='gzip', index_col=0)
print("Expression data loaded:", expression_data.shape)

metadata = pd.read_csv(metadata_file, sep="\t", compression='gzip')
print("Metadata loaded:", metadata.shape)

# Check for duplicate cell barcodes in metadata
print("Checking for duplicates in metadata cell barcodes...")
print(metadata["Sample name (24185 single cells)"].duplicated().sum(), "duplicates found in metadata.")

# Drop duplicates in metadata
metadata = metadata.drop_duplicates(subset="Sample name (24185 single cells)")
print("Duplicates removed. Metadata shape:", metadata.shape)

# Align metadata with expression data
metadata.set_index("Sample name (24185 single cells)", inplace=True)
metadata = metadata.reindex(expression_data.columns)

# Handle missing barcodes
missing_barcodes = metadata.index[metadata.isnull().any(axis=1)]
print(f"{len(missing_barcodes)} barcodes missing in metadata. Filling missing rows with NaN.")
metadata.fillna("Unknown", inplace=True)

# Use the correct column name for clusters and convert to categorical
adata = sc.AnnData(expression_data.T)
adata.obs["cluster"] = metadata["characteristics: cell cluster"].astype("category")
print("Cluster column converted to categorical.")
print("AnnData object created.")

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("Data normalized and log-transformed.")

# Check the distribution of samples in each cluster
print("Samples per cluster before filtering:")
print(adata.obs["cluster"].value_counts())

# Filter clusters with at least one sample
valid_clusters = adata.obs["cluster"].value_counts()[adata.obs["cluster"].value_counts() > 0].index
adata = adata[adata.obs["cluster"].isin(valid_clusters)]
print(f"Clusters with sufficient samples: {valid_clusters}")

# Ensure valid_clusters is not empty
if len(valid_clusters) == 0:
    raise ValueError("No valid clusters found in the dataset.")

# Run COSG
cosg.cosg(adata, groupby="cluster", n_genes_user=50)
print("**finished identifying marker genes by COSG**")

# Extract marker genes data
marker_genes = adata.uns["rank_genes_groups"]

# Initialize a list to store rows for the DataFrame
rows = []

# Extract data cluster by cluster
for cluster_idx, cluster_name in enumerate(marker_genes["names"].dtype.names):
    genes = marker_genes["names"][cluster_name]
    scores = marker_genes["scores"][cluster_name]
    logfoldchanges = marker_genes["logfoldchanges"][cluster_name]

    # Combine data for the current cluster
    for gene, score, lfc in zip(genes, scores, logfoldchanges):
        rows.append({"Cluster": cluster_name, "Gene Name": gene, "COSG Score": score, "Log Fold Change": lfc})

# Convert the list of rows into a DataFrame
marker_genes_df = pd.DataFrame(rows)

# Save the results to a CSV
marker_genes_df.to_csv("marker_genes.csv", index=False)
print("Marker genes saved to marker_genes.csv")
