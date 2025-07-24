import scanpy as sc
import pandas as pd
import os

# Define file paths
h5ad_path = "DATA_INTEGRATION/p30_R4_tangram_integrated.h5ad"
output_csv_path = "DATA_INTEGRATION/p30_R4_tangram_integrated.csv"

# Load the AnnData object
print(f"Loading AnnData from {h5ad_path}...")
adata = sc.read_h5ad(h5ad_path)
print("AnnData loaded.")

# Create a pandas DataFrame from the expression data
print("Creating DataFrame from AnnData.X...")
df = pd.DataFrame(adata.X, index=adata.obs.index, columns=adata.var.index)

# Save the DataFrame to a CSV file
print(f"Saving DataFrame to {output_csv_path}...")
df.to_csv(output_csv_path)

print(f"Conversion to CSV complete. File saved at: {output_csv_path}")