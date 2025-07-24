import os
import scanpy as sc
import tangram as tg
import pandas as pd
import numpy as np

# Set working directory
os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation")

# Define file paths
spatial_data_path = "DATA/p30-E165/R4/data.h5ad"
sc_data_path = "DATA/RNA_counts/Emx1_Ctrl/Emx1_Ctrl_processed.h5ad"
output_dir = "DATA_INTEGRATION"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "p30_R4_tangram_integrated.h5ad")

# 1. Load Data
print("Loading spatial data...")
adata_sp = sc.read_h5ad(spatial_data_path)
print("Spatial data loaded:")
print(adata_sp)

print("\nLoading single-cell data...")
adata_sc = sc.read_h5ad(sc_data_path)
print("Single-cell data loaded:")
print(adata_sc)

# Ensure spatial coordinates are in .obsm
if 'spatial' not in adata_sp.obsm_keys():
    if 'center_x' in adata_sp.obs.columns and 'center_y' in adata_sp.obs.columns:
        adata_sp.obsm['spatial'] = adata_sp.obs[['center_x', 'center_y']].values
        print("Added spatial coordinates to adata_sp.obsm['spatial']")
    else:
        raise KeyError("Spatial coordinates ('center_x', 'center_y') not found in adata_sp.obs")

# 2. Pre-process Data for Tangram
print("\nPreprocessing data for Tangram...")

# Find marker genes for training. Using highly variable genes as a proxy.
# The scRNA-seq data is already processed, so we can directly find highly variable genes.
print("Finding highly variable genes in single-cell data for training...")
sc.pp.highly_variable_genes(adata_sc, n_top_genes=1000, flavor='seurat')
training_genes = adata_sc.var_names[adata_sc.var.highly_variable].tolist()
print(f"Found {len(training_genes)} highly variable genes.")

# Pre-process AnnDatas for Tangram
tg.pp_adatas(adata_sc, adata_sp, genes=training_genes)

# 3. Map single cells to space
print("\nMapping single cells to space...")
# Check for GPU availability, otherwise use CPU
try:
    import torch
    if torch.cuda.is_available():
        device = "cuda:0"
        print("Using GPU for mapping.")
    else:
        device = "cpu"
        print("Using CPU for mapping.")
except ImportError:
    device = "cpu"
    print("PyTorch not found. Using CPU for mapping.")

ad_map = tg.map_cells_to_space(
    adata_sc=adata_sc,
    adata_sp=adata_sp,
    device=device,
    mode='cells',
    density_prior='uniform', # MERSCOPE is single-cell resolution
    num_epochs=1000,
)

print("Mapping complete.")
print("Resulting mapping AnnData:")
print(ad_map)

# 4. Project gene expression
print("\nProjecting gene expression onto spatial data...")
adata_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)
print("Gene projection complete.")
print("Projected gene expression AnnData:")
print(adata_ge)

# 5. Save the integrated data
print(f"\nSaving integrated AnnData to {output_path}")
adata_ge.write(output_path)

print("\nTangram integration script finished successfully!")
print(f"The final integrated data is saved at: {output_path}")