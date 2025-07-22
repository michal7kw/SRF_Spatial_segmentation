# %% [markdown]
# # Squidpy analysis of Vizgen data with Spatial Gene Set Scoring


# %% [markdown]
# ## 1. Imports

# %%
import numpy as np
import pandas as pd
from pathlib import Path
import shutil
import os
from matplotlib import pyplot as plt
import seaborn as sns
import re
from typing import Dict, List

import scanpy as sc
import squidpy as sq
import geopandas as gpd

# %% [markdown]
# ## 2. Setup

# %%
# Create a directory to save the results
os.chdir("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation")

# %% [markdown]
# ## 3. Load Data

# %%
data_dir = "./DATA/p0-p7/R1"
results_dir = os.path.join(data_dir, "analysis_results")
os.makedirs(results_dir, exist_ok=True)

# %%
adata = sq.read.vizgen(
    path=data_dir,
    counts_file="cell_by_gene.csv",
    meta_file="cell_metadata.csv",
)
print("Data loaded:")
print(adata)


# %% [markdown]
# ### Get Library ID

# %%
# Get the library id from the anndata object
try:
    library_id = list(adata.uns['spatial'].keys())[0]
    print(f"Using library_id: {library_id}")
except (KeyError, IndexError):
    print("Could not automatically determine library_id. Spatial plots may fail.")
    library_id = None


# %% [markdown]
# ### Load Cell Boundaries

# %%
# Load the cell boundaries and add them to the anndata object
try:
    segmentation_path = os.path.join(data_dir, "cell_boundaries.parquet")
    boundaries = gpd.read_parquet(segmentation_path)
    adata.uns['spatial'][library_id]['segmentations'] = boundaries
    print("Successfully loaded cell boundaries.")
except Exception as e:
    print(f"Could not load cell boundaries: {e}")


# %% [markdown]
# ## 4. Pre-processing and QC
#
# We perform standard pre-processing and quality control steps.

# %%
adata.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
axs[0].set_title("Total transcripts per cell")
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, ax=axs[1])
axs[1].set_title("Unique transcripts per cell")
sns.histplot(adata.obs.groupby("fov").sum()["total_counts"], kde=False, ax=axs[2])
axs[2].set_title("Transcripts per FOV")
sns.histplot(adata.obs["volume"], kde=False, ax=axs[3])
axs[3].set_title("Volume of segmented cells")
fig.tight_layout()
plt.savefig(os.path.join(results_dir, "qc_metrics_distribution.png"))
print("Saved QC metrics distribution plot.")


# %% [markdown]
# ### Filtering

# %%
# Filter cells with low expression and genes that are expressed in too few cells.
print(f"Number of cells before filtering: {adata.n_obs}")
sc.pp.filter_cells(adata, min_counts=50)
print(f"Number of cells after filtering by counts: {adata.n_obs}")

print(f"Number of genes before filtering: {adata.n_vars}")
sc.pp.filter_genes(adata, min_cells=10)
print(f"Number of genes after filtering by cells: {adata.n_vars}")


# %% [markdown]
# ### Normalization and Scaling

# %%
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
print("Normalization and scaling complete.")

# %% [markdown]
# ## 5. Dimensionality Reduction and Clustering

# %%
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0)
print("Dimensionality reduction and clustering complete.")

# %% [markdown]
# ## 6. Visualization

# %% [markdown]
# ### UMAP

# %%
sc.pl.umap(adata, color=["leiden"], size=10, show=False, save="_leiden.png")
# Move the file to the results directory
Path("figures/umap_leiden.png").rename(f"{results_dir}/umap_leiden.png")
shutil.rmtree("figures")
print("Saved UMAP plot.")


# %% [markdown]
# ### Spatial Scatter

# %%
sq.pl.spatial_scatter(
    adata,
    library_id=library_id,
    color="leiden",
    img=False,
    figsize=(15, 15),
    size=10,
    save="spatial_leiden.png"
)
# Move the file to the results directory
if os.path.exists("figures/spatial_leiden.png"):
    Path("figures/spatial_leiden.png").rename(f"{results_dir}/spatial_leiden.png")
    shutil.rmtree("figures")
    print("Saved spatial scatter plot.")
else:
    print("Could not save spatial scatter plot.")

# %% [markdown]
# ## 7. Spatial Analysis
# Spatial organization of the clustered cells.

# %% [markdown]
# ### Neighborhood Enrichment
# Clusters that are spatially co-enriched.

# %%
sq.gr.spatial_neighbors(adata, coord_type="generic", spatial_key="spatial")
sq.gr.nhood_enrichment(adata, cluster_key="leiden")
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    method="average",
    cmap="inferno",
    vmin=-50,
    vmax=100,
    figsize=(7, 7),
    save="_enrichment.png",
)
Path("figures/_enrichment.png").rename(f"{results_dir}/neighborhood_enrichment.png")
shutil.rmtree("figures")
print("Saved neighborhood enrichment plot.")

# %% [markdown]
# ### Spatial Autocorrelation (Moran's I)
# Identify genes that show a non-random spatial distribution.

# %%
sq.gr.spatial_autocorr(adata, mode="moran")
print("Top 10 spatially autocorrelated genes:")
print(adata.uns["moranI"].head(10))

# %% [markdown]
# ### Visualize top spatially autocorrelated genes

# %%
top_autocorr = adata.uns["moranI"]["I"].sort_values(ascending=False).head(4).index.tolist()
sq.pl.spatial_scatter(
    adata,
    library_id=library_id,
    color=top_autocorr,
    size=5,
    cmap="Reds",
    img=False,
    figsize=(10, 10),
    save="_top_autocorr.png"
)
if os.path.exists("figures/spatial_top_autocorr.png"):
    Path("figures/spatial_top_autocorr.png").rename(f"{results_dir}/spatial_top_autocorr.png")
    os.rmdir("figures")
    print("Saved top spatially autocorrelated genes plot.")
else:
    print("Could not save top spatially autocorrelated genes plot.")

# %% [markdown]
# ## 8. Gene Set Scoring on Spatial Grid

# %% [markdown]
# ### Load and Parse Gene Lists

# %%
def parse_gene_lists(file_path: str) -> Dict[str, List[str]]:
    """
    Parse gene lists from CSV file.
    Expected format: gene_set_name,gene1,gene2,gene3,...
    """
    gene_sets = {}
    
    with open(file_path, 'r') as f:
        content = f.read().strip()
    
    # Split by gene set names (assuming they end with comma and start new line or continue)
    # The format appears to be: SetName,gene1,gene2,...SetName2,gene1,gene2,...
    
    # Find all gene set names (capitalized words followed by comma)
    gene_set_pattern = r'([A-Z][a-zA-Z\s_]+?)(?=,)'
    matches = list(re.finditer(gene_set_pattern, content))
    
    if not matches:
        # Fallback: treat as single gene set
        genes = [gene.strip() for gene in content.split(',') if gene.strip()]
        gene_sets['All_genes'] = genes
        return gene_sets
    
    for i, match in enumerate(matches):
        set_name = match.group(1).strip()
        start_pos = match.end() + 1  # Skip the comma
        
        # Find end position (next gene set or end of string)
        if i + 1 < len(matches):
            end_pos = matches[i + 1].start()
        else:
            end_pos = len(content)
        
        # Extract genes for this set
        gene_string = content[start_pos:end_pos]
        genes = [gene.strip() for gene in gene_string.split(',') if gene.strip()]
        
        # Clean gene set name
        clean_name = re.sub(r'[^a-zA-Z0-9_]', '_', set_name)
        gene_sets[clean_name] = genes
    
    return gene_sets

# Load gene lists
gene_list_path = "./GENE_LISTS/input/bp_others.csv"
gene_sets = parse_gene_lists(gene_list_path)

print("Loaded gene sets:")
for set_name, genes in gene_sets.items():
    print(f"  {set_name}: {len(genes)} genes")
    print(f"    First 5 genes: {genes[:5]}")

# %% [markdown]
# ### Create Spatial Grid

# %%
def create_spatial_grid(adata, grid_size=(40, 40)):
    """
    Create a spatial grid and assign each cell to a grid position.
    """
    # Get spatial coordinates
    spatial_coords = adata.obsm['spatial']
    
    # Calculate grid boundaries
    x_min, x_max = spatial_coords[:, 0].min(), spatial_coords[:, 0].max()
    y_min, y_max = spatial_coords[:, 1].min(), spatial_coords[:, 1].max()
    
    # Create grid
    x_bins = np.linspace(x_min, x_max, grid_size[0] + 1)
    y_bins = np.linspace(y_min, y_max, grid_size[1] + 1)
    
    # Assign cells to grid positions
    x_indices = np.digitize(spatial_coords[:, 0], x_bins) - 1
    y_indices = np.digitize(spatial_coords[:, 1], y_bins) - 1
    
    # Ensure indices are within bounds
    x_indices = np.clip(x_indices, 0, grid_size[0] - 1)
    y_indices = np.clip(y_indices, 0, grid_size[1] - 1)
    
    # Create grid identifiers
    grid_ids = x_indices * grid_size[1] + y_indices
    
    # Add to adata
    adata.obs['grid_x'] = x_indices
    adata.obs['grid_y'] = y_indices
    adata.obs['grid_id'] = grid_ids
    
    return x_bins, y_bins

# Create spatial grid
print("Creating 40x40 spatial grid...")
x_bins, y_bins = create_spatial_grid(adata, grid_size=(40, 40))
print(f"Grid created. Cells distributed across {adata.obs['grid_id'].nunique()} grid positions.")

# %% [markdown]
# ### Perform Gene Set Scoring per Grid Cell

# %%
def score_gene_sets_per_grid(adata, gene_sets, grid_size=(40, 40)):
    """
    Score gene sets for each grid cell.
    """
    # Create a copy of adata for scoring
    adata_scoring = adata.copy()
    
    # Results storage
    grid_scores = {}
    
    for set_name, genes in gene_sets.items():
        print(f"Scoring gene set: {set_name}")
        
        # Filter genes that exist in the dataset
        available_genes = [g for g in genes if g in adata_scoring.var_names]
        missing_genes = [g for g in genes if g not in adata_scoring.var_names]
        
        print(f"  Available genes: {len(available_genes)}/{len(genes)}")
        if missing_genes:
            print(f"  Missing genes: {missing_genes[:10]}{'...' if len(missing_genes) > 10 else ''}")
        
        if len(available_genes) == 0:
            print(f"  Skipping {set_name} - no genes found in dataset")
            continue
        
        # Score gene set using scanpy
        sc.tl.score_genes(adata_scoring, available_genes, score_name=f'{set_name}_score')
        
        # Calculate mean score per grid cell
        grid_scores_for_set = []
        grid_positions = []
        
        for grid_id in range(grid_size[0] * grid_size[1]):
            cells_in_grid = adata_scoring.obs['grid_id'] == grid_id
            
            if cells_in_grid.sum() > 0:
                mean_score = adata_scoring.obs.loc[cells_in_grid, f'{set_name}_score'].mean()
                grid_scores_for_set.append(mean_score)
            else:
                grid_scores_for_set.append(np.nan)
            
            # Calculate grid position
            x_pos = grid_id // grid_size[1]
            y_pos = grid_id % grid_size[1]
            grid_positions.append((x_pos, y_pos))
        
        grid_scores[set_name] = {
            'scores': np.array(grid_scores_for_set),
            'positions': grid_positions
        }
    
    return grid_scores

# Perform gene set scoring
print("Performing gene set scoring per grid cell...")
grid_scores = score_gene_sets_per_grid(adata, gene_sets, grid_size=(40, 40))

# %% [markdown]
# ### Visualize Gene Set Scores on Spatial Grid

# %%
def plot_grid_scores(grid_scores, grid_size=(40, 40), save_dir=None):
    """
    Plot gene set scores as heatmaps on the spatial grid.
    """
    n_sets = len(grid_scores)
    if n_sets == 0:
        print("No gene sets to plot")
        return
    
    # Calculate subplot layout
    n_cols = min(3, n_sets)
    n_rows = (n_sets + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_sets == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for idx, (set_name, data) in enumerate(grid_scores.items()):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]
        
        # Reshape scores to grid
        score_grid = data['scores'].reshape(grid_size)
        
        # Create heatmap
        im = ax.imshow(score_grid, cmap='viridis', aspect='equal', origin='lower')
        ax.set_title(f'{set_name}\nGene Set Score')
        ax.set_xlabel('Grid X')
        ax.set_ylabel('Grid Y')
        
        # Add colorbar
        plt.colorbar(im, ax=ax, shrink=0.8)
    
    # Hide empty subplots
    for idx in range(n_sets, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]
        ax.set_visible(False)
    
    plt.tight_layout()
    
    if save_dir:
        plt.savefig(os.path.join(save_dir, "spatial_gene_set_scores.png"), dpi=300, bbox_inches='tight')
        print(f"Saved spatial gene set scores plot to {save_dir}")
    
    plt.show()

# Plot results
plot_grid_scores(grid_scores, grid_size=(40, 40), save_dir=results_dir)

# %% [markdown]
# ### Save Results

# %%
# Save grid scores to CSV
for set_name, data in grid_scores.items():
    # Create DataFrame with grid positions and scores
    df = pd.DataFrame({
        'grid_x': [pos[0] for pos in data['positions']],
        'grid_y': [pos[1] for pos in data['positions']],
        'score': data['scores']
    })
    
    # Save to CSV
    output_file = os.path.join(results_dir, f"grid_scores_{set_name}.csv")
    df.to_csv(output_file, index=False)
    print(f"Saved {set_name} grid scores to {output_file}")

# Save summary statistics
summary_stats = []
for set_name, data in grid_scores.items():
    scores = data['scores']
    valid_scores = scores[~np.isnan(scores)]
    
    if len(valid_scores) > 0:
        stats = {
            'gene_set': set_name,
            'n_grid_cells': len(valid_scores),
            'mean_score': np.mean(valid_scores),
            'std_score': np.std(valid_scores),
            'min_score': np.min(valid_scores),
            'max_score': np.max(valid_scores)
        }
        summary_stats.append(stats)

summary_df = pd.DataFrame(summary_stats)
summary_file = os.path.join(results_dir, "gene_set_scoring_summary.csv")
summary_df.to_csv(summary_file, index=False)
print(f"Saved summary statistics to {summary_file}")

print("\nSpatial gene set scoring analysis complete!")
print("Results saved in the 'analysis_results' directory:")
print("- spatial_gene_set_scores.png: Heatmap visualization")
print("- grid_scores_*.csv: Individual gene set scores per grid cell")
print("- gene_set_scoring_summary.csv: Summary statistics")

# %%
