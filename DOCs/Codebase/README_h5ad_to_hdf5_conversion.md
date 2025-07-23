# H5AD to HDF5 Conversion Tools

This directory contains tools to convert AnnData `.h5ad` files to standard HDF5 `.hdf5` format.

## Files

- `h5ad_to_hdf5_converter.py` - Main conversion script
- `verify_hdf5.py` - Verification and inspection script
- `README_h5ad_to_hdf5_conversion.md` - This documentation

## Conversion Results

✅ **Successfully converted**: `region_R3/202503071102_SESSA-p30-E165_VMSC10702_region_R3.h5ad` → `region_R3/202503071102_SESSA-p30-E165_VMSC10702_region_R3.hdf5`

- **Original file size**: 182 MB (.h5ad)
- **Converted file size**: 49.19 MB (.hdf5)
- **Data preserved**: 43,284 cells × 960 genes

## Usage

### Basic Conversion

```bash
# Convert the specific R3 file (automatic detection)
python h5ad_to_hdf5_converter.py

# Convert any h5ad file
python h5ad_to_hdf5_converter.py input_file.h5ad

# Specify output file
python h5ad_to_hdf5_converter.py input_file.h5ad -o output_file.hdf5

# Use different compression
python h5ad_to_hdf5_converter.py input_file.h5ad -c lzf
```

### Verification

```bash
# Inspect the converted HDF5 file structure
python verify_hdf5.py
```

## Data Structure in HDF5

The converted HDF5 file contains the following structure:

```
/
├── X                           # Expression matrix (43284 × 960)
├── obs/                        # Cell/observation metadata
│   ├── _index                  # Cell IDs
│   ├── center_x, center_y      # Spatial coordinates
│   ├── leiden_codes            # Cluster assignments
│   ├── leiden_categories       # Cluster labels
│   ├── n_counts, n_genes       # QC metrics
│   └── volume                  # Cell volumes
├── var/                        # Gene/variable metadata
│   ├── _index                  # Gene names
│   ├── mean                    # Gene expression means
│   └── std                     # Gene expression standard deviations
├── obsm/                       # Multi-dimensional observations
│   ├── X_pca                   # PCA coordinates (43284 × 50)
│   └── X_umap                  # UMAP coordinates (43284 × 2)
├── varm/                       # Multi-dimensional variables
│   └── PCs                     # Principal components (960 × 50)
├── obsp/                       # Pairwise observations (sparse matrices)
│   ├── connectivities/         # Cell-cell connectivity
│   └── distances/              # Cell-cell distances
└── uns/                        # Unstructured annotations
    ├── leiden/                 # Leiden clustering parameters
    ├── pca/                    # PCA parameters
    ├── rank_genes_groups/      # Differential expression results
    └── umap/                   # UMAP parameters
```

## Loading Data from HDF5

### Python with h5py

```python
import h5py
import numpy as np

# Open the file
with h5py.File('region_R3/202503071102_SESSA-p30-E165_VMSC10702_region_R3.hdf5', 'r') as f:
    # Load expression matrix
    expression_matrix = f['X'][:]
    
    # Load cell metadata
    cell_ids = [idx.decode('utf-8') for idx in f['obs']['_index'][:]]
    leiden_clusters = f['obs']['leiden_codes'][:]
    
    # Load spatial coordinates
    center_x = f['obs']['center_x'][:]
    center_y = f['obs']['center_y'][:]
    
    # Load UMAP coordinates
    umap_coords = f['obsm']['X_umap'][:]
    
    # Load gene names
    gene_names = [gene.decode('utf-8') for gene in f['var']['_index'][:]]
```

### Python with pandas

```python
import h5py
import pandas as pd

with h5py.File('region_R3/202503071102_SESSA-p30-E165_VMSC10702_region_R3.hdf5', 'r') as f:
    # Create cell metadata DataFrame
    cell_metadata = pd.DataFrame({
        'cell_id': [idx.decode('utf-8') for idx in f['obs']['_index'][:]],
        'center_x': f['obs']['center_x'][:],
        'center_y': f['obs']['center_y'][:],
        'leiden_cluster': f['obs']['leiden_codes'][:],
        'n_counts': f['obs']['n_counts'][:],
        'n_genes': f['obs']['n_genes'][:],
        'volume': f['obs']['volume'][:]
    })
    
    # Create gene metadata DataFrame
    gene_metadata = pd.DataFrame({
        'gene_name': [gene.decode('utf-8') for gene in f['var']['_index'][:]],
        'mean_expression': f['var']['mean'][:],
        'std_expression': f['var']['std'][:]
    })
```

## Features Preserved

✅ **Expression Matrix**: Dense matrix (43,284 cells × 960 genes)  
✅ **Cell Metadata**: Spatial coordinates, clustering, QC metrics  
✅ **Gene Metadata**: Expression statistics  
✅ **Dimensionality Reduction**: PCA (50 components), UMAP (2D)  
✅ **Clustering**: Leiden clusters (22 clusters)  
✅ **Spatial Information**: Cell coordinates and volumes  
✅ **Connectivity**: Cell-cell relationships (sparse matrices)  
✅ **Differential Expression**: Rank genes groups results  

## Compression

The HDF5 format provides excellent compression:
- **Original .h5ad**: 182 MB
- **Converted .hdf5**: 49.19 MB
- **Compression ratio**: ~73% size reduction

## Compatibility

The generated HDF5 files are compatible with:
- Python (h5py, pandas, numpy)
- R (rhdf5, hdf5r)
- MATLAB (built-in HDF5 support)
- Julia (HDF5.jl)
- Any tool that supports standard HDF5 format

## Notes

- Categorical data is stored as codes + categories for efficiency
- Sparse matrices are stored in CSR format with separate data/indices/indptr arrays
- String data is UTF-8 encoded
- Some complex nested structures in `uns` may be skipped if not HDF5-compatible
- File metadata includes conversion information and original file reference 