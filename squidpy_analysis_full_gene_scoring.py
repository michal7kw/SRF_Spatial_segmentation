# %% [markdown]
# # Squidpy Spatial Gene Set Scoring Analysis
# 
# This script performs gene set scoring across spatial regions using a 40x40 grid.

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
from scipy import stats
from sklearn.preprocessing import StandardScaler

import scanpy as sc
import squidpy as sq
import geopandas as gpd

# %% [markdown]
# ## 2. Gene List Processing Functions

# %%
def parse_gene_list(gene_file_path):
    """
    Parse the gene list from CSV file and extract individual gene names.
    
    Parameters:
    -----------
    gene_file_path : str
        Path to the gene list CSV file
        
    Returns:
    --------
    dict : Dictionary with gene set names as keys and gene lists as values
    """
    gene_sets = {}
    
    with open(gene_file_path, 'r') as f:
        content = f.read().strip()
    
    # Split by commas and process each entry
    entries = content.split(',')
    current_set = None
    current_genes = []
    
    for entry in entries:
        entry = entry.strip()
        if not entry:
            continue
            
        # Check if this looks like a gene set name (contains spaces or is at start)
        if ' ' in entry or entry in ['Epilepsy', 'SETBP1 related']:
            # Save previous set if exists
            if current_set and current_genes:
                gene_sets[current_set] = current_genes
            
            # Start new set
            current_set = entry
            current_genes = []
        else:
            # This is a gene name
            if current_set:
                current_genes.append(entry)
    
    # Save the last set
    if current_set and current_genes:
        gene_sets[current_set] = current_genes
    
    # If no clear sets found, treat all as one set
    if not gene_sets:
        all_genes = [entry.strip() for entry in entries if entry.strip()]
        gene_sets['All_Genes'] = all_genes
    
