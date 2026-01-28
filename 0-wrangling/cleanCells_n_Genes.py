#!~/anaconda3/envs/homl3/bin/python

import os
from pathlib import Path
import scanpy as sc
import sys
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import re

def clean_ctmat(adata, gene_thr=0.02, sample_thr=0.02):
    """
    Cleans an AnnData object by filtering out unwanted genes and cells.

    Parameters:
        adata (AnnData): The input AnnData object.
        gene_thr (float): The minimum fraction of cells in which a gene must be expressed to be retained.
        sample_thr (float): The minimum percentile of genes expressed per cell for a cell to be retained.

    Returns:
        AnnData: A new AnnData object with filtered genes and cells.
    """
    # Filter out genes without ENSMUS prefix
    ensmus_genes = [gene for gene in adata.var_names if 'ENSMUS' in gene]
    n_removed = adata.n_vars - len(ensmus_genes)
    if n_removed > 0:
        print(f"  Removing {n_removed} genes without ENSMUS prefix")
        adata = adata[:, ensmus_genes]
    
    # Calculate the number of cells per gene and filter genes
    n_cells_tot = adata.n_obs
    n_cells_thr = gene_thr * n_cells_tot
    gene_counts = (adata.X > 0).sum(axis=0).A1 if isinstance(adata.X, csr_matrix) else (adata.X > 0).sum(axis=0)
    genes_to_keep = gene_counts > n_cells_thr

    print(f"  Gene filtering: keeping {genes_to_keep.sum()}/{len(genes_to_keep)} genes")

    # Filter genes
    adata = adata[:, genes_to_keep]

    # Calculate the number of genes per cell and filter cells
    gene_per_cell = np.array((adata.X > 0).sum(axis=1)).flatten()
    cell_percentiles = pd.Series(gene_per_cell).rank(pct=True)
    cells_to_keep = cell_percentiles > sample_thr


    print(f"  Cell filtering: keeping {cells_to_keep.sum()}/{len(cells_to_keep)} cells")
    # Filter cells
    adata = adata[cells_to_keep, :]

    return adata

def setCellId(adata):
    # Set cell_id as index and remove from columns
    if 'cell_id' in adata.obs.columns:
        adata.obs.index = adata.obs['sample_id'].astype(str) + '-' + adata.obs['cell_id'].astype(str)
        adata.obs.index.name = None  # Remove the index name
        adata.obs = adata.obs.drop(columns=['cell_id'])  # Drop the cell_id column

    return adata


# Paths
input_path = Path("../data/0.raw_data/original") #Jun 10th, 2025 
output_path = Path("../data/0.raw_data/flt")

# Ensure the output directory exists
output_path.mkdir(parents=True, exist_ok=True)

# Loop through all h5ad files in the input path
# for file in input_path.glob("Mic_ROSMAP-Micro*.h5ad"):
for file in input_path.glob("*.h5ad"):
    output_file = output_path / file.name

    # Skip processing if the file already exists in the output directory
    if output_file.exists():
        print(f"Skipping file {file.name}: already exists in the output directory.")
        continue


    cell_type_mapping = {
        'oligodendrocyte': 'Oli',
        'central nervous system macrophage': 'Mic',
        'astrocyte': 'Ast',
        'gluatamatergic neuron': 'Exc',
        'GABAergic interneuron': 'Inh',
        'oligodendrocyte precursor cell': 'Opc'
    }
        
    print(f"Processing file: {file.name}")
    
    # Read the h5ad file
    adata = sc.read_h5ad(file)

    print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes.")

    # Set cell barcodes as obs index + make unique
    adata = setCellId(adata)

    # Filter out unwanted cell types
    adata_filtered = adata[(adata.obs['cell_type'] != 'neural stem cell') & 
                           (adata.obs['cell_type'] != 'brain vascular cell')]
    

    print(f"There were {adata.n_obs} cells. After removing NSC and vascular cells, {adata_filtered.n_obs} remain.")


    cell_types = pd.unique(adata_filtered.obs['cell_type'])
    print(f"There are now {len(cell_types)} cell types: {cell_types}")
    
    for cell_type in cell_types:
        print(f"\nProcessing cell type: {cell_type}")
        
        # Filter for this cell type
        adata_ct = adata_filtered[adata_filtered.obs['cell_type'] == cell_type].copy()
        print(f"  Initial: {adata_ct.n_obs} cells and {adata_ct.n_vars} genes")
        
        # Clean this cell type
        adata_ct_clean = clean_ctmat(adata_ct, gene_thr=0.02, sample_thr=0.02)
        print(f"  After cleaning: {adata_ct_clean.n_obs} cells and {adata_ct_clean.n_vars} genes")
        

        cell_type_abbr = cell_type_mapping.get(cell_type, cell_type.replace(' ', '_').replace('/', '_'))

        # Save the cleaned cell type data
        dataset_name = re.sub(".h5ad", "", Path(file).name)
        output_file = output_path / f"{cell_type_abbr}_{dataset_name}.h5ad"
        adata_ct_clean.write_h5ad(output_file, compression="gzip")
        print(f"  Saved to: {output_file}")

    print()

print("All files have been processed successfully.")


#chmod +x cleanCells_n_Genes.py