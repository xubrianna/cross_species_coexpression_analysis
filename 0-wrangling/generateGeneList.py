import os
import sys
import json
from pathlib import Path
import anndata as ad
import scanpy as sc

def analyze_gene_dict(gene_dict):
    """
    Processes a gene dictionary to generate:
    1. A list of genes found per cell type across all studies (stringent).
    2. A list of genes found per cell type in at least 70% of the studies (less stringent).
    2. A list of genes found per cell type in at least 1/3 of the studies (less stringent).
    3. A dictionary of pairwise intersections of "less stringent" gene lists between cell types.
    
    Args:
        gene_dict (dict): The input gene dictionary {cell_type -> {study_name -> [genes]}}.
    
    Returns:
        dict: Results with keys:
            - 'stringent': Genes found across all studies per cell type.
            - 'less_stringent': Genes found in at least 70% of studies per cell type.
            - 'pairwise_intersections': Intersections of "less stringent" lists between cell types.
    """
    from itertools import combinations

    # Initialize result dictionary
    results = {
        "stringent": {},
        "less_stringent": {},
        "pairwise_intersections": {}
    }

    # Process each cell type in the gene dictionary
    for cell_type, studies in gene_dict.items():
        study_genes = list(studies.values())  # List of gene lists for each study
        num_studies = len(study_genes)

        # Compute stringent list: Genes found in all studies
        stringent_genes = set(study_genes[0])
        for genes in study_genes[1:]:
            stringent_genes.intersection_update(genes)
        results["stringent"][cell_type] = list(stringent_genes)

        # Compute less stringent list: Genes found in at least 30% of studies (1/3)
        gene_counts = {}
        for genes in study_genes:
            for gene in genes:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        less_stringent_genes = [
            gene for gene, count in gene_counts.items()
            if count >= 0.3 * num_studies
        ]
        results["less_stringent"][cell_type] = less_stringent_genes

    # Compute pairwise intersections of less stringent gene lists between cell types
    cell_types = list(results["less_stringent"].keys())
    for cell_type1, cell_type2 in combinations(cell_types, 2):
        intersect_genes = set(results["less_stringent"][cell_type1]).intersection(
            results["less_stringent"][cell_type2]
        )
        results["pairwise_intersections"][(cell_type1, cell_type2)] = list(intersect_genes)

    return results


def extract_metadata_from_filename(filename):
    """
    Extract cell type and study name from a filename.
    Assumes the format: "CellType_StudyName_Year.h5ad".
    """
    name_parts = filename.split('_')
    cell_type = name_parts[0]  # First part before the first underscore
    study_name = name_parts[-1].split('.')[0] 
    return cell_type, study_name

def process_h5ad_files(directory):
    """
    Process all .h5ad files in the specified directory, extract metadata,
    and populate a dictionary with cell type and study information.
    """
    directory = Path(directory)  # Ensure the directory is a Path object
    gene_dict = {}  # Dictionary to store metadata information

    # Iterate through all .h5ad files in the directory
    for file_path in directory.glob("*.h5ad"):
        print(f"\nProcessing file: {file_path.name}")

        # Load the AnnData object
        adata = sc.read_h5ad(file_path)
        print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes.")

        # Extract metadata (cell type, study name)
        cell_type, study_name = extract_metadata_from_filename(file_path.name)

        # Populate the dictionary
        if cell_type not in gene_dict:
            gene_dict[cell_type] = {}  # Add cell type if not already present
        if study_name not in gene_dict[cell_type]:
            gene_dict[cell_type][study_name] = []  # Add study name if not already present

        gene_dict[cell_type][study_name] = adata.var.index.tolist()  # Store all gene names

        print(f"Stored {len(adata.var.index)} genes for {cell_type} - {study_name}.")

    return gene_dict


input_dir = Path("../data/0.raw_data/flt")
output_dir = Path('../data/1.slct-genes') 
#output_file = '/home/nelazzabi/rewiring/data/genes/genes-listsV2.json'
output_file = '../data/1.slct-genes/genes-lists.json' #June 10th, 2025 added ROSMAP-MICRO

# Utility function for JSON serialization
def convert_tuple_keys(obj):
    if isinstance(obj, dict):
        return {str(k): convert_tuple_keys(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_tuple_keys(item) for item in obj]
    else:
        return obj


# Step 1: Extract genes per study
print("\nExtracting genes...")
gene_dict = process_h5ad_files(input_dir)

# Step 2: Extract intersect sets and pairs
print("\nAnalyzing gene intersections...")
results = analyze_gene_dict(gene_dict)
results_serializable = convert_tuple_keys(results)

# Step 3: Save results
print("\nSaving results to JSON...")
with open(output_file, 'w') as f:
    json.dump(results_serializable, f, indent=4)

# Step 4: Print summary of results
print("\nSummary of stringent genes per cell type:")
for cell_type, genes in results["stringent"].items():
    print(f"{cell_type}: {len(genes)} genes")

print("\nSummary of less stringent genes per cell type:")
for cell_type, genes in results["less_stringent"].items():
    print(f"{cell_type}: {len(genes)} genes")

print("\nPairwise intersections of less stringent genes:")
for pair, intersect_genes in results["pairwise_intersections"].items():
    print(f"{pair}: {len(intersect_genes)} genes")

print("\nProcessing completed successfully.")
