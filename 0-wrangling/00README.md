

cleanCells_n_Genes.py 
- This script:
    - sets unique cell identifiers in the AnnData object as obs.index
    - filters out genes that are expressed in less than 2% of cells
    - filters out cells that show expression of less than 2% of genes
