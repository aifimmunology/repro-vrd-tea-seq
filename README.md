# Reproducible analysis of T cell VRd treatment TEA-seq

This repository houses Jupyter Notebooks that are used for reproducible analysis of T cells treated 
with bortezomib, lenalidomide, or dexamethasone and profiled using the trimodal single-cell assay TEA-seq.

These notebooks can be used within the Human Immune System Explorer (HISE) framework for traceable, 
reproducible analysis that results in a Certificate of Reproducibility.

### Authors
Lauren Okada  
Lucas T. Graybuck  

## Contents

### Common files
Files that are used across multiple analyses are stored in `common/`.

### Cell Filtering
Notebooks related to QC and cell type-based cell filtering are stored in `00-cell-filtering/`.

### Cell Type Labeling
Notebooks related to T cell type labeling are stored in `01-cell-type-labeling/`.

### MAST differential expression tests
Notebooks related to the use of the MAST framework to identify differentially expressed genes (DEGs) are stored in `02-mast-deg-testing/`.

### GSEA results
Notebooks related to Gene Set Enrichment Analysis and intersections between MAST results from different conditions are stored in `03-gsea-analysis/`.

### scATAC-seq analysis
Notebooks related to analysis of scATAC-seq data and enrichment of TF motifs near DEGs are stored in `04-atac-analysis/`.

### Differential epitope detection analysis
Notebooks related to differential detection of cell surface antibodies (Antibody-derived tags, ADTs) using linear models are stored in `05-adt-lm-testing/`.

### Flow cytometry data and analysis
Data and Notebooks related to flow cytometry are stored in `06-flow-cytometry/`. This section is partially reproducible, as much of the analysis was carried out by manual gating in FlowJo. Summary statistics and outputs from FlowJo analysis are stored in `06-flow-cytometry/data/`, which are used to generate secondary analyses and figures.

## Figures
Assembly of results into panels for use in figures are stored in the `figures/` directory. Panels and the data backing those panels are in `figures/output/figure_N` or `figures/output/supp_figure_N` where `N` is the figure number.

## Interactive Dash App
An interactive app for exploration of DEG results was implemented in the Dash framework for Python. Code used to generate this app is linked as a submodule to this repository in `repro-vrd-tea-seq-deg-app/`.
