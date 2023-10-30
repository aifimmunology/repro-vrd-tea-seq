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
Notebooks for assembly of results into panels for use in figures are stored in the `figures/` directory. Panels and the data backing those panels are in `figures/output/figure_N` or `figures/output/supp_figure_N` where `N` is the figure number.

## Tables
Notebooks for assembly of results into supplementary tables are stored in the `tables/` directory. Raw versions of supplementary tables are compressed using `gzip` and stored in `tables/output/`.

## Interactive Dash App
An interactive app for exploration of DEG results was implemented in the Dash framework for Python. Code used to generate this app is linked as a submodule to this repository in `repro-vrd-tea-seq-deg-app/`.

## Source Data

The notebooks in this repository are designed to work within the HISE system, where raw data are stored in a central repository and are accessed using the HISE SDK. However, the input data are also available for use outside of HISE on dbGaP and GEO.

### Raw data
Raw, FASTQ-level data from our experiments are deposited in the database of Genotypes and Phenotypes (dbGaP) at accession number [phs003430.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003430.v1.p1).

### Processed data
Processed data, including scRNA-seq gene counts, ADT feature counts, and scATAC-seq fragment counts, are available in the Gene Expression Omnibus at accession number [GSE236422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236422).

## Legal Information

### License

The license for this repository is available on Github in the file [LICENSE in this repository](https://github.com/aifimmunology/repro-vrd-tea-seq/blob/master/LICENSE)

### Level of Support

We are not currently supporting this code, but simply releasing it to the community AS IS but are not able to provide any guarantees of support. The community is welcome to submit issues, but you should not expect an active response.
