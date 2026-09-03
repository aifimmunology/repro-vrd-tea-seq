
This app enables exploration of genes associated with single Multiple Myeloma drug perturbations in
normal T-cells in-vitro.  

### App usage:  

1) Toggle dropdowns to change volcano plot displayed  
2) Hover over volcano plot to view gene information.  
3) Click 'Annotate Volcano' pencil icon to select gene(s) of interest to label on the volcano plot. Select from dropdown or by clicking/box-selecting from volcano plot, copying and pasting into the box, or selecting a full geneset. Click 'Update Plots' to apply feature selection to all tabs.  
4) Click individual points on the volcano plot to select a single gene to plot on the all-analysis heatmap, or select a gene from the heatmap dropdown.  
5) Customize the heatmap view by clicking 'Heatmap Options' gear icon and changing selections.  

### Analysis Details: 

**DEG Results**  
In vitro DEG's were calculated using MAST (1) for each celltype and timepoint separately. Differences are determined for a given drug relative to the DMSO control at the same timepoint. The MAST model  implemented was:  
`exp ~ treatment + CDR`  
where `CDR` is the cellular detection rate calculated as n_genes detected, mean-centered and scaled to unit variance across all cells.  

Genes exhibiting less than 5% expression in all comparison groups in a cell type were excluded from analysis on a per-celltype basis.  

**Gene Info**  
Gene information is queried in real time using the mygene.info service via the mygene python package (2).   

**DEG Heatmap**  
Heatmap displays the results of DEG analysis with color indicating effect size (log2FC) and significance indicated by '*' (based on selected p-value input).

### References

1. Finak G, McDavid A, Yajima M, Deng J, Gersuk V, Shalek AK, et al. MAST: a flexible statistical framework for assessing transcriptional changes and characterizing heterogeneity in single-cell RNA sequencing data. Genome Biol. 2015;16: 278. [doi:10.1186/s13059-015-0844-5](https://doi.org/10.1186/s13059-015-0844-5)
2. Xin J, Mark A, Afrasiabi C, Tsueng G, Juchler M, Gopal N, et al. High-performance web services for querying gene and variant annotation. Genome Biol. 2016;17: 91. [doi:10.1186/s13059-016-0953-9](https://doi.org/10.1186/s13059-016-0953-9)

### Citation and contributors

**Visualization app developers**  
Lauren Okada  
Garth Kong  
Lucas Graybuck

For additional study details and citation information, visit the study overview page at:  
https://apps.allenimmunology.org/aifi/insights/tcell-vrd/