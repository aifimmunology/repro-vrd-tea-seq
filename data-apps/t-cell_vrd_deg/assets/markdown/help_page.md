# About the Project

This app lets you explore single-cell RNA sequencing data from a longitudinal study of multiple myeloma patients profiled across diagnosis, chemotherapy, stem cell transplant, and recovery. Over 6.5 million cells from matched bone marrow and blood samples were profiled across 64 and 71 cell subtypes respectively, capturing how the tumor and its treatment reshape immune cell populations over time.

---

## Data Analysis

Bone marrow (BMMC) and peripheral blood (PBMC) samples were collected from 17 newly diagnosed multiple myeloma patients at up to five timepoints: pre-treatment, end of induction (VRd), 90 days post-transplant, and 1 and 2 years post-transplant, alongside 10 healthy bone marrow donors and 32 healthy blood donors. BMMCs were profiled using 5' CITE-seq and PBMCs using 3' scRNA-seq (10x Genomics), yielding ~1.1M and ~5.4M cells annotated into 64 and 71 subtypes respectively. Malignant plasma cells were isolated computationally and analyzed separately using copy number variation inference.

---

## Additional Resources 

Raw and processed scRNA-seq data are deposited in GEO under accessions GSE309595 (PBMC) and GSE309593 (BMMC). Analysis code is available on GitHub and will be archived on Zenodo upon publication.

Processed flow cytometry, Olink proteomics, and clinical metadata are available through the Human Immune System Explorer (HISE) portal. A companion interactive viewer for the full multimodal dataset is available [here](https://apps.allenimmunology.org/aifi/insights/ndmm/).

---

## References

1. Rajkumar, S.V. (2022). Multiple myeloma: 2022 update on diagnosis, risk stratification, and management. Am. J. Hematol. 97, 1086–1107.

2. Blimark, C.H. et al. (2025). Risk of infections in multiple myeloma: a population-based study on 8,672 patients from the Swedish Myeloma Registry. Haematologica 110, 163–172.

3. Lutz, R. et al. (2024). Multiple myeloma long-term survivors exhibit sustained immune alterations decades after first-line therapy. Nat. Commun. 15, 10396.

4. Gong, Q. et al. (2025). Multi-omic profiling reveals age-related immune dynamics in healthy adults. Nature 648, 696–706.

5. Alameh, M.-G. et al. (2021). Lipid nanoparticles enhance the efficacy of mRNA and protein subunit vaccines by inducing robust T follicular helper cell and humoral responses. Immunity 54, 2877–2892.

6. Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 550.

7. Korsunsky, I. et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nat. Methods 16, 1289–1296.

---

# Contributors

Please see the original manuscript for a full list of authors and acknowledgements.