### This repository contains the resources for the Caenorhabditis becei (CBCI) RIL panels
### and corresponds to this article: XXXXXXX

`genotypes/` contains the genotypes of RILs that have been resolved with the haplotypes (see "suppl/HaplotypeReconstruction/") + their pruned version 

`haplotypes/` contains the founders' haplotype blocks inferred for each RIL

`phenotypes/` contains the sex-specific size and growth rate measurements

`figures/` contains the figures found in the article 

`suppl/` contains the resources regarding the haplotype reconstruction and the prediction of sex for sex-specific size

# Repository structure:

```text
CBCI_RILs
├── analysis
│   ├── GenomicCharacterization
│   ├── QTLmapping
│   │   └── power_analysis
│   ├── heritability
│   │   └── haplotype_relatedness_matrix
│   └── simu_panel_deriv
│       └── scripts
├── figures
├── genotypes
├── haplotypes
├── phenotypes
└── suppl
    ├── ImageXpress_sex_xgboost
    │   ├── models
    │   └── training_datasets
    ├── genotype_stringentFiltering
    ├── haplotypeReconstruction
    │   └── functions
    └── pools
```
