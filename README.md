### This repository contains the resources for the Caenorhabditis becei (CBCI) RIL panels
### and corresponds to this article: XXXXXXX

"genotypes/" contains the genotypes of RILs that have been resolved with the haplotypes (see "suppl/HaplotypeReconstruction/") + their pruned version + some meta data

"haplotypes/" contains the founders' haplotype blocks inferred for each RIL

"phenotypes/" contains the sex-specific size and growth rate measurements

"figures/" contains the figures found in the article 

"suppl/" contains the resources regarding the haplotype reconstruction and the prediction of sex for sex-specific size



CBCI_RILs
├── README.md
├── analysis
│   ├── GenomicCharacterization
│   │   ├── LD_PCA_AFS.R
│   │   ├── cumulativeBreaks_mareyMaps.R
│   │   ├── foundersHomology.R
│   │   ├── linkage_by_geneticDistance.csv
│   │   └── linkage_by_physicalDistance.csv
│   ├── QTLmapping
│   │   ├── EMMAX_functions.R
│   │   ├── GWAS_EMMAX_LOCO_SizeSexConv.Rdata
│   │   ├── GWAS_EMMAX_LOCO_SizeSexDiv.Rdata
│   │   ├── GWAS_EMMAX_LOCO_growthrate.Rdata
│   │   ├── correlation_QTL_AFC.R
│   │   ├── growthate_EMMAX_LOCO.R
│   │   ├── growthrate_haplotype_GWAS_LOCO.R
│   │   ├── haplotypeGWAS_sizeconv_results_LOCO.Rdata
│   │   ├── haplotypeGWAS_sizediv_results_LOCO.Rdata
│   │   ├── power_analysis
│   │   │   └── qtl_mapping_power_simulations.R
│   │   ├── size_EMMAX_LOCO.R
│   │   └── size_haplotype_GWAS_LOCO.R
│   ├── heritability
│   │   ├── growthrate_gvariance_partition_chrom.csv
│   │   ├── growthrate_gvariance_partition_recDomains.csv
│   │   ├── haplotype_relatedness_matrix
│   │   │   ├── Ghap_haplotypeRelatednessMatrix.csv
│   │   │   ├── chrIII_Ghap.csv
│   │   │   ├── chrII_Ghap.csv
│   │   │   ├── chrIV_Ghap.csv
│   │   │   ├── chrI_Ghap.csv
│   │   │   ├── chrV_Ghap.csv
│   │   │   ├── chrX_Ghap.csv
│   │   │   └── haplotype_relatadness_matrix.R
│   │   ├── sizeConv_gvariance_partition_chrom.csv
│   │   ├── sizeConv_gvariance_partition_recDomains.csv
│   │   ├── sizeDiv_gvariance_partition_chrom.csv
│   │   ├── sizeDiv_gvariance_partition_recDomains.csv
│   │   └── variancePartitioning&heritability.R
│   └── simu_panel_deriv
│       ├── observed_summarydata.Rdata
│       ├── scripts
│       │   ├── becei_panels_deriv_neutral.slim
│       │   └── simu_panel_deriv.R
│       ├── simu_recombinationmaps.csv.gz
│       └── simu_summary.csv
├── figures
│   ├── Fig_mosaic.png
│   ├── Fig_power.png
│   ├── SFig_LD.png
│   ├── SFig_LD_domain.png
│   ├── SFig_founderHomology.png
│   ├── SFig_powerSinglePanel.png
│   ├── SFig_simu_deriv.png
│   ├── fig_gwas_haplotype_sizeConv+Div_LOCO.pdf
│   ├── fig_gwas_haplotype_sizeConv+Div_LOCO.png
│   ├── fig_variancePartition.png
│   ├── ~$figures_becei_rils.pptx
│   └── ~$founder_IBD_plot.pptx
├── genotypes
│   ├── III_becei_genotypes_RILs.csv.gz
│   ├── III_becei_genotypes_founders.csv.gz
│   ├── III_becei_variantInfo_founders&Rils.csv.gz
│   ├── II_becei_genotypes_RILs.csv.gz
│   ├── II_becei_genotypes_founders.csv.gz
│   ├── II_becei_variantInfo_founders&Rils.csv.gz
│   ├── IV_becei_genotypes_RILs.csv.gz
│   ├── IV_becei_genotypes_founders.csv.gz
│   ├── IV_becei_variantInfo_founders&Rils.csv.gz
│   ├── I_becei_genotypes_RILs.csv.gz
│   ├── I_becei_genotypes_founders.csv.gz
│   ├── I_becei_variantInfo_founders&Rils.csv.gz
│   ├── V_becei_genotypes_RILs.csv.gz
│   ├── V_becei_genotypes_founders.csv.gz
│   ├── V_becei_variantInfo_founders&Rils.csv.gz
│   ├── X_becei_genotypes_RILs.csv.gz
│   ├── X_becei_genotypes_founders.csv.gz
│   ├── X_becei_variantInfo_founders&Rils.csv.gz
│   ├── beceiPanels_geno_RILs_pruned0.999.csv.gz
│   └── beceiPanels_variantsInfo_pruned0.999.csv.gz
├── haplotypes
│   ├── IBD_between_founding_haploid_genomes.csv
│   ├── III_rils_foundingHaplotypesBlocks.csv
│   ├── II_rils_foundingHaplotypesBlocks.csv
│   ├── IV_rils_foundingHaplotypesBlocks.csv
│   ├── I_rils_foundingHaplotypesBlocks.csv
│   ├── README.md
│   ├── V_rils_foundingHaplotypesBlocks.csv
│   ├── X_rils_foundingHaplotypesBlocks.csv
│   └── recombination_breakpoints.csv
├── phenotypes
│   ├── README.md
│   ├── growthrates.csv
│   ├── logGrowthrates_emmeans.csv
│   ├── size_data_raw.csv.gz
│   ├── size_emmeans.csv
│   └── size_summarystat.csv
├── suppl
│   ├── ImageXpress_sex_xgboost
│   │   ├── models
│   │   │   ├── sexMale_xgb_preds
│   │   │   ├── sexMale_xgb_preds.feature_names
│   │   │   ├── wrongObjects_xgb_preds
│   │   │   └── wrongObjects_xgb_preds.feature_names
│   │   ├── predictSex.R
│   │   └── training_datasets
│   │       ├── attributeSexFromFluo.R
│   │       ├── imageXpress_sexTraining_fluoAttributed.xlsx
│   │       ├── imageXpress_sexTraining_humanAttributed.xlsx
│   │       ├── imageXpress_sexTraining_humanAttributed2.xlsx
│   │       └── train_sex_XGB.R
│   ├── README.md
│   ├── RILs_sequencing_metadata.csv
│   ├── becei_domainsAndRates.csv
│   ├── cbecei_geneticMap.txt
│   ├── genotype_stringentFiltering
│   │   ├── III_geno_becei_RILs_stringentFiltering.csv.gz
│   │   ├── III_snps_becei_RILs&Pools_stringentFiltering.csv.gz
│   │   ├── II_geno_becei_RILs_stringentFiltering.csv.gz
│   │   ├── II_snps_becei_RILs&Pools_stringentFiltering.csv.gz
│   │   ├── IV_geno_becei_RILs_stringentFiltering.csv.gz
│   │   ├── IV_snps_becei_RILs&Pools_stringentFiltering.csv.gz
│   │   ├── I_geno_becei_RILs_stringentFiltering.csv.gz
│   │   ├── I_snps_becei_RILs&Pools_stringentFiltering.csv.gz
│   │   ├── V_geno_becei_RILs_stringentFiltering.csv.gz
│   │   ├── V_snps_becei_RILs&Pools_stringentFiltering.csv.gz
│   │   ├── X_geno_becei_RILs_stringentFiltering.csv.gz
│   │   └── X_snps_becei_RILs&Pools_stringentFiltering.csv.gz
│   ├── haplotypeReconstruction
│   │   ├── 00_beceiFounders_phasing.R
│   │   ├── 01_haploCall_RILs.R
│   │   ├── 02_filterHaplotypeConsistentVariants.R
│   │   ├── 03_haploCall_RILs2.R.R
│   │   ├── 04_inpute_RILs_missing_SNPs.R
│   │   └── functions
│   │       ├── functions_beceiFounders_phasing.R
│   │       └── functions_diplosearch.R
│   ├── pools
│   │   ├── FounderGenotype_logLikelihood.R
│   │   ├── III_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz
│   │   ├── III_founder_genotype_loglikelihood.csv.gz
│   │   ├── II_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz
│   │   ├── II_founder_genotype_loglikelihood.csv.gz
│   │   ├── IV_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz
│   │   ├── IV_founder_genotype_loglikelihood.csv.gz
│   │   ├── I_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz
│   │   ├── I_founder_genotype_loglikelihood.csv.gz
│   │   ├── V_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz
│   │   ├── V_founder_genotype_loglikelihood.csv.gz
│   │   ├── X_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz
│   │   └── X_founder_genotype_loglikelihood.csv.gz
│   └── variant_filtering_metadata.xlsx
└── utils.R
