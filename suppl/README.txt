##########################
### suppl directories: ###

This supplementary directory contains pipeline used to:
1) Jointly infer the founding haplotypes and filter the SNPs.
2) Measure sex-specific size from ImageXpress pictures


***************************************************************
1) Jointly infer the founding haplotypes and filter the SNPs.
***************************************************************

The .fastq files were aligned on the C.becei QG2082 reference genome (doi: https://doi.org/10.1101/2025.05.09.653148).
The variants were called using bcftools mpileup, keeping on variant with QUAL > 20.
The chosen genotype is the one that minimize the PL (Phred-scaled Likelihoods). Genotype with GQ < 10 (likelihood: < 0.9) were set as NA.

As the RILs are a mosaic of recombined haplotypes, the strategy to filter variants here is to infer the haplotype block structure and keep the SNPs that are consistent with it.
Reconstructing founding haplotypes was done by using pool sequencing data from F1 progeny of founding cross (crossPool) and linkage data available in RILs.
For the initial backbone of founding haplotype reconstruction, a set of ver high confidence SNP was used (see genotype_stringentFiltering)
For this set of SNP, log likelihood of observing the CrossPools allelic depth was calculated for every possible combination of three founder biallelic genotype (see pools).
Founding haplotypes were then estimated by selecting combination of RILs haplotype blocks that maximize the founder genotype log likelihood (see haplotypeReconstruction).
Finally The founder identity across RILs genome was inferred. SNP that are consitent with the inferred haplotype mosaic at kept, other are filtered out (see haplotypeReconstruction). 

The number of SNP kept at each step is reported in XXXX.csv

============================
genotype_stringentFiltering
============================

Contains:
- [N]_genotypes_becei_RILs_stringentFiltering.csv: RIls genotype table [snps x RILs] with very stringent filtering

- [N]_genotypes_snps_RILs_stringentFiltering.csv: corresponding snps information

* [N] is the chromosome (i.e., I,II,III,IV,V,X)

* Stringent filtering parameters:
- Filtering based on variant quality: QUAL == 999 & MQ >= 59 & DP > 2e4 & DP < 3e4 & MQSB > 0.001 & VDB > 0.001 & RPB > 0.001
- Filtering based on genotype frequency (remove fixed, often missing and suspiciously high heterozygous): n_miss > 10 & n_het < 10 & n_ref > 0 & n_alt > 0
- Filtering based depth in crosspool: DP_crossApool >= 10 & DP_crossBpool >= 10 & DP_crossCpool >= 10

* These tables are meant to be used to construct the backbones of founding haplotypes


=============
    pools
=============

Contains:
- [N]_AllelicDepth_becei_CrossPools_stringentFiltering.csv: snps info and allelic depth [ref,alt] for pools and populations

* Note:
CrossApool are the F1 progeny of cross A (Founder A x Founder M)
CrossBpool are the F1 progeny of cross B (Founder B x Founder M)
CrossCpool are the F1 progeny of Founder M x QG2082 (reference strains)
POPA and POPB are the population A and B after the 5 generations of outcrossing, they have been sequenced in triplicate

* The snps are the same as genotype_stringentFiltering

- [N]__founder_genotype_loglikelihood.csv: log likelihood of observing allelic depth in F1 pools for each possible founders diploid genotype.

* The scenario of founders genotype is in the column names
(ex: "FA:0.5;FB:0;FM:1 is the scenario that founder A is heterozygous, founder B is homozygous for the reference allele and founder M is homozygous for the alternative allele)

- FounderGenotype_logLikelihood.R: script used to compute the log likelihood

========================
haplotypeReconstruction
========================

- 00_beceiFounders_phasing.R: the script infer six founding haplotype for the stringentFiltering SNP set. 
* Across 500 overlapping SNP windows, the script infer common haplotype in RILs and attribute them to founder in a way that maximize the founder genotype loglikelihood (estimated from pools)
* Phase between heterozygous segment interrupted by homozygous segment is chosen so it minimize the number of recombination breakpoints in the RILs.

-01_haploCall_RILs.R: The script call inferred founding haplotype blocks in RILs. 
* Briefly, for a given RILs, many genomic windows of variable size are generated and, if it exist, the founding haplotype with a 99% match is returned.
Haplotypes limits are refined by chosing breakpoints minimizing the hamming distance between inferred haplotype and observed genotype.

-02_filterHaplotypeConsistentVariants: Compatibility with haplotype structure is checked for every called SNP with QUAL > 20/
* i.e., compatibility is means that the SNP is monomorphic within a founding haplotype (otherwise, it would increases the number of founding haplotype which is fixed at six in our dataset)

-03_haploCall_RILs2.R: Call the haplotypes blocks in RILs with all the kept SNPs, as more SNPs may help refine the breakpoints.

-04_inpute_RILs_missing_SNPs: Using inferred founding haplotype identity across RILs genome, we can impute missing SNPs and spot unlikely genotypes (genotyping error).

========================
ImageXpress_sex_xgboost
========================

Sex-specific size was measured by imaging the worms on a ImageXpress instruments and detecting them using the Andersen lab pipeline.
To differentiate sexes, a XGBoost model was trained on metrics extracted by the Andersen pipeline.
For objects in the training dataset, sexes was attributed visually by human, or using a strain with a sex-specific fluorescent reporter (QG4602)

- training datasets/: contains training datasets and the script to train the XGBoost model.
- models/: contains the trained XGBoost models
- predictSex.R: script used to predict sexes of the worm objects detected by the Andersen pipeline on ImageXpress images.



