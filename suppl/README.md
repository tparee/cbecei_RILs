# Supplementary Directory

This supplementary directory contains pipelines used to:

1. Jointly infer the founding haplotypes and filter SNPs.
2. Measure sex-specific body size from ImageXpress images.

---

# 1. Joint Inference of Founding Haplotypes and SNP Filtering

The `.fastq` files were aligned to the *C. becei* QG2082 reference genome  
(DOI: https://doi.org/10.1101/2025.05.09.653148).

Variants were called using `bcftools mpileup`, retaining variants with:

- `QUAL > 20`

The selected genotype for each site is the genotype minimizing the PL value (Phred-scaled likelihood). Genotypes with:

- `GQ < 10` (genotype likelihood < 0.9)

were set to `NA`.

Because the RILs are mosaics of recombined founding haplotypes, variant filtering was performed by reconstructing the haplotype block structure and retaining only SNPs consistent with this structure.

Founding haplotypes were reconstructed using:

- Pool-sequencing data from F1 progeny of founding crosses (`crossPool`)
- Linkage information from the RILs

An initial backbone of highly confident SNPs was first generated (see `genotype_stringentFiltering`).

For these SNPs, the log-likelihood of observing the allelic depths in the CrossPools was computed for every possible combination of three founder biallelic genotypes (see `pools`).

Founding haplotypes were then estimated by selecting combinations of RIL haplotype blocks that maximize the founder genotype log-likelihood (see `haplotypeReconstruction`).

Finally, founder identity across the RIL genomes was inferred. SNPs consistent with the inferred haplotype mosaic were retained, while inconsistent SNPs were filtered out.

The number of SNPs retained at each filtering step is reported in:

```text
XXXX.csv
```

---

# `genotype_stringentFiltering`

## Contents

### `[N]_genotypes_becei_RILs_stringentFiltering.csv`

RIL genotype table:

```text
[SNPs × RILs]
```

generated using highly stringent filtering criteria.

### `[N]_genotypes_snps_RILs_stringentFiltering.csv`

Corresponding SNP annotation table.

### Chromosome notation

`[N]` corresponds to the chromosome:

```text
I, II, III, IV, V, X
```

## Stringent Filtering Parameters

### Variant-level filtering

```text
QUAL == 999
MQ >= 59
DP > 2e4
DP < 3e4
MQSB > 0.001
VDB > 0.001
RPB > 0.001
```

### Genotype frequency filtering

Removes:

- Fixed variants
- Frequently missing variants
- Suspiciously heterozygous variants

Criteria:

```text
n_miss > 10
n_het < 10
n_ref > 0
n_alt > 0
```

### CrossPool depth filtering

```text
DP_crossApool >= 10
DP_crossBpool >= 10
DP_crossCpool >= 10
```

## Purpose

These tables are intended for construction of the initial backbone of founding haplotypes.

---

# `pools`

## Contents

### `[N]_AllelicDepth_becei_CrossPools_stringentFiltering.csv`

Contains:

- SNP annotations
- Allelic depth information (`ref`, `alt`)
- Pool and population sequencing data

## Pool Definitions

- `CrossApool`: F1 progeny of Cross A (`Founder A × Founder M`)
- `CrossBpool`: F1 progeny of Cross B (`Founder B × Founder M`)
- `CrossCpool`: F1 progeny of `Founder M × QG2082` (reference strain)

`POPA` and `POPB` correspond to populations A and B after five generations of outcrossing. Each population was sequenced in triplicate.

The SNP set is identical to the one in `genotype_stringentFiltering`.

---

### `[N]_founder_genotype_loglikelihood.csv`

Contains the log-likelihood of observing allelic depths in F1 pools for every possible founder diploid genotype combination.

## Founder Genotype Encoding

Founder genotype scenarios are encoded in the column names.

Example:

```text
FA:0.5;FB:0;FM:1
```

corresponds to:

- Founder A: heterozygous
- Founder B: homozygous reference
- Founder M: homozygous alternative

---

### `FounderGenotype_logLikelihood.R`

R script used to compute founder genotype log-likelihoods.

---

# `haplotypeReconstruction`

## `00_beceiFounders_phasing.R`

Infers the six founding haplotypes using the stringent SNP set.

### Method Overview

- The genome is analyzed using overlapping windows of 500 SNPs.
- Common haplotypes in the RILs are inferred.
- Haplotypes are assigned to founders by maximizing the founder genotype log-likelihood estimated from pool sequencing data.

### Phasing Strategy

When heterozygous segments are interrupted by homozygous regions, phasing is selected to minimize the number of recombination breakpoints across RILs.

---

## `01_haploCall_RILs.R`

Calls inferred founding haplotype blocks in RILs.

### Method Overview

- Multiple genomic windows of variable size are generated for each RIL.
- If present, a founding haplotype matching at ≥99% similarity is assigned.
- Haplotype boundaries are refined by selecting breakpoints minimizing the Hamming distance between inferred haplotypes and observed genotypes.

---

## `02_filterHaplotypeConsistentVariants`

Checks compatibility between each SNP (`QUAL > 20`) and the inferred haplotype structure.

### Compatibility Criterion

A SNP is considered compatible if it is monomorphic within each founding haplotype.

Otherwise, the SNP would imply more than six founding haplotypes, which is inconsistent with the dataset design.

---

## `03_haploCall_RILs2.R`

Recomputes haplotype blocks using the complete filtered SNP set.

The increased SNP density improves breakpoint resolution.

---

## `04_inpute_RILs_missing_SNPs`

Uses inferred founding haplotype identities across the RIL genomes to:

- Impute missing SNP genotypes
- Detect unlikely genotypes consistent with genotyping errors

---

# `ImageXpress_sex_xgboost`

Sex-specific body size was measured using ImageXpress imaging and worm detection from the Andersen lab pipeline.

Sex classification was performed using an XGBoost model trained on features extracted from the Andersen pipeline.

Training labels were assigned either:

- Manually by visual inspection
- Using strain `QG4602`, which carries a sex-specific fluorescent reporter

---

## Contents

### `training datasets/`

Contains:

- Training datasets
- Scripts used to train the XGBoost model

---

### `models/`

Contains trained XGBoost models.

---

### `predictSex.R`

R script used to predict the sex of worm objects detected by the Andersen pipeline on ImageXpress images.
