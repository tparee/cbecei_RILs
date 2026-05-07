# Haplotypes Directory

---

## `IBD_between_founding_haploid_genomes.csv`

Contains genomic segments that are identical-by-descent (IBD) between pairs of founding haplotypes from `genotypes/[N]_becei_genotypes_founders.csv.gz`. Only segments with more than 20 variants, larger than 50 kb and Larger than 0.5 cM were retained.

These regions were identified using: `analysis/GenomicCharacterization/foundersHomology.R`

---

## `[N]_rils_foundingHaplotypesBlocks.csv`

`[N]` corresponds to the chromosome: I, II, III, IV, V, X. This table contains inferred founder identities across the chromosome for each RIL.

Each row corresponds to a haplotype block spanning: ( `pos1` → start position;  `pos2` → end position)

The corresponding founder assignment is stored in: `founder`. Multiple founders may be assigned to a block when founders cannot be distinguished within that genomic interval due to IBD between founding haplotypes.
Multiple founders can also occur when a heterozygous diplotype is detected `ishet == TRUE`. These correspond to genomic intervals that remain heterozygous after 25 generations of sib-mating and whose observed genotypes are well explained by a combination of founding haplotypes.

Within the `founder` column:`&` separates the two haplotypes of a heterozygous diplotype, while `;` separates multiple possible founder identities that cannot be distinguished because of IBD
Example: `FA.g1&FM.g1;FM.g2` cooresponds to first haplotype being `FA.g1` and second haplotype: either `FM.g1` or `FM.g2`.

---

## `recombination_breakpoints.csv`

Contains genomic intervals corresponding to possible recombination junctions between adjacent haplotype blocks in `[N]_rils_foundingHaplotypesBlocks.csv`.
This file was used to generate the Marey maps shown in Figure 1A.
