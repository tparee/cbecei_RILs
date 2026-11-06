library(readr)
library(data.table)
library(ggplot2)
library(sommer)

# Path
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/")
source("utils.R")

# Function to generate kinship matrix
get.K_ASV = function(GT){
  ASV = scale(GT,center=T,scale=F) %*% t(scale(GT,center=T,scale=F))
  ASV = ASV / (psych::tr(ASV)/(nrow(GT)-1))
  ASV
}


# Meta data file
meta = read_csv("suppl/RILs_sequencing_metadata.csv")

# load GWAS results
load(file = "analysis/QTLmapping/GWAS_EMMAX_LOCO_growthrate.Rdata")

# Load neutral distribution of p-values obtained by permutation
permutated_pvals = unlist(read.table("analysis/temp/permutated_pval_growthrate.txt"))
# Define the alpha error = 0.05 threshold
threshold_alpha0.05 = quantile( permutated_pvals , prob = 0.05) # 7.099286 

# Identify our four QTLs
qtl1 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'I' & res_LOCO$cm < 40 & res_LOCO$cm > 20,]
qtl1 = qtl1$id[which.min(qtl1$pval)]
qtl2 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'I' & res_LOCO$cm > 42,]
qtl2 = qtl2$id[which.min(qtl2$pval)]
qtl3 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'II' & res_LOCO$cm > 45,]
qtl3 = qtl3$id[which.min(qtl3$pval)]
qtl4 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'X' & res_LOCO$cm > 20,]
qtl4 = qtl4$id[which.min(qtl4$pval)]

# Load growthrate data
growthrates = read.csv( "phenotypes/growthrates.csv", sep = " ")
growthrates$strain = tstrsplit(growthrates$strain,"_")[[2]] # Remove the first part of the strain idenfifier to match genotype table naming convention
growthrates = subset(growthrates, strain %in% meta$rilname[meta$panel == "alpha"]) # keep panel alpha

# Load snps info and genotype matrix
snps = as.data.frame(read_csv("genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))
genotypes <- as.matrix(fread("genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
genotypes = (genotypes-0.5)*2

# Keep only strain that have both phenotype and genotype data
growthrates = subset(growthrates, strain %in% colnames(genotypes))
genotypes = genotypes[,colnames(genotypes) %in% growthrates$strain]

# remove fixed snps if any
af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})
genotypes = genotypes[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]

#  Add QTLs genotype 
growthrates$qtl1 = genotypes[snps$id==qtl1,growthrates$strain]
growthrates$qtl2 = genotypes[snps$id==qtl2,growthrates$strain]
growthrates$qtl3 = genotypes[snps$id==qtl3,growthrates$strain]
growthrates$qtl4 = genotypes[snps$id==qtl4,growthrates$strain]

# Scale QTL genotype
# Important that variance is one so random effect values are comparable
# see https://doi.org/10.1371/journal.pcbi.1009923
growthrates$qtl1 = scale(growthrates$qtl1)
growthrates$qtl2 = scale(growthrates$qtl2)
growthrates$qtl3 = scale(growthrates$qtl3)
growthrates$qtl4 = scale(growthrates$qtl4)

# Create a unique identifier per strain 
growthrates$id = paste0(growthrates$strain,'_', 1:nrow(growthrates))

# Generate our marker matrix (n x m) that match observation ID
GT = t(genotypes) # m x n => n x m
GT = GT[match(growthrates$strain,row.names(GT)),]
row.names(GT) = growthrates$id
GT[is.na(GT)]=0

# Build the kinship matrix
K=get.K_ASV(GT)

# For each QTL, estimate variance: random effect + fixed effect approaches
thisqtl = 'qtl1'
dt = growthrates
dt$focalqtl = dt[,thisqtl]

# Model focal qtl as a fixed effect
mod_qtl_fixed <- mmer(
  log_hours_to_starve ~ focalqtl,
  random = ~  block + vs(id, Gu = K),
  data = dt
)

X <- model.matrix(~ focalqtl, data = dt) 
b <- mod_qtl_fixed$Beta[,3] # regression coeficients
Vq_fixed <- var(drop(X %*% b)) # Variance of the QTL effect estimates
Vk_fixed <- summary(mod_qtl_fixed)$varcomp[2,1] # Variance explained by the kinship matrix K
PVE_fixed <- Vq_fixed/(Vq_fixed + Vk_fixed) # Variance proportions 

# Model focal qtl as a fixed effect
Kqtl <- tcrossprod(dt$focalqtl) # "QTL kinship" matrix
dimnames(Kqtl) <- list(as.character(dt$id), as.character(dt$id))
dt$id_qtl <- factor(dt$id)

modh2_qtl_random <- mmer(
  log_hours_to_starve ~ 1,
  random = ~ vs(id_qtl, Gu = Kqtl) + block + vs(id, Gu = K),
  data = dt
)

# Variance components
vc <- summary(modh2_qtl_random)$varcomp
V_QTL <- vc[1,1]
V_K <- vc[3,1]

PVE_random = V_QTL/(V_QTL + V_K) # Variance proportions 

c(PVE_fixed=PVE_fixed,  PVE_random=PVE_random)
# qtl1: 0.2433925  0.2346618
# qtl2: 0.05641927 0.05163050 
# qtl3:  0.1277792  0.1218484
# qtl4: 0.09431978 0.08475154 


s


