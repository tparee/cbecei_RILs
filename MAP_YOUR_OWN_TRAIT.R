library(sommer)
library(readr)
library(data.table)
library(ggplot2)

# Change the path to your own computer path
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/")

# Load EMMAX mapping function
source("analysis/QTLmapping/EMMAX_functions.R")

# Kinship matrix building function
# See Average semi-variance method: https://doi.org/10.1093/g3journal/jkac080
get.K_ASV = function(GT){
  ASV = scale(GT,center=T,scale=F) %*% t(scale(GT,center=T,scale=F))
  ASV = ASV / (psych::tr(ASV)/(nrow(GT)-1))
  ASV
}

# Load snps information (chrom, pos, ref, alt, qual, etc.)
snps = as.data.frame(read_csv("genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))
# Load pruned genotype (you can use unpruned matrix but mapping will be longer doing redundant association tests)
genotypes <- as.matrix(fread("genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
# Center the genotype matrix: [0,1] to [-1,1]
genotypes = (genotypes-0.5)*2 


SIMULATE_PHENOTYPE = T

if(!SIMULATE_PHENOTYPE){
  
  # Your trait measurments must me named "trait_value" (if not, change the EMMAX function accordingly)
  # You must have a phenotypic observation "ID" (different ID per row)
  # You must have a "strain" column identify the strain/RIL 
  # pheno = read.csv( "my_phenotypic_data.csv")
  
}else{
  set.seed(1234)
  
  VP = 100 # phenotypic variance
  VK = VP*0.4 # Variance explained by the polygenic background
  VQTL = VP*0.1 # Variance explained by the QTL
  VE = VP-VK-VQTL # Environmental variance
  
  n_replicate_phenotypic_measurments = 3
  nloci = 1000 # number of loci
  
  causal_loci = sample(1:nrow(genotypes), nloci) # sample loci
  qtl = causal_loci[1] # The QTL
  polygenic_background = sort(causal_loci[-1])
  polygenic_background_effects = rnorm(nloci-1) # sample effect from a normal distribution
  
  # polygenic background genotypic value for each genotype
  g_poly = apply(polygenic_background_effects*genotypes[polygenic_background,],2,sum,na.rm=T) 
  g_poly = c(scale(g_poly))*sqrt(VK)  # scale the variance to match VK
  g_qtl = c(scale(genotypes[qtl,]))*sqrt(VQTL)
  
  gtot = g_poly+g_qtl
  
  nobs = ncol(genotypes)*n_replicate_phenotypic_measurments
  pheno = data.frame(ID = factor(as.character(1:nobs)), # observation unique identifier
                     strain = rep(colnames(genotypes), each = 3), # strain identifier
                     trait_value = rep(gtot, each = 3) + rnorm(nobs,sd = sqrt(VE)) # phenotypic measurment: P = G + E
                     )
  
  #lvls = aggregate(trait_value~strain, pheno,mean)
  #lvls = lvls[order(lvls[,2]),1]
  #ggplot(pheno, aes(factor(strain, levels = lvls), trait_value))+
  #  theme_classic()+
  #  geom_point(color = 'grey')+
  #  stat_summary(size = 0.2)+
  #  theme(axis.text.x = element_blank())
}

#####################################################################
################  heritability and GWAS, NOT LOCO  ##################

# Keep only strain that both have phnotype and genotype
pheno = subset(pheno, strain %in% colnames(genotypes))
genotypes = genotypes[,colnames(genotypes) %in% pheno$strain]

# Remove fixed SNPs, if any (can happens when only a subset of strains is used)
af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})
genotypes = genotypes[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]

if("ID" %in% colnames(pheno)){
  warning("Make sure that ID is a unique row identifier and is a factor")
}else{
  pheno$ID = factor(as.character(1:nrow(pheno)))
}

# Build the lines x marker matrix with each row matching the pheno data.frame
GT = t(genotypes)
GT = GT[match(pheno$strain,row.names(GT)),]
row.names(GT) = pheno$ID
GT[is.na(GT)]=0

########################################################
###################  Heritability  #####################

K=get.K_ASV(GT)
colnames(K)
# Heritability 
modh2 <- mmer(
  trait_value ~ 1,
  random = ~  vs(ID, Gu = K),
  data = pheno
)

vcomp = summary(modh2)$varcomp[1:2,1]
h2 = vcomp[1]/sum(vcomp)
h2 

delta = vcomp[2]/vcomp[1]


########################################################
###################  GWAS  #############################

emmax_res = emmax_assoc(y=pheno$trait_value,G=GT,
            covar = NULL, kinship = K,
            delta = delta)

emmax_res = cbind(snps, emmax_res)

ggplot(emmax_res, aes(pos/1e6,-log10(pval)))+
  theme_classic()+
  facet_wrap(~chrom,nrow=1)+
  geom_point()+
  geom_point(data=emmax_res[ifelse(SIMULATE_PHENOTYPE,qtl,NA),],color = 'orange')+ # true causal SNP
  xlab("Physical position (Mb)")





#####################################################################
################  heritability and GWAS, LOCO method  ###############
################# (Leave One Chromosome Out)          ###############

# => for each chromosome, calculate K while excluding it and run gwas

# emmax_res_loco = NULL
# for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
#   print(CHR)
#   KCHR = get.K_ASV(GT[,snps$chrom != CHR])
#   GTCHR = GT[,snps$chrom == CHR]
#   
#   modh2 <- mmer(
#     trait_value ~ 1,
#     random = ~  vsr(ID, Gu = KCHR),
#     data = pheno
#   )
#   
#   deltachr = vcomp[2]/vcomp[1]
#   
#   emmax_res_chr = emmax_assoc(y=pheno$trait_value,G=GTCHR,
#                           covar = NULL, kinship = KCHR,
#                           delta = deltachr, verbose = F)
#   
#   
#   emmax_res_chr = cbind(snps[snps$chrom == CHR,], emmax_res_chr)
#   emmax_res_loco = rbind(emmax_res_loco,emmax_res_chr)
# }
# 
# 
# ggplot(emmax_res_loco, aes(pos/1e6,-log10(pval)))+
#   theme_classic()+
#   facet_wrap(~chrom,nrow=1)+
#   geom_point()+
#   geom_point(data=emmax_res[ifelse(SIMULATE_PHENOTYPE,qtl,NA),],color = 'orange')+ # true causal SNP
#   xlab("Physical position (Mb)")
