library(emmeans)
library(lme4)
library(sommer)
library(readr)
library(data.table)
library(ggplot2)
library(ggplot2)
library(gridExtra)
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/")
meta = read_csv("suppl/RILs_sequencing_metadata.csv")
source("utils.R")
source("analysis/QTLmapping/EMMAX_functions.R")

get.K_ASV = function(GT){
  ASV = scale(GT,center=T,scale=F) %*% t(scale(GT,center=T,scale=F))
  ASV = ASV / (psych::tr(ASV)/(nrow(GT)-1))
  ASV
}



###################################################################
################### GWAS with pruned genotypes ####################
###################################################################

growthrates = read.csv( "phenotypes/growthrates.csv", sep = " ")
growthrates$log_hours_to_starve = log(growthrates$hours_to_starve)
growthrates$strain = tstrsplit(growthrates$strain,"_")[[2]]
growthrates = subset(growthrates, strain %in% meta$rilname[meta$panel == "α"])

snps = as.data.frame(read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))
genotypes <- as.matrix(fread("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
genotypes = (genotypes-0.5)*2
#genotypes[genotypes == 0]=1

growthrates = subset(growthrates, strain %in% colnames(genotypes))
genotypes = genotypes[,colnames(genotypes) %in% growthrates$strain]

af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})
genotypes = genotypes[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]

growthrates$id = paste0(growthrates$strain,'_', 1:nrow(growthrates))
GT = t(genotypes)
GT = GT[match(growthrates$strain,row.names(GT)),]
row.names(GT) = growthrates$id
GT[is.na(GT)]=0

#growthrates$id = factor(growthrates$strain, levels = rownames(GT))
K=get.K_ASV(GT)

modh2 <- mmer(
  log_hours_to_starve ~ block,
  random = ~  block + vs(id, Gu = K),
  data = growthrates
)

summary(modh2)$varcomp[1:2,1]/sum(summary(modh2)$varcomp[1:2,1])

# block:0.1340243, K: 0.8659757

#############################################################################
############### SNP by EMMAX + LOCO #########################################
#############################################################################

fixedBlock = F
perm = T
if(perm==T){n=1000; permutated_pvals = matrix(NA,ncol=n,nrow=6); set.seed(123)}else{n=1}

if(fixedBlock){fixed_effect_modelMatrix = model.matrix(~ block, growthrates)}

KCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)
GTCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)
deltas.list =list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)

for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
  KCHR = get.K_ASV(GT[,snps$chrom != CHR])
  GTCHR = GT[,snps$chrom == CHR]
  
  if(!fixedBlock){
    modh2 <- mmer(
      log_hours_to_starve ~ 1,
      random = ~  block + vsr(id, Gu = KCHR),
      data = growthrates
    )
    
    vcomp = summary(modh2)$varcomp[1:3,1]
    IndicatorMatrix_block = model.matrix(~ block-1, growthrates)
    vcov_block = crossprod( t(IndicatorMatrix_block))
    KCHR = (vcomp[2]*KCHR + vcomp[1]*vcov_block)/sum(vcomp[1:2])
    deltachr = vcomp[3]/sum(vcomp[1:2])
    
  }else{
    modh2 <- mmer(
      log_hours_to_starve ~ as.factor(block),
      random = ~ vsr(id, Gu = KCHR),
      data = growthrates
    )
    
    vcomp = summary(modh2)$varcomp[1:32,1]
    deltachr = vcomp[2]/vcomp[1]
  }
  
  KCHR.list[[CHR]] = KCHR
  GTCHR.list[[CHR]] = GTCHR
  deltas.list[[CHR]] = deltachr

}

for(i in 1:n){
  if(perm == T & i %% 100 == 0){print(i)}
  res_LOCO = NULL
  for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
    KCHR = KCHR.list[[CHR]]
    GTCHR = GTCHR.list[[CHR]]
    deltachr = deltas.list[[CHR]]
    
    if(perm==T){
      permorder = match(growthrates$strain, sample(unique(growthrates$strain)))
      GTCHR = GTCHR[permorder,]
      }
    
    if(fixedBlock){
      reschr = emmax_assoc(y=growthrates$log_hours_to_starve,G=GTCHR,
                           covar = fixed_effect_modelMatrix, kinship = KCHR,
                           delta = deltachr, verbose = !perm)
    }else{
      reschr = emmax_assoc(y=growthrates$log_hours_to_starve,
                           G=GTCHR, covar = NULL, kinship = KCHR,
                           delta = deltachr, verbose = !perm)
    }
    
    
    reschr = cbind(snps[snps$chrom == CHR,], reschr)
    res_LOCO = rbind(res_LOCO,reschr)
  }
  
  if(perm==T){
    minp = min(res_LOCO$pval, na.rm = T)
    outfile = "analysis/temp/permutated_pval.txt"
    cat(minp, file = outfile, append = TRUE, sep = "\n")
  }
}


ggplot(reschr, aes(cm, -log10(pval)))+geom_point()


threshold_alpha0.05 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.05)) # 6.450801 
threshold_alpha0.1 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.1)) # 5.898521 
#save(permutated_pvals, file= "analysis/permutated_pvals_growthrate.Rdata")

threshold_alpha0.05 =  6.450801 

ggplot()+
  #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, growthrates=1, color = 'blue')+
  geom_point(data = res_LOCO, aes(cm, -log10(p)))+
  facet_grid(.~chrom, scales='free_x')+theme_classic()+
  geom_hline(yintercept = threshold_alpha0.05)


ggplot()+
  #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, growthrates=1, color = 'blue')+
  geom_point(data = res_LOCO, aes(cm, -log10(p)), size=0.3)+
  facet_grid(.~chrom, scales='free_x')+theme_classic()+
  geom_hline(yintercept = threshold_alpha0.05)


#save(res_LOCO, file = "analysis/QTLmapping/GWAS_EMMAX_LOCO_growthrate.Rdata")

load("analysis/QTLmapping/GWAS_EMMAX_LOCO_growthrate.Rdata")


pgwas = ggplot(res_LOCO)+
  facet_grid(.~chrom, scales='free_x')+theme_Publication3()+
  coord_cartesian(ylim = c(0,-log10(min(res_LOCO$p))+0.6))+
  geom_hline(yintercept = threshold_alpha0.05, color = "#CC6600")+
  #geom_hline(yintercept = threshold_alpha0.1, color = "#0073B3")+
  #geom_rect(data=data.frame(chrom="I"), aes(xmin = linked[1], xmax=linked[2],ymin=-Inf, ymax = Inf), fill = "yellow", color=NA, alpha=0.5)+
  geom_point(alpha=1, size=0.25,  aes(cm, -log10(p)))+
  geom_point( data = data.frame(chrom = unique(res_LOCO$chrom), x=0,y=0), aes(x,y),color = NA)+
  #geom_point(data = subset(snps2, -log10(pval) > threshold_alpha0.05), alpha=1, size=1.2,  aes(cm, -log10(pval)), fill = "red", shape = 21, color='black')+
  theme(panel.spacing = unit(0, "lines"))+
  scale_x_continuous(breaks = c(0,15,30, 45))+
  theme(legend.position = "none")+
  xlab("Genetic distance (cM)")


ggsave(pgwas, file="figures/fig_gwas_SNP_growthrate_LOCO.png", width=3.6, height=1.15, dpi=1200)
ggsave(pgwas, file="figures/fig_gwas_SNP_growthrate_LOCO.pdf", width=3.6, height=1.15, dpi=1200)


# Fixed block = F
#############################################################################
############### SNP by EMMAX + LOCO #########################################
#############################################################################
perm = F

if(perm==T){n=1000; permutated_pvals = matrix(NA,ncol=n,nrow=6)}else{n=1}
for(i in 1:n){
  if(perm == T & i %% 100 == 0){print(i)}
  res_LOCO = NULL
  for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
    KCHR = get.K_ASV(GT[,snps$chrom != CHR])
    GTCHR = GT[,snps$chrom == CHR]
    if(perm==T){GTCHR = GTCHR[sample(1:nrow(GTCHR), nrow(GTCHR)),] }
    reschr = emmax_assoc(y=growthrates$emmean,G=GTCHR, covar = NULL, kinship = KCHR, verbose = FALSE)
    reschr = reschr$results
    reschr = cbind(snps[snps$chrom == CHR,], reschr[,-1])
    res_LOCO = rbind(res_LOCO,reschr)
  }
  
  if(perm==T){pperm_min = aggregate(p~chrom,res_LOCO,min,na.rm=T); permutated_pvals[,i] = pperm_min$p}
}


#-log10(apply(permutated_pvals, 1, function(x){quantile(x, prob = 0.05)}))

threshold_alpha0.05 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.05))
threshold_alpha0.1 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.1))

# threshold_alpha0.1 = 3.793954
# threshold_alpha0.05 = 4.106768
ggplot()+
  #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, growthrates=1, color = 'blue')+
  geom_point(data = res_LOCO, aes(cm, -log10(p)),growthrates=1)+
  facet_grid(.~chrom, scales='free_x')+theme_classic()+geom_hline(yintercept = threshold_alpha0.05)



#############################################################################
############### SNP by EMMAX + LOCO #########################################
#############################################################################

for(i in 1:2){
  
  if(i == 1){
    GTw = GT
  }else{
    GTw = t(t(GT)*res_LOCOw$PP)
  }
  res_LOCOw = NULL
  for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
    KCHR = get.K_ASV(GTw[,snps$chrom != CHR])
    reschr = emmax_assoc(y=growthrates$emmean,G=GT[,snps$chrom == CHR], covar = NULL, kinship = KCHR, verbose = TRUE)
    reschr = reschr$results
    reschr = cbind(snps[snps$chrom == CHR,], reschr[,-1])
    res_LOCOw = rbind(res_LOCOw,reschr)
  }
  
  pi = 0.05
  res_LOCOw$PP = (pi*exp(res_LOCOw$LR))/(pi*exp(res_LOCOw$LR) + (1-pi))
  
}


ggplot()+
  geom_point(data = res_LOCO, aes(cm, -log10(p)), alpha=0.3, growthrates=1, color = 'blue')+
  geom_point(data = res_LOCOw, aes(cm, -log10(p)), alpha=0.3, growthrates=1, color = 'red')+
  facet_grid(.~chrom, scales='free_x')+theme_classic()+geom_hline(yintercept = 4)#+ylim(0,4)





