library(emmeans)
#source("./basics_TP.R")
library(ggplot2)
library(gridExtra)
library(data.table)
library(sommer)


###########################################################################
################## Sex specific size ######################################
###########################################################################


###################################################################
################### GWAS with pruned genotypes ####################
###################################################################

get.K_ASV = function(GT){
  ASV = scale(GT,center=T,scale=F) %*% t(scale(GT,center=T,scale=F))
  ASV = ASV / (psych::tr(ASV)/(nrow(GT)-1))
  ASV
}

setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/becei_rils/")
source("utils.R")

size = as.data.frame(fread("phenotypes/size_summarystat.csv"))
#size = as.data.frame(fread("phenotypes/size_emmeans_Fixed.csv"))

size = as.data.frame(fread("phenotypes/size_summarystat.csv"))

genotypes =  as.matrix(read_csv("genotypes/beceiPanels_geno_RILs_pruned0.9999.csv"))
colnames(genotypes) = tstrsplit(colnames(genotypes), "-")[[1]]
snps =  as.data.frame(fread("genotypes/beceiPanels_snps_pruned0.9999.csv"))
genotypes = (genotypes-0.5)*2
genotypes[genotypes == 0]=1
genotypes = genotypes[,grepl("A_", colnames(genotypes)) & colnames(genotypes) %in% size$strain]

size = droplevels(subset(size, strain %in% colnames(genotypes) & !is.na(weight_rotated)))
size$weight_rotated = size$weight_rotated/mean(size$weight_rotated)
size$weight_f = size$weight_f/mean(size$weight_f)
size$id = paste0(size$strain,"_", 1:nrow(size))
afs = apply(genotypes, 1, function(x){mean(x+1,na.rm=T)/2})
genotypes = genotypes[afs > 0 & afs < 1,]
snps = snps[afs > 0 & afs < 1,]

GT = t(genotypes)
GT = GT[match(size$strain,row.names(GT)),]
row.names(GT) = size$id
GT[is.na(GT)]=0
K=get.K_ASV(GT)

modh2 <- mmer(
  conv ~ 1,
  random = ~  block + vsr(id, Gu = K),
  data = size, weights = weight_rotated
)



fixedBlock = F
perm = F
if(perm==T){n=1000; permutated_pvals = matrix(NA,ncol=n,nrow=6)}else{n=1}

fixed_effect_modelMatrix = model.matrix( ~ block, size)

KCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)
GTCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)

for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
  KCHR = get.K_ASV(GT[,snps$chrom != CHR])
  GTCHR = GT[,snps$chrom == CHR]
  
  if(!fixedBlock){
    modh2 <- mmer(
      conv ~ 1,
      random = ~  block + vsr(id, Gu = KCHR),
      data = size, weights = weight_rotated
    )
    
    var_prop = summary(modh2)$varcomp[2:3,1]/sum(summary(modh2)$varcomp[2:3,1])
    IndicatorMatrix_block = model.matrix(~ block-1, size)
    vcov_block = crossprod( t(IndicatorMatrix_block))
    KCHR =var_prop[2]*KCHR + var_prop[1]*vcov_block
  }
  
  KCHR.list[[CHR]] = KCHR
  GTCHR.list[[CHR]] = GTCHR
  
}

for(i in 1:n){
  if(perm == T & i %% 100 == 0){print(i)}
  res_LOCO = NULL
  for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
    KCHR = get.K_ASV(GT[,snps$chrom != CHR])
    GTCHR = GT[,snps$chrom == CHR]
    if(perm==T){
      permorder = match(size$strain, sample(unique(size$strain), length(unique(size$strain))))
      GTCHR = GTCHR[permorder,] }
    
    if(fixedBlock){
      reschr = emmax_assoc(y=size$conv,G=GTCHR, covar = fixed_effect_modelMatrix, kinship = KCHR, verbose = !perm)#, weights = size$weight_rotated)
    }else{
      reschr = emmax_assoc(y=size$conv,G=GTCHR, covar = NULL, kinship = KCHR, verbose = !perm)#, weights = size$weight_rotated)
    }
    
    reschr = reschr$results
    reschr = cbind(snps[snps$chrom == CHR,], reschr[,-1])
    res_LOCO = rbind(res_LOCO,reschr)
  }
  
  if(perm==T){pperm_min = aggregate(p~chrom,res_LOCO,min,na.rm=T); permutated_pvals[,i] = pperm_min$p}
}
threshold_alpha0.05 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.05)) # 4.454741 
threshold_alpha0.1 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.1)) # 4.063658 
#save(permutated_pvals, file= "analysis/permutated_pvals_sizeConv.Rdata")

ggplot()+
  #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, size=1, color = 'blue')+
  geom_point(data = res_LOCO, aes(cm, -log10(p)),size=1)+
  facet_grid(.~chrom, scales='free_x')+theme_classic()+
  geom_hline(yintercept = threshold_alpha0.05)


res_LOCO_conv = res_LOCO


fixedBlock = F
perm = F
if(perm==T){n=1000; permutated_pvals = matrix(NA,ncol=n,nrow=6)}else{n=1}

fixed_effect_modelMatrix = model.matrix( ~ block, size)

KCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)
GTCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)

for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
  KCHR = get.K_ASV(GT[,snps$chrom != CHR])
  GTCHR = GT[,snps$chrom == CHR]
  
  if(!fixedBlock){
    modh2 <- mmer(
      div ~ 1,
      random = ~  block + vsr(id, Gu = KCHR),
      data = size, weights = weight_rotated
    )
    
    var_prop = summary(modh2)$varcomp[2:3,1]/sum(summary(modh2)$varcomp[2:3,1])
    IndicatorMatrix_block = model.matrix(~ block-1, size)
    vcov_block = crossprod( t(IndicatorMatrix_block))
    KCHR =var_prop[2]*KCHR + var_prop[1]*vcov_block
  }
  
  KCHR.list[[CHR]] = KCHR
  GTCHR.list[[CHR]] = GTCHR
  
}

for(i in 1:n){
  if(perm == T & i %% 100 == 0){print(i)}
  res_LOCO = NULL
  for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
    KCHR = get.K_ASV(GT[,snps$chrom != CHR])
    GTCHR = GT[,snps$chrom == CHR]
    if(perm==T){
      permorder = match(size$strain, sample(unique(size$strain), length(unique(size$strain))))
      GTCHR = GTCHR[permorder,] }
    
    if(fixedBlock){
      reschr = emmax_assoc(y=size$div,G=GTCHR, covar = fixed_effect_modelMatrix, kinship = KCHR, verbose = !perm)#, weights = size$weight_rotated)
    }else{
      reschr = emmax_assoc(y=size$div,G=GTCHR, covar = NULL, kinship = KCHR, verbose = !perm)#, weights = size$weight_rotated)
    }
    
    reschr = reschr$results
    reschr = cbind(snps[snps$chrom == CHR,], reschr[,-1])
    res_LOCO = rbind(res_LOCO,reschr)
  }
  
  if(perm==T){pperm_min = aggregate(p~chrom,res_LOCO,min,na.rm=T); permutated_pvals[,i] = pperm_min$p}
}
threshold_alpha0.05 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.05)) # 3.914705 
threshold_alpha0.1 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.1)) # 3.595249 
#save(permutated_pvals, file= "analysis/permutated_pvals_sizeDiv.Rdata")

ggplot()+
  #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, size=1, color = 'blue')+
  geom_point(data = res_LOCO, aes(cm, -log10(p)),size=1)+
  facet_grid(.~chrom, scales='free_x')+theme_classic()+
  geom_hline(yintercept = threshold_alpha0.05)


res_LOCO_div = res_LOCO



pgwas = ggplot(res_LOCO_div)+
  facet_grid(.~chrom, scales='free_x')+theme_Publication3()+
  coord_cartesian(ylim = c(0,-log10(min( rbind(res_LOCO_div, res_LOCO_conv)$p))+0.6))+
  geom_hline(yintercept = 3.914705 , color =  "#B56516")+
  geom_hline(yintercept = 4.454741 , color = "#F8D1A2")+
  #geom_hline(yintercept = threshold_alpha0.1, color = "#0073B3")+
  #geom_rect(data=data.frame(chrom="I"), aes(xmin = linked[1], xmax=linked[2],ymin=-Inf, ymax = Inf), fill = "yellow", color=NA, alpha=0.5)+
  geom_point(data = res_LOCO_conv, size=0.25,  aes(cm, -log10(p)), color = 'grey')+
  geom_point(data = res_LOCO_div, size=0.25,  aes(cm, -log10(p)), color = 'black')+
  geom_point( data = data.frame(chrom = unique(res_LOCO$chrom), x=0,y=0), aes(x,y),color = NA)+
  #geom_point(data = subset(snps2, -log10(pval) > threshold_alpha0.05), alpha=1, size=1.2,  aes(cm, -log10(pval)), fill = "red", shape = 21, color='black')+
  theme(panel.spacing = unit(0, "lines"))+
  scale_x_continuous(breaks = c(0,15,30, 45))+
  theme(legend.position = "none")+
  xlab("Genetic distance (cM)")

 


ggsave(pgwas, file="figures/fig_gwas_SNP_size_LOCO.png", width=3.6, height=1.15, dpi=1200)
ggsave(pgwas, file="figures/fig_gwas_SNP_size_LOCO.pdf", width=3.6, height=1.15, dpi=1200)




# perm = F
# if(perm==T){n=1000; permutated_pvals = matrix(NA,ncol=n,nrow=6)}else{n=1}
# for(i in 1:n){
#   if(perm == T & i %% 100 == 0){print(i)}
#   res_LOCO = NULL
#   for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
#     KCHR = get.K_ASV(GT[,snps$chrom != CHR])
#     GTCHR = GT[,snps$chrom == CHR]
#     if(perm==T){GTCHR = GTCHR[sample(1:nrow(GTCHR), nrow(GTCHR)),] }
#     reschr = emmax_assoc(y=size2$emmean_conv,G=GTCHR, covar = NULL, kinship = KCHR, verbose = !perm)
#     reschr = reschr$results
#     reschr = cbind(snps[snps$chrom == CHR,], reschr[,-1])
#     res_LOCO = rbind(res_LOCO,reschr)
#   }
#   
#   if(perm==T){pperm_min = aggregate(p~chrom,res_LOCO,min,na.rm=T); permutated_pvals[,i] = pperm_min$p}
# }
# 
# 
# #-log10(apply(permutated_pvals, 1, function(x){quantile(x, prob = 0.05)}))
# 
# #threshold_alpha0.05 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.05))
# #threshold_alpha0.1 = -log10(quantile( apply(permutated_pvals, 2, min) , prob = 0.1))
# 
# # threshold_alpha0.1 = 3.793954
# # threshold_alpha0.05 = 4.106768
# ggplot()+
#   #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, size=1, color = 'blue')+
#   geom_point(data = res_LOCO, aes(cm, -log10(p)),size=1)+
#   facet_grid(.~chrom, scales='free_x')+theme_classic()#+geom_hline(yintercept = threshold_alpha0.05)
# 
# 
# 
# 
# 
