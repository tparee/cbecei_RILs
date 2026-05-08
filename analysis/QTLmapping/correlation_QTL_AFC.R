library(readr)
library(data.table)
library(ggplot2)
source("utils.R")

# 6.641217
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/")
load(file = "analysis/QTLmapping/GWAS_EMMAX_LOCO_growthrate.Rdata")

FGT = NULL
for(CHR in unique(res_LOCO$chrom)){
  print(CHR)
  fgt = as.matrix(read_csv(paste0("genotypes/",CHR,"_becei_genotypes_founders.csv.gz")))
  snps = as.data.frame(read_csv(paste0("genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz")))
  fgt = fgt[match(subset(res_LOCO, chrom == CHR)$id,snps$ID),]
  FGT = rbind(FGT,fgt)
}


growthrates = read.csv( "phenotypes/growthrates.csv", sep = " ")
growthrates$log_hours_to_starve = log(growthrates$hours_to_starve)

# snps =  as.data.frame(fread("genotypes/beceiPanels_snps_pruned0.9999.csv"))
# cmcor = NULL
# for(CHR in unique(res_LOCO$chrom)){
#   print(CHR)
#   snpsunpruned = as.data.frame(read_csv(paste0("genotypes/chr",CHR,"_beceiPanels_snps_founders&Rils_haplotypeResolved_V1.csv.gz")))
#   snpsunpruned = snpsunpruned[match(subset(snps, chrom == CHR)$id,snpsunpruned$ID),]
#   cmcor = c(cmcor,snpsunpruned$cM)
# }


genotypes =  as.matrix(read_csv("genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
snps =  as.data.frame(fread("genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))
genotypes = (genotypes-0.5)*2
growthrates = subset(growthrates, grepl("A_", strain))
genotypes = genotypes[,colnames(genotypes) %in% tstrsplit(growthrates$strain,"_")[[2]]]
af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})
genotypes = genotypes[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]

af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})


FGT = FGT[,-c(3,4)]
FGT = abs(FGT - ifelse(af>0.5, 1,0))
genotypes = genotypes * ifelse(af>0.5, -1,1)
res_LOCO$t = res_LOCO$beta/res_LOCO$SE
res_LOCO$t = res_LOCO$t*ifelse(af>0.5, -1,1)
res_LOCO$b = res_LOCO$b*ifelse(af>0.5, -1,1)

res_LOCO$founder = apply(FGT, 1, function(x){ paste(colnames(FGT)[which(x==1)],collapse = ';')})
res_LOCO$nfounder = apply(FGT, 1, function(x){ length(which(x==1))})

#ggplot()+
#  #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, growthrates=1, color = 'blue')+
#  geom_point(data = res_LOCO, aes(cm, -log10(pval), color = founder), size=1)+
#  facet_grid(.~chrom, scales='free_x')+theme_classic()


res_LOCO$allelefreq = apply( genotypes+1, 1, mean,na.rm=T)/2
res_LOCO$afc = res_LOCO$allelefreq - res_LOCO$nfounder*0.25
res_LOCO$s = (qlogis(res_LOCO$allelefreq) - qlogis(res_LOCO$nfounder*0.25))/5


isinf = which(is.infinite(res_LOCO$s))
cor.test(res_LOCO$t[-isinf], res_LOCO$s[-isinf], use = "complete.obs" )

ggplot(res_LOCO[-isinf,], aes(s, t*-1))+geom_point(size = 0.5, shape=1)


permutated_pvals = unlist(read.table("analysis/temp/permutated_pval_growthrate.txt"))
threshold_alpha0.05 = quantile( permutated_pvals , prob = 0.05) # 7.099286 


qtl1 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'I' & res_LOCO$cm < 40 & res_LOCO$cm > 20,]
qtl1 = qtl1$id[which.min(qtl1$pval)]

qtl2 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'I' & res_LOCO$cm > 42,]
qtl2 = qtl2$id[which.min(qtl2$pval)]

qtl3 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'II' & res_LOCO$cm > 45,]
qtl3 = qtl3$id[which.min(qtl3$pval)]

qtl4 = res_LOCO[res_LOCO$pval < threshold_alpha0.05  & res_LOCO$chrom == 'X' & res_LOCO$cm > 20,]
qtl4 = qtl4$id[which.min(qtl4$pval)]

res_LOCO$leadingsnp = ifelse(res_LOCO$id %in% c(qtl1, qtl2, qtl3, qtl4), T,F)

res_LOCO$asso = NA
for(thisqtl in c(qtl1, qtl2, qtl3, qtl4)){
  print(thisqtl)
  CHR = subset(res_LOCO, id == thisqtl)$chrom
  qtlgt = genotypes[snps$id == thisqtl,]
  asso = unlist(lapply(snps$id[snps$chrom == CHR], function(thislocus){
    r2 = cor(genotypes[snps$id == thislocus,],qtlgt)^2
    if(r2 > 0.8){return(thislocus)}else{return(NULL)}
  }))
  
  res_LOCO$asso[res_LOCO$id %in% asso] = thisqtl
}

subset(res_LOCO, t < -5)

ggplot(data = res_LOCO, aes(s, t*-1))+theme_Publication3()+
  geom_hline(yintercept = 0,color = 'grey30')+geom_vline(xintercept = 0, color = 'grey30')+
  geom_point(data = subset(res_LOCO, is.na(asso)), size = 0.8, color = 'grey', shape=1, alpha= 0.5)+
  geom_point(data = subset(res_LOCO, !is.na(asso)), aes(color= as.factor(asso)), size = 1.5)+
  geom_point(data = subset(res_LOCO, leadingsnp == T), aes(fill= as.factor(asso)), size = 3, shape = 21, color = "black")+
  geom_smooth(data = res_LOCO, method='lm', linetype = 'dashed', se=F, color = 'black')+
  scale_color_manual(values = c("#D2605E","#488A8B", "#3D7EB6","#DB8539"), name = 'QTL for growth rate')+
  scale_fill_manual(values = c("#D2605E","#488A8B", "#3D7EB6", "#DB8539"), name = 'QTL for growth rate')+
  ylab("Minor allele's t-value for growth rate")+
  xlab("Selection coeficient during panel derivation")+
  theme(legend.key.size = unit(0.7,"line"))


ggsave(p, file="figures/fig_tvalue~s.png", width=3, height=1.85, dpi=1200)
