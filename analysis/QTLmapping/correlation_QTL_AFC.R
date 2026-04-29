library(readr)
library(data.table)
library(ggplot2)
source("utils.R")

load(file = "analysis/QTLmapping/GWAS_EMMAX_LOCO_growthrate.Rdata")

FGT = NULL
for(CHR in unique(res_LOCO$chrom)){
  print(CHR)
  fgt = as.matrix(read_csv(paste0("genotypes/chr",CHR,"_beceiPanels_geno_founders_haplotypeResolved_V1.csv.gz")))
  snps = as.data.frame(read_csv(paste0("genotypes/chr",CHR,"_beceiPanels_snps_founders&Rils_haplotypeResolved_V1.csv.gz")))
  fgt = fgt[match(subset(res_LOCO, chrom == CHR)$id,snps$ID),]
  FGT = rbind(FGT,fgt)
}



setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/becei_rils/")
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


genotypes =  as.matrix(read_csv("genotypes/beceiPanels_geno_RILs_pruned0.9999.csv"))
snps =  as.data.frame(fread("genotypes/beceiPanels_snps_pruned0.9999.csv"))
genotypes = (genotypes-0.5)*2
genotypes[genotypes == 0]=1



growthrates = subset(growthrates, grepl("A_", strain))

colnames(genotypes) = tstrsplit(colnames(genotypes), "-")[[1]]
growthrates = subset(growthrates, strain %in% colnames(genotypes))
genotypes = genotypes[,colnames(genotypes) %in% growthrates$strain]

af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})
genotypes = genotypes[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]

af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})


FGT = FGT[,-c(3,4)]
FGT = abs(FGT - ifelse(af>0.5, 1,0))
genotypes = genotypes * ifelse(af>0.5, -1,1)
res_LOCO$t = res_LOCO$t*ifelse(af>0.5, -1,1)
res_LOCO$b = res_LOCO$b*ifelse(af>0.5, -1,1)

res_LOCO$founder = apply(FGT, 1, function(x){ paste(colnames(FGT)[which(x==1)],collapse = ';')})
res_LOCO$nfounder = apply(FGT, 1, function(x){ length(which(x==1))})

ggplot()+
  #geom_point(data = res, aes(cm, -log10(p)), alpha=0.3, growthrates=1, color = 'blue')+
  geom_point(data = res_LOCO, aes(cm, -log10(p), color = founder), size=1)+
  facet_grid(.~chrom, scales='free_x')+theme_classic()




res_LOCO$allelefreq = apply( genotypes+1, 1, mean,na.rm=T)/2
res_LOCO$afc = res_LOCO$allelefreq - res_LOCO$nfounder*0.25
res_LOCO$s = (qlogis(res_LOCO$allelefreq) - qlogis(res_LOCO$nfounder*0.25))/5



isinf = which(is.infinite(res_LOCO$s))
cor.test(res_LOCO$t[-isinf], res_LOCO$s[-isinf], use = "complete.obs" )

ggplot(res_LOCO, aes(s, b*-1))+geom_point(size = 0.5, shape=1)

res_LOCO$qtlgroup = NA
res_LOCO$qtlgroup[res_LOCO$p < 10^-6.450801 & res_LOCO$chrom == 'I' & res_LOCO$cm < 17] = "I:0.9cM"
res_LOCO$qtlgroup[res_LOCO$p < 10^-6.450801 & res_LOCO$chrom == 'I' & res_LOCO$cm > 19] = "I:30.1cM"
res_LOCO$qtlgroup[res_LOCO$p < 10^-6.450801 & res_LOCO$chrom == 'X'] = "X:23.6cM"

res_LOCO[res_LOCO$p < 10^-6.450801 & res_LOCO$chrom == 'X',] 
ggplot(data = res_LOCO, aes(s, t*-1))+theme_classic()+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(data = subset(res_LOCO,p > 10^-6.450801), size = 1, color = 'grey', shape=1)+
  geom_point(data = subset(res_LOCO,p < 10^-6.450801), aes(color=founder), size = 1.5)+
  geom_smooth(data = res_LOCO, method='lm', linetype = 'dashed', se=F, color = 'black')+
  scale_color_manual(values = c("#455180", "#4D8BA6", "#016250", "#386124","#A4C595"), breaks = c("FM.g2","FM.g1;FM.g2", "FA.g2;FM.g2", "FA.g2","FA.g1"))

p=ggplot(data = res_LOCO, aes(s, t*-1))+theme_Publication3()+
  geom_hline(yintercept = 0,color = 'grey30')+geom_vline(xintercept = 0, color = 'grey30')+
  geom_point(data = subset(res_LOCO,p > 10^-6.450801), size = 0.8, color = 'grey', shape=1)+
  geom_point(data = subset(res_LOCO,p < 10^-6.450801), aes(color= as.factor(qtlgroup)), size = 1.5)+
  geom_smooth(data = res_LOCO, method='lm', linetype = 'dashed', se=F, color = 'black')+
  scale_color_manual(values = c("#BC7E61","#6C8992", "#D0C5A1"), name = 'QTL for growth rate')+
  ylab("Minor allele's t-value for growth rate")+
  xlab("Selection coeficient during panel derivation")+
  theme(legend.key.size = unit(0.7,"line"))


ggsave(p, file="figures/fig_tvalue~s.png", width=3, height=1.85, dpi=1200)

res_LOCO[which.min(res_LOCO$p),]$s

qlogis(0.1768293)-qlogis(0.25)

ggplot(data = res_LOCO, aes(cm, -log10(p)))+theme_classic()+facet_wrap(~chrom, nrow=1)+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(data = subset(res_LOCO,p > 10^-6.450801), size = 1, color = 'grey', shape=1)+
  geom_point(data = subset(res_LOCO,p < 10^-6.450801), aes(color=founder), size = 1)+
  scale_color_manual(values = c("#455180", "#4D8BA6", "#016250", "#386124","#A4C595"), breaks = c("FM.g2","FM.g1;FM.g2", "FA.g2;FM.g2", "FA.g2","FA.g1"))


table(res_LOCO[qtlregion,]$founder)
qtlregion = res_LOCO$chrom == 'I' & res_LOCO$cm > 17 & res_LOCO$cm < 45 & -log10(res_LOCO$p) > 7


snpsqtl = snps[qtlregion,]
nsnps = sum(qtlregion)

GT = t(genotypes[qtlregion,])
GT[is.na(GT)]=0
r2 = cor(GT)^2
r2 = reshape2::melt(r2)

r2$xmin = c(min(snpsqtl$cm),snpsqtl$cm[-nsnps])[r2$Var1]
r2$xmax = snpsqtl$cm[r2$Var1]

r2$ymin = c(min(snpsqtl$cm),snpsqtl$cm[-nsnps])[r2$Var2]
r2$ymax = snpsqtl$cm[r2$Var2]

ggplot(r2, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile()+scale_fill_viridis()

library(viridis)
ggplot(r2, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax, fill = value), color = NA)+
  geom_rect()+scale_fill_viridis()


