library(readr)
library(data.table)
library(ggplot2)
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/")
source("utils.R")
meta = read_csv("suppl/RILs_sequencing_metadata.csv")

# GWAS 
load(file = "analysis/QTLmapping/GWAS_EMMAX_LOCO_growthrate.Rdata")


# Founder haplotypes
FGT = NULL
for(CHR in unique(res_LOCO$chrom)){
  print(CHR)
  fgt = as.matrix(read_csv(paste0("genotypes/",CHR,"_becei_genotypes_founders.csv.gz")))
  snps = as.data.frame(read_csv(paste0("genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz")))
  fgt = fgt[match(subset(res_LOCO, chrom == CHR)$id,snps$ID),]
  if(CHR == "X"){fgt[,5] = NA}
  FGT = rbind(FGT,fgt)
}


# Growth rate
#growthrates = read.csv( "phenotypes/growthrates.csv", sep = " ")
#growthrates$log_hours_to_starve = log(growthrates$hours_to_starve)


genotypes =  as.matrix(read_csv("genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
snps =  as.data.frame(fread("genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))

#growthrates = subset(growthrates, grepl("A_", strain))
#genotypes = genotypes[,colnames(genotypes) %in% tstrsplit(growthrates$strain,"_")[[2]]]

genotypes = genotypes[,colnames(genotypes) %in% subset(meta, panel == 'alpha')$rilname]
genotypes = genotypes[match(res_LOCO$id,snps$id),]
snps = snps[match(res_LOCO$id,snps$id),]

FGT = FGT[,-c(3,4)] # Remove founder beta

#af = apply(genotypes,1,mean, na.rm = T)
#genotypes = genotypes[af > 0 & af < 1,]
#snps = snps[af > 0 & af < 1,]

# Polarize so 1 is the minor allele
af = apply(genotypes,1,mean, na.rm = T)
FGT = abs(FGT - ifelse(af>0.5, 1,0)) 
genotypes = abs(genotypes - ifelse(af>0.5, 1,0)) 

res_LOCO$beta = res_LOCO$beta*-1 # The change of sign is because the GWAS was done on log(hours_to_starve), which is inversely proportional to growthrate
res_LOCO$beta = res_LOCO$beta*ifelse(af>0.5, -1,1) # polarize
res_LOCO$t = res_LOCO$beta/res_LOCO$SE # t-value

#res_LOCO$founder = apply(FGT, 1, function(x){ paste(colnames(FGT)[which(x==1)],collapse = ';')})
#res_LOCO$nfounder = apply(FGT, 1, function(x){ length(which(x==1))})

res_LOCO$rils_allelefreq = apply( genotypes, 1, mean,na.rm=T)
res_LOCO$founders_allelefreq = apply(FGT, 1, mean,na.rm=T)
res_LOCO$afc = res_LOCO$rils_allelefreq - res_LOCO$founders_allelefreq # allele frequency change
res_LOCO$s = (qlogis(res_LOCO$rils_allelefreq) - qlogis(res_LOCO$founders_allelefreq)) # selection coefficient
isinf = which(is.infinite(res_LOCO$s))




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

res_LOCO = res_LOCO[-isinf,]
FGT = FGT[-isinf,]

#qtl_snps = subset(res_LOCO, !is.na(asso))
#save(qtl_snps,file = "analysis/temp/qtl_snps.Rdata")

p=ggplot(data = res_LOCO, aes(afc, t))+theme_Publication3()+
  geom_hline(yintercept = 0,color = 'grey30')+geom_vline(xintercept = 0, color = 'grey30')+
  geom_point(data = subset(res_LOCO, is.na(asso)), size = 0.8, color = 'grey', shape=1, alpha= 0.5)+
  geom_point(data = subset(res_LOCO, !is.na(asso)), aes(color= as.factor(asso)), size = 1.5)+
  geom_point(data = subset(res_LOCO, leadingsnp == T), aes(fill= as.factor(asso)), size = 1.9, shape = 21, color = "black")+
  geom_smooth(data = res_LOCO, method='lm', linetype = 'dashed', se=F, color = 'black')+
  scale_color_manual(values = c("#D2605E","#488A8B","#DB8539","#3D7EB6"), name = 'QTL for growth rate')+
  scale_fill_manual(values = c("#D2605E","#488A8B",  "#DB8539","#3D7EB6"), name = 'QTL for growth rate')+
  ylab("Minor allele's t-value for growth rate")+
  xlab("Allele frequency change during panel derivation")+
  theme(legend.key.size = unit(0.7,"line"))


#ggsave(p, file="figures/Fig_tvalue~s.png", width=3, height=1.85, dpi=1200)



p2 = ggplot(data = res_LOCO, aes(afc, t))+theme_Publication3()+
  geom_hline(yintercept = 0,color = 'grey30')+geom_vline(xintercept = 0, color = 'grey30')+
  geom_point(data = subset(res_LOCO, is.na(asso)), size = 0.8, color = 'grey', shape=1, alpha= 0.5)+
  geom_point(data = subset(res_LOCO, !is.na(asso)), aes(color= as.factor(asso)), size = 1.5)+
  geom_point(data = subset(res_LOCO, leadingsnp == T), aes(fill= as.factor(asso)), size = 1.9, shape = 21, color = "black")+
  geom_smooth(data = res_LOCO, method='lm', linetype = 'dashed', se=F, color = 'black')+
  scale_color_manual(values = c("#D2605E","#488A8B","#DB8539","#3D7EB6"), name = 'QTL for growth rate')+
  scale_fill_manual(values = c("#D2605E","#488A8B",  "#DB8539","#3D7EB6"), name = 'QTL for growth rate')+
  ylab("Minor allele's t-value for growth rate")+
  xlab("Allele frequency change during panel derivation")+
  theme(legend.key.size = unit(0.7,"line"))+
  facet_wrap(~chrom)

ggsave(p2, file="figures/SFig_tvalue~afc_chrom.png", width=5, height=3, dpi=1200)



# Statistics
# Null distribution by permutating founder haplotype frequencies in the RILs
# Yielding a null distribution of allele frequency change

# founder haplotype frequencies
foundersfreq = read_csv("suppl/founder_haplotype_frequencies.csv.gz")
#foundersfreq = foundersfreq[paste0(foundersfreq$chrom, foundersfreq$pos) %in% paste0(res_LOCO$chrom, res_LOCO$pos),]
foundersfreq = subset(foundersfreq, cross == 'alpha')

x=split(foundersfreq, foundersfreq$chrom)[[1]]
foundersfreq = lapply(split(foundersfreq, foundersfreq$chrom), function(x){
  
  x2 <- dcast(
    as.data.table(x),
    pos ~ founder,
    value.var = "freq",
    fun.aggregate = mean
  )
  
  x2 = x2[match(subset(res_LOCO, chrom == x$chrom[1])$pos, x2$pos),]
  cbind(chrom = x$chrom[1] , x2)
})






qtlchrom =  (res_LOCO$chrom %in% c("I","II","X"))
robs = cor(res_LOCO$t, res_LOCO$afc, use = "complete.obs" )
robs_withoutqtl = cor(res_LOCO$t[!qtlchrom], res_LOCO$afc[!qtlchrom], use = "complete.obs" )

perm_r = do.call(rbind,lapply(1:10000, function(n){
  if(n %% 100 == 0){print(n)}
  permutated_afc = do.call(c,lapply(foundersfreq, function(X){
   
    X = as.data.frame(X)
    CHR = X$chrom[1]
    X = X[,3:ncol(X)]
    cnames = colnames(X)
    colnames(X) = cnames
    
    FGTperm = FGT[res_LOCO$chrom == CHR,sample(cnames)]
    
    colnames(FGTperm) = cnames
    
    fouders_af_perm = apply(FGTperm, 1, mean, na.rm=T)
   
    
    #plot(out-fouders_allelefreq,
    #     res_LOCO$rils_allelefreq[res_LOCO$chrom == CHR]-res_LOCO$founders_allelefreq[res_LOCO$chrom == CHR])
    
    rilaf_perm = apply(X * FGTperm, 1, sum)
    #s_perm = qlogis(rilaf_perm)-qlogis(fouders_af_perm)
    afc_perm = rilaf_perm - fouders_af_perm
    
    return(afc_perm)
  }))
 
  rperm = cor(res_LOCO$t, permutated_afc, use = "complete.obs" )
  rperm_withoutqtl = cor(res_LOCO$t[!qtlchrom], permutated_afc[!qtlchrom], use = "complete.obs" )
  
  data.frame(rperm,rperm_withoutqtl)
 
}))

#save(perm_r, file = "analysis/temp/null_QTL_afc_correlation.Rdata")
load("analysis/temp/null_QTL_afc_correlation.Rdata")

quantile(perm_r$rperm, prob = 0.95)
quantile(perm_r$rperm_withoutqtl, prob = 0.95)

robs_pvalue = mean(perm_r$rperm >= robs) # prop. of permutation where correlation is higher 
robs_withoutqtl_pvalue = mean(perm_r$rperm_withoutqtl >= robs_withoutqtl) # prop. of permutation where correlation is higher 



library("ggridges")
library(patchwork)

p_perm  = ggplot(perm_r, aes(x = rperm, y = 1, fill = stat(quantile))) +
  theme_Publication3()+
  stat_density_ridges(
    geom = "density_ridges_gradient",
    #calc_ecdf = TRUE,
    quantile_lines = TRUE,
    quantiles = c(0.95), color = "black", linewidth = 0.02)+
  scale_fill_manual(values = c('grey90',"#F6E093"))+
  xlab("Pearson Correleation (r)")+
  ylab("Density")+
  theme(legend.position = "none")+
  coord_cartesian(expand = 0)+
  geom_vline(xintercept = robs, size = 0.35, linetype = "dashed")+
  theme(axis.title = element_text(size = 5))+
  theme(axis.text = element_text(size = 4))+
  theme(plot.margin = unit(c(1,1,2,1), "pt"))+
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )



p = ggplot(data = res_LOCO, aes(afc, t))+theme_Publication3()+
  geom_hline(yintercept = 0,color = 'grey30')+geom_vline(xintercept = 0, color = 'grey30')+
  geom_point(data = subset(res_LOCO, is.na(asso)), size = 0.8, color = 'grey', shape=1, alpha= 0.5)+
  geom_point(data = subset(res_LOCO, !is.na(asso)), aes(color= as.factor(asso)), size = 1.5)+
  geom_point(data = subset(res_LOCO, leadingsnp == T), aes(fill= as.factor(asso)), size = 1.9, shape = 21, color = "black")+
  geom_smooth(data = res_LOCO, method='lm', linetype = 'dashed', se=F, color = 'black')+
  scale_color_manual(values = c("#D2605E","#488A8B","#DB8539","#3D7EB6"), name = 'QTL for growth rate')+
  scale_fill_manual(values = c("#D2605E","#488A8B",  "#DB8539","#3D7EB6"), name = 'QTL for growth rate')+
  ylab("Minor allele's t-value for growth rate")+
  xlab("Allele frequency change during panel derivation")+
  theme(legend.key.size = unit(0.7,"line"))+
  xlim(-0.152,0.19)+ylim(-11,5.463)+
  inset_element(
    p_perm,
    left = 0.64,
    bottom = 0.01,
    right = 1.06,
    top = 0.45
  )
  

ggsave(p, file="figures/Fig_tvalue~afc.png", width=3.5, height=2.2, dpi=1200)

robs


# ggplot(perm_r, aes(x = rperm_withoutqtl, y = 1, fill = stat(quantile))) +
#   theme_Publication3()+
#   stat_density_ridges(
#     geom = "density_ridges_gradient",
#     #calc_ecdf = TRUE,
#     quantile_lines = TRUE,
#     quantiles = c(0.95), color = NA)+
#   scale_fill_manual(values = c('grey',"#E69900"))+
#   xlab("Pearson Correleation (r)")+
#   ylab("Density")+
#   theme(legend.position = "none")+
#   coord_cartesian(expand = 0)+
#   geom_vline(xintercept = robs_withoutqtl, size = 1, linetype = "dashed")





