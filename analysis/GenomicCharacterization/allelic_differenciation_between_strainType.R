library(data.table)
library(parallel)
setwd("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/")
source("utils.R")
geno = fread("genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz")
geno = as.matrix(geno)

####################
### Permutations ###

meta = fread("suppl/RILs_sequencing_metadata.csv")
meta = meta[order(meta$panel),]

nperm = 500

for(n in 1:nperm){
  
  if(n %% 50 == 0){print(n)}
  
  # permute strain type
  meta$strain_type = unlist(lapply(split(meta$strain_type, meta$panel), sample))
  
  # extract allelic count per strain type: [SNV, c(alt, total)]
  count_continuous = t(apply(geno[,colnames(geno) %in% meta$rilname[meta$strain_type == "ContinuousSibmating"]], 1, function(x) c(sum(x, na.rm=T), sum(!is.na(x)))))
  count_restored = t(apply(geno[,colnames(geno) %in% meta$rilname[meta$strain_type == "RestoredFromBackup"]], 1, function(x) c(sum(x, na.rm=T), sum(!is.na(x)))))
  
  # compute p-value for difference between allelic count at each SNP
  pvals = unlist(mclapply(seq_len(nrow(count_continuous)), mc.cores = 3 , function(i) {
    prop.test(
      x = c(count_continuous[i, 1]+1, count_restored[i, 1]+1),
      n = c(count_continuous[i,2]+1, count_restored[i, 2]+1),
      correct = FALSE
    )$p.value
  }))
  
  # extract minimum p-value
  minp = min(pvals, na.rm = T)
  outfile = "analysis/temp/permutated_pval_Continuous_vs_restoredFromBackup.txt"
  cat(minp, file = outfile, append = TRUE, sep = "\n")
}

################
### Observed ###
snps = fread("genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz")
meta = fread("suppl/RILs_sequencing_metadata.csv")
# extract allelic count per strain type: [SNV, c(alt, total)]
count_continuous = t(apply(geno[,colnames(geno) %in% meta$rilname[meta$strain_type == "ContinuousSibmating"]], 1, function(x) c(sum(x, na.rm=T), sum(!is.na(x)))))
count_restored = t(apply(geno[,colnames(geno) %in% meta$rilname[meta$strain_type == "RestoredFromBackup"]], 1, function(x) c(sum(x, na.rm=T), sum(!is.na(x)))))

# compute p-value for difference between allelic count at each SNP
snps$p = unlist(mclapply(seq_len(nrow(count_continuous)), mc.cores = 3 , function(i) {
  prop.test(
    x = c(count_continuous[i, 1]+1, count_restored[i, 1]+1),
    n = c(count_continuous[i,2]+1, count_restored[i, 2]+1),
    correct = FALSE
  )$p.value
}))

threshold0.05 = -log10(quantile(unlist(read.table("analysis/temp/permutated_pval_Continuous_vs_restoredFromBackup.txt")),prob = 0.05))


ggplot(snps, aes(cm, -log10(p)))+
  theme_classic()+
  geom_point()+
  geom_hline(yintercept = threshold0.05, color = 'orange')+
  facet_wrap(~chrom, nrow=1)

p=ggplot(snps)+
  facet_grid(.~chrom, scales='free_x')+theme_Publication3()+
  coord_cartesian(ylim = c(0,-log10(min(snps$p))+0.6))+
  geom_hline(yintercept = threshold0.05, color = "#CC6600")+
  geom_point(alpha=1, size=0.25,  aes(cm, -log10(p)))+
  geom_point( data = data.frame(chrom = unique(snps$chrom), x=0,y=0), aes(x,y),color = NA)+
  theme(panel.spacing = unit(0, "lines"))+
  scale_x_continuous(breaks = c(0,15,30, 45))+
  theme(legend.position = "none")+
  xlab("Genetic distance (cM)")


ggsave(p, file="figures/sfig_allelicDifferenciation_strainType.png", width=3.6, height=1.15, dpi=1200)
ggsave(p, file="figures/sfig_allelicDifferenciation_strainType.pdf", width=3.6, height=1.15, dpi=1200)






