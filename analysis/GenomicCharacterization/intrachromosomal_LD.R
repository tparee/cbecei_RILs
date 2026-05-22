library(data.table)
library(readr)
setwd("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/")

THISPANEL = "beta"
print(THISPANEL)
meta = fread("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/suppl/RILs_sequencing_metadata.csv")
geno = as.matrix(fread("beceiPanels_geno_RILs_pruned0.8.csv.gz"))
snps = as.data.frame(fread("beceiPanels_variantsInfo_pruned0.8.csv.gz"))
geno = geno[,colnames(geno) %in% subset(meta, panel == THISPANEL)$rilname]
af = apply(geno,1,mean)
geno = geno[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]
geno = split(as.data.frame(geno), snps$chrom)

chromcomb = expand.grid(1:length(geno),1:length(geno))
chromcomb = chromcomb[chromcomb[,1] < chromcomb[,2],]


perm = F
if(perm){nitt = 1000}else{nitt = 1}


for(n in 1:nitt){
  
  if(n %% 100 == 0){print(n)}
  
  if(perm){
    GT = lapply(geno, function(X){
      X = X[,sample(1:ncol(X))]
      return(X)
    })
  }else{
    GT = geno
  }
  
  
  interchrom_r2 = do.call(rbind, lapply(1:nrow(chromcomb), function(i){
    chrom1 = chromcomb[i,1]; chrom2 = chromcomb[i,2]
    X = GT[[chrom1]]
    Y = GT[[chrom2]]
    cor_mat <- cor(t(X),t(Y))^2
    nobs = nrow(cor_mat)*ncol(cor_mat)
    meanr2 = mean(cor_mat)
    data.frame(chrom1, chrom2, meanr2, nobs, n)
  }))
  
  if(perm){
    write_csv(interchrom_r2,
              paste0("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/analysis/temp/interchrom_r2_", THISPANEL,"_null.csv"), append = TRUE)
  }else{
    write_csv(interchrom_r2, paste0("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/analysis/temp/interchrom_r2_", THISPANEL,".csv"))
  }
  
  
  
  
}


DATA = list(alpha = NULL, beta = NULL)
for(THISPANEL in c("alpha", "beta")){
  
  obsdata = fread(paste0("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/analysis/temp/interchrom_r2_", THISPANEL,".csv"))
  permdata = fread(paste0("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/analysis/temp/interchrom_r2_", THISPANEL,"_null.csv"))
  colnames(permdata) = c("chrom1", "chrom2", "meanr2", "nobs", "n")
  
  obsdata$panel = THISPANEL; permdata$panel = THISPANEL
  
  pvalue = do.call(rbind, lapply(1:nrow(chromcomb), function(i){
    null = subset(permdata, chrom1 == chromcomb[i,1] & chrom2 == chromcomb[i,2])
    
    
    obs = subset(obsdata, chrom1 == chromcomb[i,1] & chrom2 == chromcomb[i,2])$meanr2
    
    
    data.frame(chrom1 = chromcomb[i,1], chrom2 = chromcomb[i,2], p=  mean(null$meanr2 > obs), panel = THISPANEL)
  }))
  
  
  if(THISPANEL == "beta"){
    colnames(pvalue)[1:2] = colnames(pvalue)[2:1]
    colnames(permdata)[1:2] = colnames(permdata)[2:1]
    colnames(obsdata)[1:2] = colnames(obsdata)[2:1]
  }
  
  DATA[[THISPANEL]] = list(permdata = permdata, pvalue = pvalue, obsdata = obsdata)
}

save(DATA, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/analysis/GenomicCharacterization/iterchrom_LD.Rdata")
pvalue_threhold = 0.05/nrow(chromcomb) # bonferoni correction



obsdata = do.call(rbind, lapply(DATA, function(x){ as.data.frame(x$obsdata)}))
permdata = do.call(rbind, lapply(DATA, function(x){ as.data.frame(x$permdata)}))
pvalue = do.call(rbind, lapply(DATA, function(x){ as.data.frame(x$pvalue)}))

library("ggridges")
source("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/utils.R")

p=ggplot(permdata, aes(x = meanr2, y = 1, fill = panel))+
  theme_Publication3()+
  stat_density_ridges(
    geom = "density_ridges_gradient",
    #calc_ecdf = TRUE
    color = NA)+
  scale_fill_manual(values = c('#396224','#DAB855'))+
  geom_vline(data=obsdata, aes(xintercept = meanr2))+
  geom_text(data=pvalue, aes(x=0.011, y = 1100, label = paste0("p=", round(p, digits = 3))), size = 5/.pt)+
  geom_vline(xintercept =  -Inf, color = 'darkgrey')+
  geom_hline(yintercept =  -Inf, color = 'darkgrey')+
  theme(axis.text.y = element_blank())+
  facet_grid(factor(chrom1, labels = c("I","II","III","IV", "V","X"))~factor(chrom2,  labels = c("I","II","III","IV", "V","X")))+
  coord_cartesian(expand = 0)+
  ylab("")+xlab(expression(r^2))



ggsave(p, file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/figures/SFig_interchrom_linkage.png", width=7.2, height=6, dpi=1200)



