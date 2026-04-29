CHR = "I"

setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/")
source("utils.R")

meta = read_csv("suppl/RILs_sequencing_metadata.csv")

################################################################################
###################### Compute LD ##############################################

param = list(#list(keepdomain = c("arm", "center"), keeppanel = "α"),
                #list(keepdomain = c("arm"), keeppanel = "α"),
                # list(keepdomain = c("arm", "center"), keeppanel = "β"),
               # list(keepdomain = "arm", keeppanel = "β"),
                 # list(keepdomain = "center", keeppanel = "β"),
                 # list(keepdomain = c("arm", "center"), keeppanel = c("α","β")),
                 # list(keepdomain = "arm", keeppanel = c("α","β")),
                  list(keepdomain = "center", keeppanel = c("α","β")))


gbins = seq(0,5,0.1) # genetic distance bins
gbins = cbind(gbins[-length(gbins)], gbins[-1])
pbins = seq(1,1e6,2e4) # physical distance bins
pbins = cbind(pbins[-length(pbins)], pbins[-1])

colnames(pbins) = c("ppos1","ppos2")
colnames(gbins) = c("gpos1","gpos2")

# Iterate across chromosome and calculate LD as r2 = cor(t(geno))^2
# To avoid dealing with huge correlation matrix, caluclate LD per group of 4000 variants randomly sampled genome-wide
# Not every pairs of variants will be tested by it is a good sampling (each variant tested with 3999 other variants across the whole chromosome)
# Keep only the average distance per bin
for(thisparam in param){
  for(CHR in c("I","II",'III',"IV",'V',"X")){
    print(CHR)
    geno = fread(paste0(file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/",CHR,"_becei_genotypes_RILs.csv.gz"))
    geno = as.matrix(geno)
    snps = fread(paste0(file="~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz"))
    
    keepi = snps$domain %in% thisparam$keepdomain
    keepj = colnames(geno) %in% subset(meta, panel %in% thisparam$keeppanel)$rilname
    snps = snps[keepi,]
    geno = geno[keepi,keepj]
    
    fixed = apply(geno, 1, mean, na.rm=T) %in% c(0,1) 
    snps = snps[!fixed,]; geno = geno[!fixed,]
    
    windows = get.win(nrow(snps), winsize = 4000, minsize = 2000)
    set.seed(123)
    sampledorder = sample(1:nrow(snps))
    windows = lapply(windows, function(win){sort(sampledorder[win])})
    
    binnedr2 = lapply(windows, function(win){
      #win = windows[[1]]
      gt = as.matrix(geno[win,])
      r2 = cor(t(gt))^2
      pd = abs(outer(snps$POS[win], snps$POS[win], "-"))
      gd = abs(outer(snps$cM[win], snps$cM[win], "-"))
      gd[lower.tri(gd,diag = T)] = -1
      pd[lower.tri(pd,diag = T)] = -1
      
      r2pbin = sapply(1:nrow(pbins), function(i){
        out = r2[pd >= pbins[i,1]  & pd < pbins[i,2]]
        data.frame(sumr2 = sum(out), nobs = length(out))
      })
      
      r2gbin = sapply(1:nrow(gbins), function(i){
        out = r2[gd >= gbins[i,1]  & gd < gbins[i,2]]
        data.frame(sumr2 = sum(out), nobs = length(out))
      })
      
      r2pbin = cbind(pbins, t(r2pbin))
      r2gbin = cbind(gbins, t(r2gbin))
      
      return(list(r2pbin, r2gbin))
      
    })
    
    R2pbin = Reduce("+", lapply(binnedr2, function(x){ matrix( unlist(x[[1]][,3:4]), ncol=2) }))
    R2gbin = Reduce("+", lapply(binnedr2, function(x){ matrix( unlist(x[[2]][,3:4]), ncol=2) }))
    R2pbin[,1] = R2pbin[,1]/R2pbin[,2]
    R2gbin[,1] = R2gbin[,1]/R2gbin[,2]
    
    colnames(R2gbin) = c("mean.r2", "nobs")
    colnames(R2pbin) = c("mean.r2", "nobs")
    
    R2pbin = cbind(pbins, R2pbin)
    R2gbin = cbind(gbins, R2gbin)
    
    R2gbin = as.data.frame(R2gbin)
    R2pbin = as.data.frame(R2pbin)
    
    save(R2gbin, file = paste0("analysis/temp/",CHR, "_","LD_gbins", "_panel", paste(thisparam$keeppanel, collapse='&'), "_domain", paste(thisparam$keepdomain, collapse = "&"),".Rdata"))
    save(R2pbin, file = paste0("analysis/temp/",CHR, "_","LD_pbins", "_panel", paste(thisparam$keeppanel, collapse='&'), "_domain", paste(thisparam$keepdomain, collapse = "&"),".Rdata"))
  }
}


ldfiles = list.files("analysis/temp/")
ldfiles = ldfiles[grepl("_LD_",ldfiles) & grepl("pbins",ldfiles)]
thisfile = ldfiles[1]
linkage_physicalDistance = NULL
for(thisfile in ldfiles){
  load(paste0("analysis/temp/",thisfile))
  R2pbin$chrom = tstrsplit(thisfile,"_")[[1]]
  R2pbin$panel = sub(".*panel(.*?)_domain.*", "\\1", thisfile)
  R2pbin$domain = sub(".*domain(.*?).Rdata.*", "\\1", thisfile)
  linkage_physicalDistance = rbind(linkage_physicalDistance, R2pbin)
}

ldfiles = list.files("analysis/temp/")
ldfiles = ldfiles[grepl("_LD_",ldfiles) & grepl("gbins",ldfiles)]
thisfile = ldfiles[1]
linkage_geneticDistance = NULL
for(thisfile in ldfiles){
  load(paste0("analysis/temp/",thisfile))
  R2gbin$chrom = tstrsplit(thisfile,"_")[[1]]
  R2gbin$panel = sub(".*panel(.*?)_domain.*", "\\1", thisfile)
  R2gbin$domain = sub(".*domain(.*?).Rdata.*", "\\1", thisfile)
  linkage_geneticDistance = rbind(linkage_geneticDistance, R2gbin)
}


#write_csv(linkage_physicalDistance, file = "analysis/GenomicCharacterization/linkage_by_physicalDistance.csv")
#write_csv(linkage_geneticDistance, file = "analysis/GenomicCharacterization/linkage_by_geneticDistance.csv")

################################################################################
###################### PLOT LD per chrom #######################################

source("utils.R")
linkage_physicalDistance = read_csv(file = "analysis/GenomicCharacterization/linkage_by_physicalDistance.csv")
linkage_geneticDistance = read_csv(file = "analysis/GenomicCharacterization/linkage_by_geneticDistance.csv")
linkage_physicalDistance$panel = factor(linkage_physicalDistance$panel, levels = c("α" , "β" ,"α&β" ))
linkage_geneticDistance$panel = factor(linkage_geneticDistance$panel, levels = c("α" , "β" ,"α&β" ))

pg=ggplot(subset(linkage_geneticDistance, domain == "arm&center"), aes((gpos1+gpos2)/2, mean.r2,
                                             color = panel))+
  geom_smooth(size=0.8, span = 0.2, se = F)+
  #geom_line(size=1)+
  facet_wrap(.~chrom, scales = 'free', nrow=1)+
  theme_Publication3()+
  coord_cartesian(ylim=c(0,1),xlim=c(0,5))+
  scale_color_manual(values=c('#396224','#DAB855','#7CB6CE'))+
  xlab("Genetic position (cM)")+
  ylab(expression(r^2))+
  theme(legend.key.height = unit(0.3, "cm"))

pp=ggplot(subset(linkage_physicalDistance, domain == "arm&center"), aes((ppos1+ppos2)/2e3, mean.r2,
                                                                       color = panel))+
  geom_smooth(size=0.8, span = 0.2, se=F)+
  #geom_line(size=1)+
  facet_wrap(.~chrom, scales = 'free', nrow=1)+
  theme_Publication3()+
  coord_cartesian(ylim=c(0,1),xlim=c(0,10e2))+
  scale_color_manual(values=c('#396224','#DAB855','#7CB6CE'))+
  xlab("Physical position (Kb)")+
  ylab(expression(r^2))+
  theme(legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"))+
  guides(color = guide_legend(override.aes = list(color = "white", fill = "white")))+
  scale_x_continuous(breaks = c(0,500,1000))


pLD = grid.arrange(pg, pp, nrow=2)



################################################################################
###################### PLOT PCA ################################################

genotypes =  as.matrix(fread("genotypes/beceiPanels_geno_RILs_pruned0.999.csv.gz"))
#genotypes = genotypes[,grepl("A_",colnames(genotypes))]
afs = apply(genotypes, 1, function(x){mean(x[x==0 | x==1],na.rm=T)})
snps =  as.data.frame(fread("genotypes/beceiPanels_variantsInfo_pruned0.999.csv.gz"))
ix = (afs > 0 & afs < 1)
genotypes = genotypes[ix,]
genotypes[is.na(genotypes)]=0.5

pca = prcomp(t(genotypes), center = TRUE, scale. = TRUE)
pca = pca$x[, 1:2]
pca = as.data.frame(pca)
pca$id = rownames(pca)
pca$panel = meta$panel[match(pca$id, meta$rilname)]
pPCA = ggplot(pca, aes(PC1,PC2, color=panel))+geom_point(size=2/.pt)+
  scale_color_manual(values=c('#396224','#DAB855'))+
  theme_Publication3()+
  theme(legend.position = 'none',legend.direction = "vertical", legend.key.height = unit(0.2, "cm"))+
  theme(plot.margin = unit(c(12,5,5,10), 'pt'))


#pout = grid.arrange(pPCA, pLD, ncol=2, widths = c(1,5))


################################################################################
###################### PLOT AFS ################################################
afs = apply(genotypes,1,mean, na.rm=T)
afs = abs(afs - ifelse(afs > 0.5, 1, 0))
pafs = ggplot(data = data.frame(afs=afs), aes(x=afs))+
  theme_Publication3()+
  geom_density(aes(y = after_stat(count / sum(count))), fill = 'grey', color = NA)+
  #geom_histogram(aes(y = after_stat(count / sum(count))), bins = 50)+
  coord_cartesian(expand=0,xlim=c(0,0.6))+
  ylab("Density")+
  xlab("Minor allele frequency")+
  theme(plot.margin = unit(c(12,5,5,5), 'pt'))#+
  #scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = as.character(c(0, 0.1, 0.2, 0.3, 0.4, 0.5)))


pleft = grid.arrange(pafs, pPCA, ncol=1, heights = c(1,1))

pout = grid.arrange(pleft, pLD, ncol=2, widths = c(1,4))

ggsave(pout, file="figures/Fig_PCA&Linkage.png", width=7.2, height=2.22, dpi=1200)
ggsave(pout, file="figures/Fig_PCA&Linkage.pdf", width=7.2, height=2.22, dpi=1200)




LDdom = subset(linkage_physicalDistance, domain %in% c("arm", "center"))
LDdom = do.call(rbind, lapply(split(LDdom, paste0(LDdom$ppos1,LDdom$panel, LDdom$domain )), function(x){
  out = x[1,-which(colnames(x)=='chrom')]
  out$mean.r2[is.nan(out$mean.r2)] = 0
  out$mean.r2 = sum(x$nobs*x$mean.r2)/sum(x$nobs)
  out
}))

LDdom$panel = factor(LDdom$panel, levels = c("α" , "β" ,"α&β" ))
pdom=ggplot(LDdom, aes((ppos1+ppos2)/2000,mean.r2, linetype = domain, color = panel))+
  theme_Publication3()+
  geom_line(size=0.8)+facet_wrap(~panel, scales = 'free_y')+ylim(0,1)+
  scale_color_manual(values=c('#396224','#DAB855','#7CB6CE'), name = "Panel" )+
  xlab("Physical position (kb)")+
  ylab(expression(r^2))+
  theme(legend.key.height = unit(0.3, "cm"))

