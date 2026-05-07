library(readxl)
library(readr)
library(reshape2)
library(data.table)
library(ggplot2)
library(sommer)

PATH = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/"
setwd(PATH)

meta = read_csv("suppl/RILs_sequencing_metadata.csv")

source("utils.R")

WinSize_cM = 1


permutation = F
# Import data 

Ghap = as.matrix(read_csv("analysis/heritability/haplotype_relatedness_matrix/Ghap_haplotypeRelatednessMatrix.csv"))
size = as.data.frame(fread("phenotypes/size_summarystat.csv"))
size = subset(size, strain %in% subset(meta, panel == "alpha")$rilname & strain %in% colnames(Ghap))
size$rilname =size$strain
rownames(Ghap) = colnames(Ghap)
Ghap = Ghap[rownames(Ghap) %in% size$strain,rownames(Ghap) %in% size$strain]


if(permutation){
  niter = 1000
}else{
  niter = 1
}


GHAPCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)
HAPCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)

for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
  
  snps = as.data.frame(read_csv(paste0("genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz")))
  fgt = as.matrix(read_csv(paste0("genotypes/",CHR,"_becei_genotypes_founders.csv.gz")))
  fgt = fgt[,c("FA.g1","FA.g2","FM.g1","FM.g2")]
  
  haplo = as.data.frame(read_csv(paste0("haplotypes/",CHR,"_rils_foundingHaplotypesBlocks.csv")))
  haplo = subset(haplo, rilname %in% size$rilname)
  
  # Add x cM bins
  cuthere = seq(min(snps$cM), max(snps$cM),WinSize_cM)
  cuthere[1] = min(cuthere)-0.1
  cuthere[length(cuthere)] = max(cuthere)+0.1
  snps$win = cut(snps$cM, cuthere)
  
  haplomatrix = do.call(rbind, lapply(split(snps, snps$win), function(win){
    
    #win = split(snps, snps$win)[[45]]
    print(win$win[1])
    win = snps$ID %in% win$ID
    x = fgt[win,]
    afs = apply(x,1,mean, na.rm=T)
    x = x[afs < 1 & afs > 0,]
    
    if(length(x) == ncol(fgt)){return(NULL)}
    if(nrow(x) == 0){return(NULL)}
    
    minpos = min(snps$POS[win],na.rm=T)
    maxpos = max(snps$POS[win],na.rm=T)
    midppos = mean(range(snps$POS[win],na.rm=T))
    midgpos = mean(range(snps$cM[win],na.rm=T))
    
    hapgroups =  cutree(hclust(dist(t(x))), h = 0)
    groupnames = names(table(hapgroups))
    
    
    haplowin = haplo[haplo$pos1 <= midppos & haplo$pos2 > midppos,]
    
    haplowin$overlap = (pmin(maxpos, haplowin$pos2) - pmax(minpos,haplowin$pos1))/(maxpos-minpos)
    
    dt = do.call(rbind, lapply(split(haplowin, haplowin$rilname), function(hap){
      
      #hap = split(haplowin, haplowin$rilname)[[7]]
      #hap = subset(haplowin, rilname == "A_QG3119")
      hap = hap[which.max(hap$overlap),]
      foundernames = unlist(strsplit(hap$founder,";"))
      foundergroup = unique(hapgroups[foundernames])
      
      if(length(foundergroup) > 1){
        foundergroup=NA
      }
      
      data.frame(rilname = hap$rilname[1], haplo=foundergroup)
      
      
    }))
    
    dt$haplo[match(size$rilname,dt$rilname)]
    
  }))
  
  
  HAPCHR.list[[CHR]] = haplomatrix
  
  
  otherchrom = c("I",'II', 'III', 'IV', 'V','X') != CHR
  
  weights = c(rep(50,5),25)[otherchrom]
  weights = weights/sum(weights)
  Ghap=NULL
  for(CHR2 in c("I","II","III",'IV',"V",'X')[otherchrom]){
    g = fread(file = paste0("analysis/heritability/haplotype_relatedness_matrix/chr",CHR2, "_Ghap.csv"))
    g = as.matrix(g)
    rownames(g)=colnames(g)
    g=g*weights[ which(c("I","II","III",'IV',"V",'X')[otherchrom]==CHR2) ]
    if(is.null(Ghap)){Ghap=g}else{
      Ghap=Ghap+g
    }
  }
  
  GHAPCHR.list[[CHR]] = Ghap
  
}



for(n in 1:niter){
  
  print(n)
  hapgwas = NULL
  
  for(CHR in c("I","II","III",'IV',"V",'X')){
    
    haplomatrix = HAPCHR.list[[CHR]]
    Ghap = GHAPCHR.list[[CHR]]
    Ghap = Ghap + diag( rep(1e-4, nrow(Ghap)))
    gpos = sapply(regmatches(rownames(haplomatrix), gregexpr("-?\\d+\\.?\\d*", rownames(haplomatrix), perl = TRUE)), function(x){mean(as.numeric(x))})
    
    
    if(permutation){
      permorder = match(size$strain, sample(unique(size$strain), length(unique(size$strain))))
      haplomatrix = haplomatrix[,permorder] }
    
    
    
    
    res = lapply(1:nrow(haplomatrix), function(i){
      
      dt = cbind(size, haplo = haplomatrix[i,])
      dt$haplo =  as.factor(dt$haplo)
      dt$id = dt$rilname
      dt = na.omit(dt)
      dt$weight_rotated = 1/(dt$SE_rotated^2)
      # Using mmes to fitr the LMM
      # Null model
      fit0 <- mmes(
        div ~ 1,
        random = ~ block +  vsm(ism(id), Gu=Ghap),
        data = dt, verbose=F, W = diag(dt$weight_rotated)
      )
      
      # Model with haplo
      fit <- mmes(
        conv ~ haplo,
        random = ~ block +  vsm(ism(id), Gu=Ghap),
        data = dt, verbose = F, W = diag(dt$weight_rotated)
      )
      
      pval = anova.mmes(fit0, fit)[2,"PrChisq"]
      pval = as.numeric(tstrsplit(pval," ")[[1]])
      
      #anova(lm(emmean ~ haplo, dt))
      
      ### Aletrnatively, use  lme4breeding::lmebreed
      
      
      #fit0 = lmebreed(div ~ 1 + (1|id),
      #               relmat = list(id = Ghap ),
      #               verbose = FALSE, data=dt)

      #fit = lmebreed(div ~ haplo + (1|id),
      #         relmat = list(id = Ghap ),
      #         verbose = FALSE, data=dt, weights = weight_rotated)

      #anova(fit0,fit, test='LRT')
      
      data.frame(chrom=CHR,gpos=gpos[i],pval=pval)
      
    })
    
    res = do.call(rbind, res)
    #ggplot(hapgwas, aes(gpos, -log10(pval)))+geom_line()+facet_wrap(~chrom)
    
    
    hapgwas = rbind(hapgwas,res)
    
  }
  
  if(permutation){
    minp = min(hapgwas$pval)
    write(minp, file = 'analysis/QTLmapping/haplotypeGWAS_permutation_results_sizeConv_LOCO.txt', append = TRUE)
    
  }else{
    save(hapgwas, file = "analysis/QTLmapping/haplotypeGWAS_sizeConv_results_LOCO.Rdata")
  }
  
  
}






load("analysis/QTLmapping/haplotypeGWAS_sizeDiv_results_LOCO.Rdata")

hapgwas = do.call(rbind, lapply(split(hapgwas, hapgwas$chrom), function(x){
  
  CHR = x$chrom[1]
  snps = as.data.frame(read_csv(paste0("genotypes/chr",CHR,"_beceiPanels_snps_founders&Rils_haplotypeResolved_V1.csv.gz")))
  
  y = x[c(1,1),]
  
  #y$ppos = range(snps$POS)
  y$gpos = range(snps$cM)
  y$pval = 1
  
  x=rbind(x,y)
  x=x[order(x$gpos),]
  x
  
}))

hapgwas$colorchrom = ifelse(hapgwas$chrom %in% c('II','IV','X'), 'grey30', 'black')


perm = unlist(read.table("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/becei_rils/analysis/QTLmapping/haplotypeGWAS_permutation_results_sizeDiv_LOCO.txt", quote="\"", comment.char=""))

threshold_alpha0.1 = -log10(quantile(perm, 0.1))
threshold_alpha0.05 = -log10(quantile(perm, 0.05)) # 4.20412

#qtllim = range(subset(hapgwas, -log10(pval)>threshold_alpha0.05)$gpos)

p= ggplot(hapgwas)+
  #geom_hline(yintercept = threshold_alpha0.1, color = "#0073B3")+
  geom_hline(yintercept = threshold_alpha0.05, color = "#CC6600")+
  #geom_rect(data=data.frame(chrom="I"), aes(xmin = qtllim[1], xmax=qtllim[2],ymin=-Inf, ymax = Inf), fill = "#CC6600", color=NA, alpha=0.2)+
  geom_line(aes(gpos, -log10(pval)))+
  facet_grid(.~chrom, scales='free_x')+
  #scale_color_manual(values = unique(hapgwas$colorchrom), breaks =  unique(hapgwas$colorchrom))+
  theme_Publication3()+
  theme(panel.spacing = unit(0, "lines"))+
  coord_cartesian(ylim = c(0,-log10(min(hapgwas$pval))+0.02))+
  scale_x_continuous(breaks = c(0,15,30, 45))+
  theme(legend.position = "none")+
  xlab("Genetic distance (cM)")


#hapgwas_div = hapgwas
ggsave(p, file="figures/fig_gwas_haplotype_sizeDiv_LOCO.png", width=3.6, height=1.15, dpi=1200)
ggsave(p, file="figures/fig_gwas_haplotype_sizeDiv_LOCO.pdf", width=3.6, height=1.15, dpi=1200)



#################################################################################
#-------------------------------------------------------------------------------
#------------------------  SEXUALLY CONVERGENT ---------------------------------
#-------------------------------------------------------------------------------
################################################################################

library(readxl)
library(readr)
library(reshape2)
library(data.table)
library(ggplot2)
library(sommer)

PATH = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/becei_rils/"
setwd(PATH)

source("utils.R")

WinSize_cM = 1


permutation = F
# Import data 

Ghap = as.matrix(read_csv("analysis/heritability/haplotypeRelatednessMatrix.csv"))

size = as.data.frame(fread("phenotypes/size_summarystat.csv"))
size = subset(size, grepl("A_", strain) & strain %in% colnames(Ghap))
size$rilname =size$strain
rownames(Ghap) = colnames(Ghap)
Ghap = Ghap[rownames(Ghap) %in% size$strain,rownames(Ghap) %in% size$strain]


if(permutation){
  niter = 100
}else{
  niter = 1
}


GHAPCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)
HAPCHR.list = list(I=NULL, II=NULL, III=NULL, IV = NULL, V=NULL, X=NULL)

for(CHR in c("I",'II', 'III', 'IV', 'V','X')){
  
  snps = as.data.frame(read_csv(paste0("genotypes/chr",CHR,"_beceiPanels_snps_founders&Rils_haplotypeResolved_V1.csv.gz")))
  fgt = as.matrix(read_csv(paste0("genotypes/chr",CHR,"_beceiPanels_geno_founders_haplotypeResolved_V1.csv.gz")))
  fgt = fgt[,c("FA.g1","FA.g2","FM.g1","FM.g2")]
  
  haplo = as.data.frame(read_csv(paste0("haplotypes/chr",CHR,"_RILs_founderHaplotypeBlocks_V1.csv")))
  haplo = subset(haplo, grepl("A_",haplo$rilname))
  
  # Add x cM bins
  cuthere = seq(min(snps$cM), max(snps$cM),WinSize_cM)
  cuthere[1] = min(cuthere)-0.1
  cuthere[length(cuthere)] = max(cuthere)+0.1
  snps$win = cut(snps$cM, cuthere)
  
  haplomatrix = do.call(rbind, lapply(split(snps, snps$win), function(win){
    
    #win = split(snps, snps$win)[[45]]
    print(win$win[1])
    win = snps$ID %in% win$ID
    x = fgt[win,]
    afs = apply(x,1,mean, na.rm=T)
    x = x[afs < 1 & afs > 0,]
    
    if(length(x) == ncol(fgt)){return(NULL)}
    if(nrow(x) == 0){return(NULL)}
    
    minpos = min(snps$POS[win],na.rm=T)
    maxpos = max(snps$POS[win],na.rm=T)
    midppos = mean(range(snps$POS[win],na.rm=T))
    midgpos = mean(range(snps$cM[win],na.rm=T))
    
    hapgroups =  cutree(hclust(dist(t(x))), h = 0)
    groupnames = names(table(hapgroups))
    
    
    haplowin = haplo[haplo$pos1 <= midppos & haplo$pos2 > midppos,]
    
    haplowin$overlap = (pmin(maxpos, haplowin$pos2) - pmax(minpos,haplowin$pos1))/(maxpos-minpos)
    
    dt = do.call(rbind, lapply(split(haplowin, haplowin$rilname), function(hap){
      
      #hap = split(haplowin, haplowin$rilname)[[7]]
      #hap = subset(haplowin, rilname == "A_QG3119")
      hap = hap[which.max(hap$overlap),]
      foundernames = unlist(strsplit(hap$founders,";"))
      foundergroup = unique(hapgroups[foundernames])
      
      if(length(foundergroup) > 1){
        foundergroup=NA
      }
      
      data.frame(rilname = hap$rilname[1], haplo=foundergroup)
      
      
    }))
    
    dt$haplo[match(size$rilname,dt$rilname)]
    
  }))
  
  
  HAPCHR.list[[CHR]] = haplomatrix
  
  
  otherchrom = c("I",'II', 'III', 'IV', 'V','X') != CHR
  
  weights = c(rep(50,5),25)[otherchrom]
  weights = weights/sum(weights)
  Ghap=NULL
  for(CHR2 in c("I","II","III",'IV',"V",'X')[otherchrom]){
    g = fread(file = paste0("analysis/heritability/haplotypeRelatednessMatrix/chr",CHR2, "_Ghap.csv"))
    g = as.matrix(g)
    rownames(g)=colnames(g)
    g=g*weights[ which(c("I","II","III",'IV',"V",'X')[otherchrom]==CHR2) ]
    if(is.null(Ghap)){Ghap=g}else{
      Ghap=Ghap+g
    }
  }
  
  GHAPCHR.list[[CHR]] = Ghap
  
}



for(n in 1:niter){
  
  print(n)
  hapgwas = NULL
  
  for(CHR in c("I","II","III",'IV',"V",'X')){
    
    haplomatrix = HAPCHR.list[[CHR]]
    Ghap = GHAPCHR.list[[CHR]]
    Ghap = Ghap + diag( rep(1e-4, nrow(Ghap)))
    gpos = sapply(regmatches(rownames(haplomatrix), gregexpr("-?\\d+\\.?\\d*", rownames(haplomatrix), perl = TRUE)), function(x){mean(as.numeric(x))})
    
    
    if(permutation){
      permorder = match(size$strain, sample(unique(size$strain), length(unique(size$strain))))
      haplomatrix = haplomatrix[,permorder] }
    
    
    
    
    res = lapply(1:nrow(haplomatrix), function(i){
      
      dt = cbind(size, haplo = haplomatrix[i,])
      dt$haplo =  as.factor(dt$haplo)
      dt$id = dt$rilname
      dt = na.omit(dt)
      dt$weight_rotated = 1/(dt$SE_rotated^2)
      # Using mmes to fitr the LMM
      # Null model
      fit0 <- mmes(
        conv ~ 1,
        random = ~ block +  vsm(ism(id), Gu=Ghap),
        data = dt, verbose=F, W = diag(dt$weight_rotated)
      )
      
      # Model with haplo
      fit <- mmes(
        conv ~ haplo,
        random = ~ block +  vsm(ism(id), Gu=Ghap),
        data = dt, verbose = F, W = diag(dt$weight_rotated)
      )
      
      pval = anova.mmes(fit0, fit)[2,"PrChisq"]
      pval = as.numeric(tstrsplit(pval," ")[[1]])
      
      #anova(lm(emmean ~ haplo, dt))
      
      ### Aletrnatively, use  lme4breeding::lmebreed
      
      
      #fit0 = lmebreed(div ~ 1 + (1|id),
      #               relmat = list(id = Ghap ),
      #               verbose = FALSE, data=dt)
      
      #fit = lmebreed(div ~ haplo + (1|id),
      #         relmat = list(id = Ghap ),
      #         verbose = FALSE, data=dt, weights = weight_rotated)
      
      #anova(fit0,fit, test='LRT')
      
      data.frame(chrom=CHR,gpos=gpos[i],pval=pval)
      
    })
    
    res = do.call(rbind, res)
    #ggplot(hapgwas, aes(gpos, -log10(pval)))+geom_line()+facet_wrap(~chrom)
    
    
    hapgwas = rbind(hapgwas,res)
    
  }
  
  if(permutation){
    minp = min(hapgwas$pval)
    write(minp, file = 'analysis/QTLmapping/haplotypeGWAS_permutation_results_sizeConv_LOCO.txt', append = TRUE)
    
  }else{
    save(hapgwas, file = "analysis/QTLmapping/haplotypeGWAS_sizeConv_results_LOCO.Rdata")
  }
  
  
}











load("analysis/QTLmapping/haplotypeGWAS_sizeConv_results_LOCO.Rdata")

hapgwas = do.call(rbind, lapply(split(hapgwas, hapgwas$chrom), function(x){
  
  CHR = x$chrom[1]
  snps = as.data.frame(read_csv(paste0("genotypes/chr",CHR,"_beceiPanels_snps_founders&Rils_haplotypeResolved_V1.csv.gz")))
  
  y = x[c(1,1),]
  
  #y$ppos = range(snps$POS)
  y$gpos = range(snps$cM)
  y$pval = 1
  
  x=rbind(x,y)
  x=x[order(x$gpos),]
  x
  
}))

hapgwas$colorchrom = ifelse(hapgwas$chrom %in% c('II','IV','X'), 'grey30', 'black')


perm = unlist(read.table("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/becei_rils/analysis/QTLmapping/haplotypeGWAS_permutation_results_sizeConv_LOCO.txt", quote="\"", comment.char=""))

threshold_alpha0.1 = -log10(quantile(perm, 0.1))
threshold_alpha0.05 = -log10(quantile(perm, 0.05)) # 4.296709

#qtllim = range(subset(hapgwas, -log10(pval)>threshold_alpha0.05)$gpos)

p= ggplot(hapgwas)+
  #geom_hline(yintercept = threshold_alpha0.1, color = "#0073B3")+
  geom_hline(yintercept = threshold_alpha0.05, color = "#CC6600")+
  #geom_rect(data=data.frame(chrom="I"), aes(xmin = qtllim[1], xmax=qtllim[2],ymin=-Inf, ymax = Inf), fill = "#CC6600", color=NA, alpha=0.2)+
  geom_line(aes(gpos, -log10(pval)))+
  facet_grid(.~chrom, scales='free_x')+
  #scale_color_manual(values = unique(hapgwas$colorchrom), breaks =  unique(hapgwas$colorchrom))+
  theme_Publication3()+
  theme(panel.spacing = unit(0, "lines"))+
  coord_cartesian(ylim = c(0,threshold_alpha0.05+0.02))+
  scale_x_continuous(breaks = c(0,15,30, 45))+
  theme(legend.position = "none")+
  xlab("Genetic distance (cM)")



ggsave(p, file="figures/fig_gwas_haplotype_sizeConv_LOCO.png", width=3.6, height=1.15, dpi=1200)
ggsave(p, file="figures/fig_gwas_haplotype_sizeConv_LOCO.pdf", width=3.6, height=1.15, dpi=1200)



#hapgwas_conv = hapgwas










p=ggplot(hapgwas_div)+
  #geom_hline(yintercept = threshold_alpha0.1, color = "#0073B3")+
  geom_hline(yintercept = 4.20412 , color =  "#B56516")+
  geom_hline(yintercept = 4.296709 , color = "#F8D1A2")+
  #geom_rect(data=data.frame(chrom="I"), aes(xmin = qtllim[1], xmax=qtllim[2],ymin=-Inf, ymax = Inf), fill = "#CC6600", color=NA, alpha=0.2)+
  geom_line(aes(gpos, -log10(pval)))+
  geom_line(data = hapgwas_conv, aes(gpos, -log10(pval)), color = 'grey')+
  facet_grid(.~chrom, scales='free_x')+
  #scale_color_manual(values = unique(hapgwas$colorchrom), breaks =  unique(hapgwas$colorchrom))+
  theme_Publication3()+
  theme(panel.spacing = unit(0, "lines"))+
  coord_cartesian(ylim = c(0,-log10(min(hapgwas_div$pval))+0.02))+
  scale_x_continuous(breaks = c(0,15,30, 45))+
  theme(legend.position = "none")+
  xlab("Genetic distance (cM)")


ggsave(p, file="figures/fig_gwas_haplotype_sizeConv+Div_LOCO.png", width=3.6, height=1.15, dpi=1200)
ggsave(p, file="figures/fig_gwas_haplotype_sizeConv+Div_LOCO.pdf", width=3.6, height=1.15, dpi=1200)



