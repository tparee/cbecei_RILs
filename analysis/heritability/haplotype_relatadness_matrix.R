library(readxl)
library(readr)
library(reshape2)
library(data.table)
library(ggplot2)

PATH = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/CBCI_RILs/"
setwd(PATH)
WinSize_cM = 1

################################################################
################# HAPLOTYPE-BASED GRM ##########################
### Hickey et al., 2013 ( https://doi.org/10.1111/jbg.12020 )###

phasedGT_from_Haplo = function(haplo, phasedfounderhaplotypes){
  
  # haplo = data.frame containing a single RIL haplotype block (ex: subset(chrI_RILs_founderHaplotypeBlocks_v01.csv, rilname == 'xxxx'))
  # This works when the data.frame correponds homozygous diplotype ("g12") or the two phased haploid genome ('g1' & 'g2'), which are generated when heterozygsoity is detected 
  
  # phasedfounderhaplotypes = founder haplotypes (ex: chrI_beceiPanels_geno_founders_haplotypeResolved_v2.csv)
  
  # replaceMissingByObs:
  # The function is going to generate snps data from haplotype. Some SNPs can be NA because 
  # (i) they are missing in the phasedfounderhaplotypes, (ii) their is two possible founder at that position and their SNP differ
  # If, replaceMissingByObs is set as TRUE, the observed genotype will be added
  # In that case, the observed ril genotype must be provided in rilgenotype agrument (vector of genotype)
  
  HGT = matrix(NA, ncol=2, nrow = nrow(phasedfounderhaplotypes))
  for(i in 1:nrow(haplo)){
    x = haplo[i,]
    foudersx = unlist(strsplit( x$founder, "&"))
    win = (x$whichsnp1):(x$whichsnp2)
    
    hgti = do.call(cbind, lapply(foudersx, function(fx){
      fx = unlist(strsplit(fx, ";")) # founders names corresponding to this block
      fgtx = phasedfounderhaplotypes[win,fx] # founders genotypes
      
      if(length(fx)>1){ # of more than one possible founder, keep only shared snps
        div = apply(fgtx,1,sd)>0 # divergent snp
        fgtx=fgtx[,1]
        fgtx[div]=NA # shared snps
      }
      
      return(matrix(fgtx, ncol=1))
    }))
    
    if(ncol(hgti) == 1){hgti = hgti[,c(1,1)]}
    
    if(sum(HGT[win,] - hgti, na.rm=T) > sum(HGT[win,] - hgti[,c(2,1)], na.rm=T)){
      hgti = hgti[,c(2,1)]
    }

    # this part of the code is necessary because blocks overlap
    # if there is divergence with SNP inferred from previously estimated block, put as NA
    div = which(HGT[win,] != hgti)
    if(length(div)>0){hgti[div]=NA}
    HGT[win,] = hgti
  }
  
  return(HGT)
}


calculateHapSimilarity1 = function(x){
  # H1 method in Hickey et al., 2013
  # x = phasedrils[win,]
  
  afs = apply(x,1,mean, na.rm=T)
  pna = apply(x,1, function(i){mean(is.na(i))})
  pna2 = apply(x,2, function(i){mean(is.na(i))})
 
  rilnames = colnames(x)
  
  x = x[afs < 1 & afs > 0 & pna < 0.1,pna2 < 0.10]
  
  hapgroups =  cutree(hclust(dist(t(x))), h = 0)
  groupnames = names(table(hapgroups))
  
  
  hapLib = do.call(cbind, lapply(groupnames, function(thisgroup){
    gt =  x[,hapgroups == thisgroup]
    if(is.null(ncol(gt))){return(gt)}
    gt = round(apply(gt,1, mean, na.rm=T))
  }))
  
  hapLib[is.nan(hapLib)]=NA
  
  
  H1 = (nrow(hapLib) - as.matrix(dist(t(hapLib), method = "manhattan")))/nrow(hapLib)
  
  
  hapgroups = hapgroups[match(rilnames,names(hapgroups))]
  names(hapgroups) = rilnames
  list(hapgroups,H1)
}



# calculateHapSimilarity2 = function(x){
#   # H2 method in Hickey et al., 2013
#   # x = phasedrils[win,]
#   
#   afs = apply(x,1,mean, na.rm=T)
#   pna = apply(x,1, function(i){mean(is.na(i))})
#   pna2 =  apply(x,2, function(i){mean(is.na(i))})
#   
#   excludedRils = colnames(x)[which(pna2 >= 0.25)]
#   
#   x = x[afs < 1 & afs > 0 & pna < 0.1,pna2 < 0.25]
#   
#   hapgroups =  cutree(hclust(dist(t(x))), h = 0)
#   groupnames = names(table(hapgroups))
#   
#   hapLib = do.call(cbind, lapply(groupnames, function(thisgroup){
#     
#     gt =  x[,hapgroups == thisgroup]
#     if(is.null(ncol(gt))){return(gt)}
#     gt = round(apply(gt,1, mean, na.rm=T))
#   }))
#   
#   hapLib[is.nan(hapLib)]=NA
#   
#   coords = expand.grid(1:ncol(hapLib), 1:ncol(hapLib) )
#   coords = coords[coords$Var1 <= coords$Var2,]
#   
#   coords$hapsim = unlist(lapply(1:nrow(coords), function(i){
#     
#     hap1 = hapLib[,coords[i,1]]
#     hap2 = hapLib[,coords[i,2]]
#     
#     isna = is.na(hap1) | is.na(hap2)
#     hap1 = hap1[!isna]
#     hap2 = hap2[!isna]
#     
#     segments = rle(as.numeric(hap1 == hap2))
#     segments_size = segments$lengths[segments$values == 1]
#     sum(segments_size^2)
#   }))
#   
#   
#   H2 = xtabs(hapsim ~ Var1 + Var2, data = coords)
#   H2[lower.tri(H2)] <- t(H2)[lower.tri(H2)]
#   H2 = as.matrix(H2)
#   
#   attr(H2, "class") <- NULL
#   attr(H2, "call") <- NULL
#   
#   H2 = sqrt(H2 / max(H2, na.rm=TRUE))
#   
#   list(hapgroups,H2)
#   
# }



buildG = function(hapsim){
  
  hapgroups = hapsim[[1]]
  hapsim = hapsim[[2]]
  
  rilnames = unique(tstrsplit(names(hapgroups), ".g")[[1]])
  rilspairs = expand.grid(rilnames, rilnames, stringsAsFactors = FALSE)
  rilspairs = rilspairs[rilspairs$Var1 <= rilspairs$Var2, ]
  
  rilssimilarity = unlist(lapply(1:nrow(rilspairs), function(i){
    
    thispair = rilspairs[i,]
    ril1hap = hapgroups[names(hapgroups) %in% paste0(thispair[1], c(".g1",".g2"))]
    ril2hap = hapgroups[names(hapgroups) %in% paste0(thispair[2], c(".g1",".g2"))]
    
    S = (hapsim[ ril1hap[1], ril2hap[1] ] +
           hapsim[ ril1hap[2], ril2hap[2] ] +
           hapsim[ ril1hap[1], ril2hap[2] ] + 
           hapsim[ ril1hap[2], ril2hap[1] ])/2
    
    S
    
  }))
  
  rilspairs = as.data.frame(rilspairs)
  
  rilspairs$similarity  = rilssimilarity
  
  G = xtabs(similarity ~ Var1 + Var2, data = rilspairs)
  G[lower.tri(G)] <- t(G)[lower.tri(G)]
  G = as.matrix(G)
  
  attr(G, "class") <- NULL
  attr(G, "call") <- NULL
  
  G
  
}


for(CHR in c("I","II","III",'IV',"V",'X')){
  print(CHR)
  snpsx = as.data.frame(fread(paste0("genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz")))
  fgt = as.matrix(fread(paste0("genotypes/",CHR,"_becei_genotypes_founders.csv.gz")))
  rils = as.matrix(fread(paste0("genotypes/",CHR,"_becei_genotypes_RILs.csv.gz")))
  haplo = as.data.frame(fread(paste0("haplotypes/",CHR,"_rils_foundingHaplotypesBlocks.csv")))
  
  haplo$whichsnp1 = sapply(haplo$pos1, function(p){which.min(abs(snpsx$POS-p))})
  haplo$whichsnp2 = sapply(haplo$pos2, function(p){which.min(abs(snpsx$POS-p))})
  
  # Add x cM bins
  cuthere = seq(min(snpsx$cM), max(snpsx$cM),WinSize_cM)
  cuthere[1] = min(cuthere)-0.1
  cuthere[length(cuthere)] = max(cuthere)+0.1
  snpsx$win = cut(snpsx$cM, cuthere)

  phasedrils = do.call(cbind, lapply(unique(haplo$rilname), function(thisril){
    #print(thisril)
    hap = subset(haplo, rilname == thisril)
    gt_fromHap = phasedGT_from_Haplo(haplo=hap, phasedfounderhaplotypes=fgt)
    colnames(gt_fromHap) = paste0(thisril, c('.g1','.g2'))
    gt_fromHap
  }))
  
  # Across windows, calculate the haplotype similarity matrix
  # Then use them to build relationship G matrix between lines
  GMATS = lapply(split(snpsx, snpsx$win), function(win){
    #print(win$win[1])
    win = snpsx$ID %in% win$ID 
    # Find all of the unique haplotypes at that window (including recombinant ones)
    # and caluclate their similarity (proportion of shared allele)
    # The function return the haplotype similarity matrix (H mat) and the haplotype identity for each line
    Hmat = calculateHapSimilarity1(x=phasedrils[win,]) 
    
    # For each pair of lines, caluclate their relatedness based on their two haplotypes
    # That is "Gmat" in Hickey et al. 2013.
    Gmat = buildG(hapsim=Hmat)
    Gmat
  })
  
  # Average Gmat across all windows
  Garr <- abind::abind(GMATS, along = 3)    # stack into 3D array
  Ghap <- apply(Garr, c(1, 2), mean, na.rm=T) 
  
  write.csv(Ghap, file = paste0("analysis/heritability/haplotype_relatedness_matrix/chr",CHR, "_Ghap.csv"), row.names=F)
}





# weights = c(rep(50,5),25)
# weights = weights/sum(weights)
# Ghap=NULL
# for(CHR in c("I","II","III",'IV',"V",'X')){
#   g = fread(file = paste0("analysis/GWAS/chr",CHR, "_Ghap.csv"))
#   g = as.matrix(g)
#   rownames(g)=colnames(g)
#   g=g*weights[ which(c("I","II","III",'IV',"V",'X')==CHR) ]
#   if(is.null(Ghap)){Ghap=g}else{
#     Ghap=Ghap+g
#   }
# }
# 
# write.csv(Ghap, file = paste0("analysis/GWAS/haplotypeRelatednessMatrix.csv"),row.names = F)

Ghap = as.matrix(read_csv("analysis/GWAS/haplotypeRelatednessMatrix.csv"))
rownames(Ghap) = colnames(Ghap)

         
Ghap = reshape2::melt(Ghap)
Ghap$cross = paste0(tstrsplit(Ghap$Var1, '_')[[1]],tstrsplit(Ghap$Var2, '_')[[1]])
Ghap$cross[Ghap$cross == "BA"] = 'AB'
Ghap$Var1 = as.character(Ghap$Var1)
Ghap$Var2 = as.character(Ghap$Var2)
Ghap = Ghap[Ghap$Var1 <= Ghap$Var2,]

ggplot(Ghap, aes(cross,value, color = cross))+geom_boxplot()




genotypes =  as.matrix(read_csv("~/Downloads/beceiPanels_geno_RILs_pruned0.9999.csv"))
colnames(genotypes) = tstrsplit(colnames(genotypes),'-')[[1]]
afs = apply(genotypes, 1, function(x){mean(x[x==0 | x==1],na.rm=T)})
snps =  as.data.frame(read_csv("~/Downloads/beceiPanels_snps_pruned0.9999.csv"))

genotypes = (genotypes-0.5)*2
genotypes[genotypes == 0]=1


GT = t(genotypes)
GT[is.na(GT)]=0


A = sommer::A.mat(GT)


A = reshape2::melt(A)
A$cross = paste0(tstrsplit(A$Var1, '_')[[1]],tstrsplit(A$Var2, '_')[[1]])
A$cross[A$cross == "BA"] = 'AB'
colnames(A)[3] = 'gsim_VR' # Van raden
colnames(Ghap)[3] = 'gsim_hap' # Van raden

ggplot(A[A$Var1 != A$Var2,], aes(cross, value, color = cross))+geom_boxplot()
test = merge(A,Ghap)
ggplot(test[test$Var1 != test$Var2,], aes(gsim_VR, gsim_hap, color = cross))+geom_point(shape=1, alpha=0.5)

mean(test$gsim_VR)
mean(test$gsim_hap)



                