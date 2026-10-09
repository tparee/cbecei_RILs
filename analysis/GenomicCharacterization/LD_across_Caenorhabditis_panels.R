library(readxl)
library(data.table)
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cebcei_RILs/")
source("utils.R")

doLD = function(gt,cM, cM_range = c(0.9,1.1)){
  if(length(cM) > 5000){
    set.seed(123)
    ix = sort(sample(1:nrow(gt),5000))
    gt = gt[ix,]
    cM = cM[ix]
  }
  # mean impute missing snps
  gt = do.call(rbind, lapply(1:nrow(gt), function(i){
    out = as.numeric(gt[i,])
    out[is.na(out)] = mean(out, na.rm=T)
    out
  }))
  
  gdis = as.matrix(dist(cM))
  r2 = cor(t(gt))^2
  meanr2 = mean(r2[which(gdis > cM_range[1] & gdis < cM_range[2], arr.ind = TRUE)] )
  return(meanr2)
}

###########################################
########## Linkage in the beRILs ##########
berils = NULL
for(CHR in c("I","II",'III',"IV",'V',"X")){
  geno = fread(paste0(file="genotypes/",CHR,"_becei_genotypes_RILs.csv.gz"))
  geno = as.matrix(geno)
  snps = fread(paste0(file="genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv.gz"))
  r2_1cM = doLD(gt=as.matrix(geno[snps$CHR == CHR,]),
                cM = snps$cM[snps$CHR == CHR],
                cM_range = c(0.9,1.1)) 
  
  berils = rbind(berils, data.frame(chrom = CHR, r2_1cM))
}

#chrom    r2_1cM
#1     I 0.2078993
#2    II 0.2426550
#3   III 0.2260950
#4    IV 0.3268441
#5     V 0.1881509
#6     X 0.2713615

mean(berils$r2_1cM) # 0.2438343


######################################
########## Linkage in CeMEE ##########

geno = fread("genotypes/private/otherPanels/CeMEEv2_RIL_geno.csv.gz")
snps = fread("genotypes/private/otherPanels/CeMEEv2_RIL_snps_ws245.csv.gz")
snps = do.call(rbind, lapply(split(snps, snps$chrom), function(x){
  x$cM = ifelse(x$chrom[1] == "X", 25, 50)*x$cM/max(x$cM)
  x
  }))


cemee = do.call(rbind, lapply(unique(snps$chrom), function(CHR){
  
  r2_1cM = doLD(gt=as.matrix(geno[snps$chrom == CHR,]),
                cM = snps$cM[snps$chrom == CHR],
                cM_range = c(0.45,0.55)) 
  
  data.frame(chrom = CHR, r2_1cM)
}))

#chrom     r2_0.5cM
#1     I 0.09436691
#2    II 0.06570146
#3   III 0.05008462
#4    IV 0.10007950
#5     V 0.08180788
#6     X 0.02963895

#chrom      r2_1cM
#1     I 0.022997133
#2    II 0.015197399
#3   III 0.015229639
#4    IV 0.015141915
#5     V 0.023412352
#6     X 0.007901381

mean(cemee$r2_1cM) # 0.01664664

###################################################
########## Linkage in Snoek,2019, mpRILs ##########

# Load C. elegans linkage maps
load("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/genotypes/private/otherPanels/cemap.RData")
#cemap$chrom = c("I","II","III","IV","V","X")[match(cemap$chrom, 1:6)]
#cemap = do.call(rbind, lapply(split(cemap, cemap$chrom), function(x){
#  x$genetic = ifelse(x$chrom[1] == "X", 25, 50)*x$genetic/max(x$genetic)
#  x
#}))
#save(cemap, file = "/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/genotypes/private/otherPanels/cemap.RData")

snps <- as.data.frame(read_excel("genotypes/private/otherPanels/12915_2019_642_MOESM1_ESM.xlsx", sheet = 4))
geno <- as.data.frame(read_excel("genotypes/private/otherPanels/12915_2019_642_MOESM1_ESM.xlsx", sheet = 5))
geno = geno[c(-1,-2),-1]
snps = do.call(rbind, lapply(split(snps, snps$ChrID), function(x){
  map = subset(cemap, chrom == x$ChrID[1])
  x$genetic = approx(x=map$cpos, y=map$genetic, xout = x$Position)$y
  x
}))


snoek = do.call(rbind, lapply(unique(snps$ChrID), function(CHR){
 
  r2_1cM = doLD(gt=as.matrix(geno[snps$ChrID == CHR,]),
                cM = snps$genetic[snps$ChrID == CHR],
                cM_range = c(0.9,1.1)) 
  
  data.frame(chrom = CHR, r2_1cM)
}))

#chrom    r2_1cM
#1     I 0.6727515
#2    II 0.8070851
#3   III 0.4327150
#4    IV 0.6732177
#5     V 0.6812258
#6     X 0.5904542

mean(snoek$r2_1cM) # 0.6429082



#########################################
###### Linkage in Brady, 2019 ###########
# devtools::install_github("AndersenLab/linkagemapping")
library("linkagemapping")
# Cannot load the full cross
#load_cross_obj("N2xCB4856cross_full") # cannot open URL 'https://storage.googleapis.com/linkagemapping/data/N2xCB4856cross_full2.Rda': HTTP status was '404 Not Found'

# So using the reduced marker version; should be enough for our goal
load_cross_obj("N2xCB4856cross")
load_cross_obj("N2xCB4856map")
load_cross_obj("N2xCB4856markers")

N2xCB4856markers$position = as.numeric(N2xCB4856markers$position)

brady = do.call(rbind, lapply(split(N2xCB4856markers, N2xCB4856markers$chr.roman), function(x){
  CHR = x$chr.roman[1]
  bradymap = unlist(N2xCB4856map[[CHR]])
  x$genetic = bradymap[match(x$marker, names(bradymap) )]
  x$genetic = x$genetic-min(x$genetic,na.rm=T)
  x$genetic = ifelse(CHR == "X", 25, 50)*x$genetic/max(x$genetic,na.rm=T)
  
  x = x[order(x$genetic),]
  
  #ggplot(x, aes(position, genetic))+geom_point()
  
  gt = N2xCB4856cross$geno[[CHR]]$data
  gt = t(gt[,match(x$marker,colnames(gt))])
  
  r2_1cM = doLD(gt=gt,
                cM = x$genetic,
                cM_range = c(0.9,1.1)) 
  data.frame(chrom = CHR, r2_1cM)
}))

#chrom    r2_1cM
#I       I 0.7742110
#II     II 0.8112385
#III   III 0.8084757
#IV     IV 0.7896830
#V       V 0.7972049
#X       X 0.7068866

mean(brady$r2_1cM) # 0.7812833


###########################################
###### Linkage in Stevens, 2022 ###########

load("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/genotypes/private/otherPanels/CB_genetic_map.Rda")
cbmap = fread("genotypes/private/otherPanels/c_briggsae_genetic_distances.txt")
colnames(cbmap) = c("marker", "chr","genetic")


stevens = do.call(rbind, lapply(1:6, function(CHR){
  
  x = subset(cbmap, chr == CHR)
  x$genetic = x$genetic-min(x$genetic,na.rm=T)
  x$genetic = ifelse(CHR == 6, 25, 50)*x$genetic/max(x$genetic,na.rm=T)
  x$position = as.numeric(tstrsplit(x$marker, "_")[[2]])
  x = x[order(x$genetic),]
  
  #ggplot(x, aes(position, genetic))+geom_point()
  
  gt = CB_newmap_manual$geno[[CHR]]$data
  gt = t(gt[,match(x$marker,colnames(gt))])
  
  r2_1cM = doLD(gt=gt,
                cM = x$genetic,
                cM_range = c(0.9,1.1)) 
  data.frame(chrom = CHR, r2_1cM)
}))

#chrom    r2_1cM
#1     1 0.8527137
#2     2 0.8926016
#3     3 0.9153223
#4     4 0.9037254
#5     5 0.8785468
#6     6 0.8575648

mean(stevens$r2_1cM) # 0.8834124


########################################
###### Linkage in Ross, 2011 ###########

ross = read_xlsx("genotypes/private/otherPanels/pgen.1002174.s004.xlsx")
ross = ross[-1,]
snps = ross[,c(4,8,9)]
colnames(snps) = c("chrom","pos","genetic")
snps$genetic = as.numeric(snps$genetic)
snps$pos = as.numeric(snps$pos)
geno = ross[,13:ncol(ross)]
geno = as.matrix(geno)
geno[geno=="A"] = 0
geno[geno=="B"] = 1
geno[geno=="H"] = 0.5

ross = do.call(rbind, lapply(unique(snps$chrom)[-7], function(CHR){
  cM = snps$genetic[snps$chrom == CHR]
  cM = ifelse(CHR == 6, 25, 50) * cM/max(cM)
  
  r2_1cM = doLD(gt=as.matrix(geno[snps$chrom == CHR,]),
                cM = cM,
                cM_range = c(0.9,1.1)) 
  
  data.frame(chrom = CHR, r2_1cM)
}))

#chrom    r2_1cM
#1     1 0.8568589
#2     2 0.8839434
#3     3 0.8306192
#4     4 0.8403623
#5     5 0.8604835
#6     6 0.8609557

mean(ross$r2_1cM) # 0.8555372



########################################
###### Noble, him-5 rils, 2011 ###########
load("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/genotypes/private/otherPanels/cemap.RData")
him5 = read_xlsx("genotypes/private/otherPanels/mmc2.xlsx", sheet = 4)
him5=him5[,-1]
him5 = as.data.frame(him5)
snps = data.frame(chrom = as.numeric(him5[1,]),pos = as.numeric(colnames(him5)))
snps$chrom = c("I","II",'III',"IV",'V',"X")[snps$chrom]

snps = do.call(rbind, lapply(split(snps, snps$chrom), function(x){
  map = subset(cemap, chrom == x$chrom[1])
  x$genetic = approx(x=map$cpos, y=map$genetic, xout = x$pos)$y
  x
}))

geno = him5[c(-1,-2),]
geno[geno=="CB4856"] = 1
geno[geno=="AB2"] = 0
geno=t(as.matrix(geno))

him5 = do.call(rbind, lapply(unique(snps$chrom), function(CHR){
  
  r2_1cM = doLD(gt=as.matrix(geno[snps$chrom == CHR,]),
                cM = snps$genetic[snps$chrom == CHR],
                cM_range = c(0.9,1.1)) 
  
  data.frame(chrom = CHR, r2_1cM)
}))

#chrom    r2_1cM
#1     I 0.8319088
#2    II 0.8606032
#3   III 0.7844588
#4    IV 0.9213172
#5     V 0.8519228
#6     X 0.9273619


mean(him5$r2_1cM) # 0.8629288


########################################
###### Noble, 2021, tropicalis #########

load("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/cbecei_RILs/genotypes/private/otherPanels/c_tropicalis_flye_geneticMap.rda")

tropi = do.call(rbind, lapply(1:6, function(CHR){
  
  cM = unlist(c(cross[["geno"]][[CHR]][["map"]]))
  cM = ifelse(CHR == 6, 25, 50) * cM/max(cM)
 
  
  gt = cross[["geno"]][[CHR]][["data"]]
  gt = t(gt[,match(names(cM),colnames(gt))])
  
  r2_1cM = doLD(gt=gt,
                cM = cM,
                cM_range = c(0.9,1.1)) 
  data.frame(chrom = CHR, r2_1cM)
}))

#chrom    r2_1cM
#1     1 0.9044025
#2     2 0.9534557
#3     3 0.8213583
#4     4 0.9259777
#5     5 0.8783667
#6     6 0.9190738

mean(tropi$r2_1cM) # 0.9004391
