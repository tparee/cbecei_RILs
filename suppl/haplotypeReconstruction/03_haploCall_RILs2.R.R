library(parallel)
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/")
meta = read_csv("suppl/RILs_sequencing_metadata.csv")
###################################
##### Parameters & functions ######

CHR = "III" # Chromosome of interest

source("suppl/haplotypeReconstruction/functions/functions_diplosearch.R")
for(CHR in c("X")){
##########################
##### Files path #########



phasedfounderhaplotypes = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_genotypes_founders.csv")))
snps =  as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv")))
snps$pos = snps$POS
geno = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/private/alignment&variantCalling/VCF/RILs/unfiltered/",CHR,"_becei_variants_GT.csv.gz")))
geno = geno[match(snps$POS,geno$POS),which(colnames(geno) %in% meta$rilname)]
geno=as.matrix(geno)
geno[geno == '1/1'] = "1"
geno[geno == '0/0'] = "0"
geno[geno == '0/1'] = "0.5"
rils = do.call(cbind, lapply(1:ncol(geno), function(j){
  x= as.numeric(geno[,j])
  matrix(x, ncol=1)
}))

colnames(rils) = colnames(geno)
rm(geno)




thisril = "QG3161" # "QG3272"


rilhaplotypes = lapply(colnames(rils), function(thisril){
  print(thisril)
  if(subset(meta, rilname == thisril)$panel == "α") foundernames = c("FA.g1","FA.g2", "FM.g1", "FM.g2")
  if(subset(meta, rilname == thisril)$panel == "β") foundernames = c("FB.g1","FB.g2", "FM.g1", "FM.g2")
  
  ril = rils[,thisril]
  keep = which(!is.na(ril))
  founders =phasedfounderhaplotypes[keep,foundernames]
  ril = ril[keep]
  
  diplo = try(diplosearch(ril = ril, founders=founders))
  if(class(diplo) == 'try-error'){
    print(paste0("failed: ",thisril))
    return(NULL)
  }
  
  runs = rle(as.numeric(as.factor(diplo$founder)))$lengths
  runs = unlist(lapply(1:length(runs), function(r){rep(r,runs[r])}))
  diplo = do.call(rbind, lapply(split(diplo, runs), function(x){
    if(nrow(x) == 1){return(x)}else{
      out = x[1,]
      out[,'whichsnp2'] = x$whichsnp2[nrow(x)]
      return(out)
    }
  }))
  
  diplo$rilname = thisril
  #diplo
  #length(ril)
  colnames(diplo)[1:2] = c("pos1","pos2")
  diplo$pos1 = snps$pos[keep][diplo$pos1]
  diplo$pos2 = snps$pos[keep][diplo$pos2]
  return(diplo)
})


rilhaplotypes = do.call(rbind, rilhaplotypes)
readr::write_csv(rilhaplotypes,file=paste0("haplotypes/",CHR,"_rils_foundingHaplotypesBlocks.csv"))
}




# I: QG3363
# X: QG3363 

# CHR = 'V'
# 
# for(CHR in c("I", "II", "III","IV","V","X")){
#   rilhaplotypes = as.data.frame(readr::read_csv(paste0("haplotypes/",CHR,"_rils_foundingHaplotypesBlocks.csv")))
#   sort(table(rilhaplotypes$rilname))
#   subset(rilhaplotypes, rilname == "QG3500")
#   rilhaplotypes = do.call(rbind, lapply(split(rilhaplotypes, rilhaplotypes$rilname), function(diplo){
#     runs = rle(as.numeric(as.factor(diplo$founder)))$lengths
#     runs = unlist(lapply(1:length(runs), function(r){rep(r,runs[r])}))
#     diplo = do.call(rbind, lapply(split(diplo, runs), function(x){
#       if(nrow(x) == 1){return(x)}else{
#         out = x[1,]
#         out[,'pos2'] = x$pos2[nrow(x)]
#         return(out)
#       }
#     }))
#     return(diplo)
#   }))
#   rilhaplotypes$chrom = CHR
#   readr::write_csv(rilhaplotypes,file=paste0("haplotypes/",CHR,"_rils_foundingHaplotypesBlocks.csv"))
# }

