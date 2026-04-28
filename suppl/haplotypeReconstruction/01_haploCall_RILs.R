library(parallel)
library(readr)
setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/")
meta = read_csv("RILs_sequencing_metadata.csv")
###################################
##### Parameters & functions ######

source("haplotypeReconstruction/functions/functions_diplosearch.R")

for(CHR in c("I",'III', 'IV','V', 'X')){
print(CHR)
##########################
##### Files path #########

# The phased founder haplotype .Rdata file obtained through beceiFounders_phasing.R
phasedfoundersFile = paste0("haplotypeReconstruction/temp/phased_founders_base/",CHR,"_beceiFounderPhasedHaplotypes_base.Rdata")

# The stringent filtering genotype table used to infer phased founder haplotype
GTfile = paste0("genotype_stringentFiltering/",CHR,"_geno_becei_RILs_stringentFiltering.csv.gz")
infofile = paste0("genotype_stringentFiltering/",CHR,"_snps_becei_RILs&Pools_stringentFiltering.csv.gz")

###################################################################
######## STEP 1 : Haplotyping of founder blocks in RILs ###########
###################################################################

####### Import data #######
load(phasedfoundersFile) # Founder phased haplotypes
rils <- read_csv(GTfile) # RILs genotype Tables
rils = as.matrix(rils)
snps = read_csv(infofile) # snps info
####### Format data #######
mm = match(phasedfounderhaplotypes$POS, snps$POS) # Match positions 
rils=rils[mm,]
snps = phasedfounderhaplotypes[,1:5] # snps info
colnames(snps)=tolower(colnames(snps))

phasedfounderhaplotypes = phasedfounderhaplotypes[,6:11] # Keep only genotype info
foudernames = colnames(phasedfounderhaplotypes) # extract foundernames
phasedfounderhaplotypes = matrix(as.numeric(unlist(c(phasedfounderhaplotypes))), ncol=ncol(phasedfounderhaplotypes)) # format as numeric matrix
colnames(phasedfounderhaplotypes)=foudernames

####### haplotyping #######

#3093
#3094
rils = rils[,colnames(rils) %in% meta$rilname]

# "QG3143"
# "QG3144"
# "QG3146"
# "QG3157"
# "QG3172"
# "QG3204"
# "QG3220"
# "QG3227"
# "QG3249"
# "QG3251"
#  "QG3254"
# "QG3257"
# "QG3265"

thisril = "QG3500" # "QG3272"

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
  
  diplo$rilname = thisril
  colnames(diplo)[1:2] = c("pos1","pos2")
  diplo$pos1 = snps$pos[keep][diplo$pos1]
  diplo$pos2 = snps$pos[keep][diplo$pos2]
  return(diplo)
})


rilhaplotypes = do.call(rbind, rilhaplotypes)
save(rilhaplotypes, file=paste0("haplotypeReconstruction/temp/rils_halplotypes_stringentFiltering/",CHR,"_rils_haplotypes_stringentFiltering.Rdata"))
}
sort(table(rilhaplotypes$rilname))

subset(rilhaplotypes,  rilname == "QG3500")
table(rilhaplotypes$rilname)[which.max(table(rilhaplotypes$rilname))]
