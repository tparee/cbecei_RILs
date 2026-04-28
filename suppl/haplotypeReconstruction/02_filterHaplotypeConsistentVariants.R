library(data.table)
genodistMax = 3
max_pNA = 0.2

CHR = 'X'

for(CHR in c("I",'III', 'IV','V', 'X')){
  print(CHR)
###########################
### IMPORT DATA  ##########

setwd("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/")
load(file=paste0("haplotypeReconstruction/temp/rils_halplotypes_stringentFiltering/",CHR,"_rils_haplotypes_stringentFiltering.Rdata"))
load(paste0("haplotypeReconstruction/temp/phased_founders_base/",CHR,"_beceiFounderPhasedHaplotypes_base.Rdata")) # Founder phased haplotypes
snps = phasedfounderhaplotypes[,1:4] # snps info
snps$pos=snps$POS

genoall = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/private/alignment&variantCalling/VCF/RILs/unfiltered/",CHR,"_becei_variants_GT.csv.gz")))
snpsall=genoall[,1:6]
genoall=as.matrix(genoall[,-1:-6])

#table(genoall[snpsall$POS == 308200,])

###########################
### Little format change ##

rilhaplotypes = subset(rilhaplotypes, !(ishet))
rilhaplotypes=rilhaplotypes[,-3]
rilhaplotypes = do.call(rbind, lapply(split(rilhaplotypes , rilhaplotypes$rilname), function(x){
  
  #print(x$rilname[1])
  #x = split(rilhaplotypes , rilhaplotypes$rilname)[[3]]
  if(nrow(x) == 1){return(x)}
  
  xsplit = NULL
  for(i in 1:nrow(x)){
    #print(i)
    h = i-1
    j = i+1
    xi_cut = x[i,]
    xi_cut[,"pos1"] = max(c(x[h,"pos2"],xi_cut[,"pos1"]),na.rm=T)
    xi_cut[,"pos2"] = min(c(x[j,"pos1"],xi_cut[,"pos2"]), na.rm=T)
    if(i < nrow(x)){
      xij_overlap = data.frame(pos1 = x[j,]$pos1, pos2 = x[i,]$pos2,
                               founder =  paste(sort(unique(unlist(sapply(x[c(i,j),]$founder, function(y){unlist(strsplit(y, ";"))})))), collapse=";"),
                               rilname = x$rilname[1])
      
      if( sign(diff(unlist(xij_overlap[,1:2]))) == -1){xij_overlap = NULL}
    }else{
      xij_overlap = NULL
    }
    
    
    xi_cut =rbind(xi_cut,xij_overlap)
    xsplit = rbind(xsplit, xi_cut)
    
  }
  
  return(xsplit)
  
}))


################################################################################
#### Add any genotypes that is consitent with haplotype structure in RILs ######
#### i.e., Monomorphic within founding haplotype indentity at that position ####

supplFounderGeno = do.call(rbind, lapply(1:nrow(snpsall), function(i){
  
  # for each snps that were previously filtered out
  # we are going to see if there are consistent with the known linkage blocks
  
  #if(i %% 1000 == 0) print(i)
  #print(i)
  
  ppos = snpsall[i,]$POS
  
  startpos = min(rilhaplotypes$pos1)
  endpos = max(rilhaplotypes$pos2)
  if(ppos<startpos){ppos = startpos}
  if(ppos>endpos){ppos = endpos}
  
  # All the haplotype blocks among rils that encompass the target position
  blocks = rilhaplotypes[(rilhaplotypes$pos1 <= ppos & rilhaplotypes$pos2 >= ppos),]
  
  # Add the genotype observed in the rils at the target position
  blocks$genotype = unlist(c(genoall[i,match(blocks$rilname, colnames(genoall))]))
  blocks$genotype[is.na(blocks$genotype)] = "./."
  if( mean(blocks$genotype == "./.") > 0.5){return(NULL)}
  # Verify there is not too much NA
  #if( sum(is.na(blocks$genotype))/nrow(blocks) > fNAmax ) return(NULL)
  
  # Verify that there is only one allele per founder haplotype 
  nperhap = table(blocks$founder, blocks$genotype)
  if( !("./." %in% colnames(nperhap)) ){nperhap = cbind(`./.` = rep(0, nrow(nperhap)), nperhap)}
  founders = strsplit(rownames(nperhap),";")
  #genovalue = as.numeric(colnames(nperhap))
  genovalue = colnames(nperhap)
  
  # Ensure there is variation 0/1 at this snp (if only 0/0.5, not good)
  #if(sum(c(0,1) %in% genovalue)<2) return(NULL) 
  isuniquefounder = unlist(lapply(founders, length))==1
  wmiss = which(genovalue == "./.")
  # A
  foundersgeno = do.call(rbind, lapply(which(isuniquefounder), function(f){
    
    x=nperhap[f,] # genotype distribution for this founder
    if(sum(x[-wmiss]) == 0 & x[wmiss] > 0){
      return(data.frame(founders=rownames(nperhap)[f], genotype = "./.", genodist=0, nmiss = 0, ntot = sum(x)))
    }
    
    missing = as.numeric(x[wmiss])
    wgeno = which.max(x[-wmiss]) # The most frequent genotype = the genotype of this founder
    gdist= sum(x[-wmiss][-wgeno]) # number of individual with a different genotype than expected
    missing = x[wmiss]
    #whet = which(names(x)==0.5)
    #nhet = ifelse(length(whet)==1, x[whet], 0)
    data.frame(founders=rownames(nperhap)[f], genotype = genovalue[-wmiss][wgeno], genodist=gdist, nmiss = missing, ntot = sum(x))
    
  }))
  
  
  #possibility that two or more founders are the same and thus missing from unique founders
  foundernames = c("FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2")
  wmissingfounders = which(!(foundernames %in% foundersgeno$founders))
  
  if(length(wmissingfounders)>1){
    
    missingfounders = foundernames[wmissingfounders]
    
    corresp = lapply(1:nrow(nperhap), function(f){
      x = unlist(strsplit(rownames(nperhap)[f],";"))
      
      if(sum(x %in% missingfounders) == length(x)){
        return(list(missingfounders = missingfounders[missingfounders %in% x],
                    whichrow = f))
      }
    })
    
    corresp = corresp[!unlist(lapply(corresp, is.null))]
    
    if(sum(missingfounders %in% unlist(lapply(corresp, function(x){x$missingfounders}))) != length(missingfounders)){stop(i)}
    
    
    for(f in corresp){
      x = nperhap[f$whichrow,]
      nf = length(f$missingfounders)
      
      if(sum(x[-wmiss]) == 0 & x[wmiss] > 0){
        foundersgeno = rbind(foundersgeno,
                             data.frame(founders=f$missingfounders, genotype = "./.", genodist=0, nmiss = 0, ntot = sum(x)/nf))
        isuniquefounder[f$whichrow]=T
      }else{
        missing = as.numeric(x[wmiss])
        wgeno = which.max(x[-wmiss]) # The most frequent genotype = the genotype of this founder
        gdist= sum(x[-wmiss][-wgeno]) # number of individual with a different genotype than expected
        
        #whet = which(names(x)==0.5)
        #nhet = ifelse(length(whet)==1, x[whet], 0)
        foundersgeno = rbind(foundersgeno,
                             data.frame(founders=f$missingfounders, genotype = genovalue[-wmiss][wgeno], genodist=gdist/nf, nmiss = missing/nf, ntot = sum(x)/nf))
        isuniquefounder[f$whichrow]=T
        
      }
    }
  }
  # Sum the genodist observed in haplotype blocks corresponding to a unique founder
  # and the ones from the blocks corresponding to blocks which can correspond do several founders
  genotypedistance = try(sum(foundersgeno$genodist) + sum(unlist(lapply(which(!isuniquefounder), function(f){
    
    x=nperhap[f,]
    thesefounders = founders[[f]]
    thesefounders = foundersgeno[match(thesefounders, foundersgeno$founders),]
    thesefounders$genotype
    isnotfoundergeno = !(names(x) %in% as.character(thesefounders$genotype))
    isnotfoundergeno[wmiss] = F
    
    x[isnotfoundergeno]
    
  }))))
  
  
  
  if(class(genotypedistance)[1]=='try-error'){return(NULL)}
  
  # if genotype distance over threshold, snp is bad
  if(genotypedistance > genodistMax){ return(NULL)}
  foundersgeno$genotype[(foundersgeno$nmiss / foundersgeno$ntot) > max_pNA] = NA
  # if snp is good, return the ganotype for each founder
  foundersgeno = foundersgeno$genotype[match(c("FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2"),foundersgeno$founders)]
  if( sum(is.na(foundersgeno)) > 2){return(NULL)}
  if( mean(foundersgeno == "0/0", na.rm = T) == 1){return(NULL)}
  ishet = grepl("0/",foundersgeno) & !grepl("/0",foundersgeno)
 
  keep = sum(foundersgeno %in% c("0/0",'1/1'),na.rm=T) == length(foundersgeno) &
    sum(is.na(foundersgeno)) == 0 & sum(ishet | foundersgeno == "./.", na.rm = T) == 0 &
    sum(foundersgeno == "1/1", na.rm = T) < 6 & sum(foundersgeno == "0/0", na.rm = T) < 6
    
  
  #if(sum(is.na(foundersgeno))>0){return(NULL)}
  #if(sum(foundersgeno) == 0 | sum(foundersgeno==0.5)>0){return(NULL)}
  
  foundersgeno = c(i, keep , foundersgeno)
  return(foundersgeno)
  
}))


supplFounderGeno=as.data.frame(supplFounderGeno)
colnames(supplFounderGeno) = c("whichsnp","keep", "FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2")
supplFounderGeno = cbind(snpsall[as.numeric(supplFounderGeno$whichsnp),], supplFounderGeno[,2:8])
supplFounderGeno$keep = as.logical(supplFounderGeno$keep) & !grepl(",",supplFounderGeno$ALT)

nhapconsistent = nrow(supplFounderGeno)
nkeep = sum(as.logical(supplFounderGeno$keep))
nmultiallelic = sum(grepl(",",supplFounderGeno$ALT))
ndeletion =  sum(apply(supplFounderGeno[,c("FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2")],1,function(foundersgeno){sum(foundersgeno == "./.", na.rm = T) > 0}))
nfixed = sum(apply(supplFounderGeno[,c("FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2")],1,function(foundersgeno){sum(foundersgeno == "1/1", na.rm = T) == 6 | sum(foundersgeno == "0/0", na.rm = T) == 6}))

diagnosis = data.frame(chrom=CHR,
                       nhapconsistent,
                       nkeep,
                       nmultiallelic,
                       ndeletion,
                       nfixed)

print(diagnosis)

#######################
####### SAVE ##########

write.csv(supplFounderGeno, paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/private/haplotypeConsistentGenotypeTables/",CHR,"_becei_haplotypeConsistentVariants.csv"),row.names = F)
#supplFounderGeno = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/private/haplotypeConsistentGenotypeTables/",CHR,"_becei_haplotypeConsistentVariants.csv")))
supplFounderGeno2 = supplFounderGeno[as.logical(supplFounderGeno$keep), -7]

snps = supplFounderGeno2[,1:6]
geno = supplFounderGeno2[,c("FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2")]
geno[geno == '1/1'] = "1"
geno[geno == '0/0'] = "0"
geno[geno == '0/1'] = "0.5"
geno = do.call(cbind, lapply(1:ncol(geno), function(j){
  x= as.numeric(geno[,j])
  matrix(x, ncol=1)
}))
colnames(geno) = c("FA.g1","FA.g2","FB.g1","FB.g2","FM.g1","FM.g2")

map = as.data.frame(read.table('/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/cbecei_geneticMap.txt'))

snps$cM = approx(subset(map, chrom == CHR)$ppos, subset(map, chrom == CHR)$gpos, xout = snps$POS)$y
mininferred = min(which(!is.na(snps$cM)))
maxinferred = max(which(!is.na(snps$cM)))
if(maxinferred > 1){snps$cM[1:mininferred] = snps$cM[mininferred]}
if(maxinferred < nrow(snps)){snps$cM[maxinferred:nrow(snps)] = snps$cM[maxinferred]}

write.csv(geno, paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_genotypes_founders.csv"),row.names = F)
write.csv(snps, paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv"),row.names = F)
}
