library(data.table)
library(readr)

get.win = function(size, winsize, minsize){
  
  if(winsize < size){
    
    nblocks = size/winsize
    fullblocks = floor(nblocks)
    
    ## split column numbers into 'nblocks' groups
    SPLIT <- split(1:(fullblocks*winsize), rep(1:fullblocks, each = winsize))
    if(nblocks>fullblocks){ 
      dblock = (nblocks-fullblocks)*winsize
      SPLIT[[length(SPLIT)+1]] = 1:dblock + max(fullblocks*winsize)
    }
    
    if(length(SPLIT[[length(SPLIT)]]) < minsize & length(SPLIT) > 1){
      SPLIT[[length(SPLIT)-1]] = c(SPLIT[[length(SPLIT)-1]],SPLIT[[length(SPLIT)]])
      SPLIT = SPLIT[1:(length(SPLIT)-1)]
    }
    
    
  }else{
    
    SPLIT = list(1:size)
  }
  
  return(SPLIT)
}

#CHR = "I"
setwd('/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/')

for(CHR in c("I","II","III","IV","V","X")){
  print(CHR)
rilhaplotypes = as.data.frame(readr::read_csv(paste0("haplotypes/",CHR,"_rils_foundingHaplotypesBlocks2.csv")))
PL <- as.data.frame(fread(paste0("suppl/private/alignment&variantCalling/VCF/RILs/unfiltered/",CHR,"_becei_variants_PL.csv.gz")))
#GT <- as.data.frame(fread(paste0("suppl/private/alignment&variantCalling/VCF/RILs/unfiltered/",CHR,"_becei_variants_GT.csv.gz")))
PL$ID = paste0(PL$CHROM, ":", PL$POS, ":", PL$REF, ":",PL$ALT)

phasedfounderhaplotypes = as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_genotypes_founders.csv")))
snps =  as.data.frame(fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_variantInfo_founders&Rils.csv")))
snps$ID = paste0(snps$CHROM, ":", snps$POS, ":", snps$REF, ":",snps$ALT)
keep = match(snps$ID, PL$ID)
PL = PL[keep,]
#GT = GT[keep,]


rilnames = sort(unique(rilhaplotypes$rilname))

ngt = unlist(lapply(strsplit(PL[,paste0(rilnames[1],"_PL")], ","), length))
if(sum(ngt != 3) > 0){stop("code is only robust for biallelic snps")}

geno_imputed = matrix(NA, nrow = nrow(snps), ncol = length(rilnames))
colnames(geno_imputed) = rilnames

#GT$ID = paste0(GT$CHROM, ":", GT$POS, ":", GT$REF, ":",GT$ALT)
#GT = GT[match(snps$ID, GT$ID),]
#View(cbind(rilgt_imputed, GT[,thisril]))

for(thisril in rilnames){
  print(thisril)
  plmatrix = do.call(rbind, lapply(strsplit(PL[,paste0(thisril,"_PL")], ","), as.numeric))
  plmatrix = plmatrix-apply(plmatrix,1,min)
  #  PL (Phred-scaled Likelihoods) = -10 x log10 P(Genotype | Data) + C
  #  Can be transformed back into likelohood => GL = 10^(PL/10); GL = GL/sum(GL) (i.e., needs to be standardized because of the constant C)
  glmatrix = 10^(-plmatrix / 10)
  glmatrix = glmatrix / apply(glmatrix,1,sum)
  haplo = subset(rilhaplotypes, rilname == thisril)
  
  expected_geno = rep(NA, nrow(snps))
  for(i in 1:nrow(haplo)){
    #print(i)
    pos1 = haplo[i,]$pos1
    pos2 = haplo[i,]$pos2
    diplo_id = haplo[i,]$founder
    win = which(snps$POS >= pos1 & snps$POS <= pos2)
    
    diplo_id = unlist(strsplit(diplo_id, "&"))
    gt = do.call(cbind, lapply(1:length(diplo_id), function(n){
      x = diplo_id[n]
      x = phasedfounderhaplotypes[win,unlist(strsplit(x, ";"))]
      if(!is.null(dim(x))){ x = apply(x, 1, mean)}
      x[x > 0 & x < 1] = NA
      x
    }))
    
    gt = apply(gt, 1, mean)
    
    div = which(expected_geno[win] != gt)
    expected_geno[win] = gt
    if(length(div) > 0){gt[div] = NA}
  }
  
  selected_PL = plmatrix[cbind(1:nrow(snps),match(expected_geno,c(0,0.5,1)))]
  
  
  windows = get.win(size = nrow(snps), winsize=1000, minsize=500)
  
  rilgt_imputed = do.call(c,lapply(windows, function(win){
    # For a window, we estimate the probability of the expected geno based on inferred haplotype
    # I do this by averaging the PL of expected geno, which will increases if the observed data diverge from the expected geno
    log10L_per_snp <- -mean(selected_PL[win], na.rm=T) / 10 # PL to log10
    
    # Even when expected genotype perfectly matches the reality, some divergence is expected by chance
    # This probability is given my the genotype likelihood, e.g., 0/0: 0.1; 1/1: 0.9 means that there is 10% chance that the most supported genotype is wrong
    # So the expected PL is given by sum(GL*PL)
    expected_log10L_per_snp = - mean(apply(plmatrix[win,]*glmatrix[win,], 1, sum), na.rm=T) / 10
    # Transform back to likelihood
    L_per_snp <- 10^( min( c(log10L_per_snp - expected_log10L_per_snp, 0)))
    if(is.nan(L_per_snp)){L_per_snp=1/3}
    Lmax = min(L_per_snp,0.999) # add a ceiling to the maximal confidence we have in out inferred haplotype so it is never 100%
    
    # Prior (from haplotype) is the Lmax for the expected genotype and (1-Lmax)/2 for the two other genotypes
    prior = matrix((1-Lmax)/2, nrow=length(win),ncol=3)
    prior[cbind(1:length(win), match(expected_geno[win],c(0,0.5,1)))] = Lmax
    prior = prior/apply(prior,1,sum)
    
    # Compute the posterior genotype likelihood based on prior and observation
    log_prior <- log(prior)
    log_likelihood <- -plmatrix[win,] / 10 * log(10)
    log_post <- log_prior + log_likelihood
    log_post <- log_post - max(log_post)
    posterior <- exp(log_post)
    posterior <- posterior / apply(posterior, 1, sum)
    
    # output the most likely genotypes for the window
    bestgeno = c(0,0.5,1)[apply(posterior,1, which.max)]
    return(bestgeno)
    
  }))
  
  geno_imputed[,thisril] = rilgt_imputed
}

write_csv(as.data.frame(geno_imputed), paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/genotypes/",CHR,"_becei_genotypes_RILs.csv"))

}

#cbind(GT[,thisril], rilgt_imputed)
