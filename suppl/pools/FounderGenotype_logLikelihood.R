library(data.table)
library(parallel)
library(readr)

#for(CHR in c("I","II","III",'IV',"V",'X')){
for(CHR in c("III")){
  print(CHR)
# ====================================
# MULTI-SAMPLE FREQUENCIES LIKELIHOOD 

# Function to compute log-likelihood for binomial
loglik_binom <- function(alt_count, total_count, p){
  p <- pmax(pmin(p, 1-1e-12), 1e-12) # avoid p=0 or 1
  ll <- lchoose(total_count, alt_count) + alt_count*log(p) + (total_count-alt_count)*log(1-p)
  return(ll)
}

# dp = 10
# obs_p=0
# nalt = round(dp*obs_p)
# p_theos = seq(0.01,0.99,0.01)
# lls = sapply(p_theos, function(p_theo){loglik_binom(alt_count=nalt, total_count=dp, p=p_theo)})
# plot(p_theos, lls)

# Function that compute the log-likelihood for multiple samples and sum it. 
# alt_vector, ref_vector, expected_p should all have the same length
# We assume a little sequencing error
multisample_loglik_binom = function(alt_vector, ref_vector, expected_p, error=0.01){
  n_samples = length(alt_vector)
  ll_total = sum(unlist(sapply(1:n_samples, function(n){
    alt_count <- alt_vector[n]
    ref_count <- ref_vector[n]
    dp <- alt_count + ref_count # depth
    # sample alt frequency assuming an error (i.e., sequencing error)
    theta <- expected_p[n] 
    p <- theta*(1-error) + (1-theta)*error 
    loglik_binom(alt_count, dp, p)
  })))
  return(ll_total)
}


scenario_call <- function(REF, ALT, EXP_FREQ, error=0.01, ncores = 3){
  # REF, ALT: matrices of observed counts (rows=SNPs, cols=samples)
  # EXP_FREQ: matrix of expected alt freq per scenario (rows=scenario, cols=samples)
  # error: sequencing error
  # returns: for each SNP, best scenario and its likelihood
  
  n_snps <- nrow(REF)
  n_scenarios <- nrow(EXP_FREQ)
  
  results = do.call(rbind, mclapply(1:n_snps,mc.cores = ncores, function(i){
    if(i %% 1000 == 0){print( paste0(round(100*i/n_snps,digits = 2),'% done') )}
    
    loglik_scenario = unlist(lapply(1:n_scenarios, function(s){multisample_loglik_binom(alt_vector = ALT[i,] , ref_vector = REF[i,], expected_p = EXP_FREQ[s,], error=error)}))
    
    #loglik_scenario <- numeric(n_scenarios)
    #for(s in 1:n_scenarios){loglik_scenario[s] = multisample_loglik_binom(alt_vector = ALT[i,] , ref_vector = REF[i,], expected_p = EXP_FREQ[s,], error=error)}
    #best_s <- which.max(loglik_scenario)
    #return( data.frame(best_scenario =  best_s, max_loglik = loglik_scenario[best_s]) )
    
    return(loglik_scenario)
  }))
  
  return(results)
  
}


###

# Compute every possible genotype combination in founders
if(CHR != 'X'){
  fgeno = expand.grid(c(0,0.5,1),c(0,0.5,1),c(0,0.5,1)) # ref (0), het(0.5), alt (1), for three founders
}else{
  fgeno = expand.grid(c(0,0.5,1),c(0,0.5,1),c(0,1)) # If X chromosome, FM (male) is haploid so cannot be heterozygous
}
colnames(fgeno) = c("FA","FB", "FM")
fgeno = fgeno[!(apply(fgeno,1,sum) == 3 | apply(fgeno,1,sum) == 0 ),] # remove fixed
fgeno = as.data.frame(fgeno)

# Add expected allele frequency in pool (assuming no deviation)
# probably a wrong assumption in POP but should not matter if allelic frequency change are moderate during the 5 generations of outcrossing
fgeno$CrossApool = (fgeno$FA + fgeno$FM)/2
fgeno$CrossBpool = (fgeno$FB + fgeno$FM)/2
fgeno$CrossCpool = fgeno$FM/2

founders_possible_genotypes = fgeno[,1:3]
expected_pool_frequencies = fgeno[,-1:-3]

scenario_names = unlist(lapply(1:nrow(founders_possible_genotypes), function(i){
  paste(apply(cbind( unlist(colnames(founders_possible_genotypes)),unlist(founders_possible_genotypes[i,])),1, paste, collapse=':'), collapse=";")
}))


# Import pool allelic count
AD = fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/pools/",CHR,"_AllelicDepth_becei_CrossPools_stringentFiltering.csv.gz"))
#AD = fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/private/alignment&variantCalling/VCF/CrossPool/",CHR,"_CrossPool_variants.filtered.csv"))
AD = as.data.frame(AD)
snps = as.data.frame(AD[,1:4])
AD = as.matrix(AD[,grepl("Cross", colnames(AD))])



REF = do.call(cbind, lapply(1:ncol(AD), function(j){as.numeric(tstrsplit(AD[,j],",")[[1]])}))
ALT = do.call(cbind, lapply(1:ncol(AD), function(j){as.numeric(tstrsplit(AD[,j],",")[[2]])}))
colnames(REF) = colnames(ALT) = tstrsplit(colnames(AD),"_AD")[[1]]


expected_pool_frequencies = expected_pool_frequencies[,colnames(REF)]


founder_genotype_loglikelihood = scenario_call(REF = REF,
                     ALT= ALT,
                     EXP_FREQ = expected_pool_frequencies,
                     error=0.001)



colnames(founder_genotype_loglikelihood) = scenario_names
founder_genotype_loglikelihood = cbind(snps, as.data.frame(founder_genotype_loglikelihood))

readr::write_csv(founder_genotype_loglikelihood, file = paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/pools/",CHR,"_founder_genotype_loglikelihood.csv"))
}



#test = fread(paste0("/Users/tomparee/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RILs/suppl/pools/",CHR,"_founder_genotype_loglikelihood.csv"))
#test[9,]  
