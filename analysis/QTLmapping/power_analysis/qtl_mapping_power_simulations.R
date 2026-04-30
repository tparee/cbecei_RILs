# Load necessary library
setwd('/Users/Sol/Documents/power_analysis_gwas/')
library(MASS)
library(sommer)
library(lme4)
library(data.table)
require("psych")
source("EMMAX_functions.R")
library(parallel)

nc = parallel::detectCores(logical=T) -1

IDOUT = 'panelA&B'

get.K_ASV = function(GT){
  ASV = scale(GT,center=T,scale=F) %*% t(scale(GT,center=T,scale=F))
  ASV = ASV / (psych::tr(ASV)/(nrow(GT)-1))
  ASV
}

calculate.genotype.value = function(gt, effects){
  c(t(gt) %*% effects)
}

calculate.scaling.factor = function(gt, effects, aimedVariance){
  ratioSD = sqrt(aimedVariance)/sd(calculate.genotype.value(gt=gt, effects=effects))
  ratioSD
}


simulate_phenotypes = function(qtl_variance_ratio, QTLpos, genotypes, heritability=0.5, phenotypic_variance = 1, n_loci = 1000){
  # Determine variance components
  qtl_variance <- phenotypic_variance * qtl_variance_ratio
  polygenic_variance <- (phenotypic_variance * heritability) - qtl_variance
  residual_variance <- phenotypic_variance * (1-heritability)
  polygenic_background = sort(sample(1:nrow(genotypes), n_loci))
  
  # Simulate effects for polygenic loci
  polygenic_effects <- rnorm(length(polygenic_background))
  polygenic_effects = polygenic_effects * calculate.scaling.factor(gt=genotypes[polygenic_background,],
                                                                   effects=polygenic_effects,
                                                                   aimedVariance=polygenic_variance)
  qtl_effect = calculate.scaling.factor(gt=matrix(genotypes[QTLpos,], nrow=1),
                                        effects=1,
                                        aimedVariance=qtl_variance)
  
  
  
  # Calculate genotypic value
  genotypic_value <- calculate.genotype.value(gt=genotypes[polygenic_background,], effects = polygenic_effects) + genotypes[QTLpos,]*qtl_effect
  
  # Simulate residuals
  residuals <- rnorm(n_individuals, mean = 0, sd = sqrt(residual_variance))
  
  # Calculate phenotypic value
  phenotype <- genotypic_value+residuals
  #phenotype = rnorm(n_individuals)
  phenotype  = data.frame(id = colnames(genotypes), value = phenotype)
  return(phenotype)
}

mapping = function(phenotype, GT, K, perm=F){
  modh2 <- mmer(
    value ~ 1,
    random = ~ vsr(id, Gu = K),
    data = phenotype
  )
  vcomp = summary(modh2)$varcomp[1:2,1]
  delta = vcomp[2]/vcomp[1]
  if(perm){ix = sample(1:nrow(GT))}else{ix = 1:nrow(GT)}
  res = emmax_assoc(y=phenotype$value,
                    G=GT[ix,], covar = NULL, kinship = K,
                    delta = delta, verbose = F)
  return(res)
}



gwasstats = function(gwas, QTLpos, threshold){
 
  sign.snps = which(gwas$p < threshold)
  isSignificant = QTLpos %in% sign.snps
  
  if(isSignificant){
    
    highestHit = which.min(gwas$pval)
    
    distance_cM = abs(diff(gwas$cm[c(highestHit, QTLpos)]))
    distance_bp = abs(diff(gwas$pos[c(highestHit, QTLpos)]))
    
    QTLsize_cM=diff(range(gwas$cm[sign.snps]))
    QTLsize_bp=diff(range(gwas$pos[sign.snps]))
    
    stats = data.frame(QTLpos,
                       isSignificant,
                       distance_cM,
                       distance_bp,
                       QTLsize_cM,
                       QTLsize_bp)
  }else{
    
    stats = data.frame(QTLpos,isSignificant,
                       distance_cM=NA,
                       distance_bp=NA,
                       QTLsize_cM=NA,
                       QTLsize_bp=NA)
  }
  
  stats
}



# Set simulation parameters
nsimulation = 500 # number of run per parameter
n_loci <- 1000       # Number of polygenic loci
heritability <- 0.5           # Total heritability
QTL_variance_ratio <- c(0.05, 0.9, 0.12, 0.15, 0.18, 0.25)   # Proportion of phenotypic variance explained by QTL
phenotypic_variance <- 1 # Phenotypic variance


# Import genotypes (pruned, LD > 0.9999)
genotypes =  as.matrix(fread("beceiPanels_geno_RILs_pruned0.999.csv.gz"))
genotypes = (genotypes-0.5)*2
snps =  as.data.frame(fread("beceiPanels_variantsInfo_pruned0.999.csv.gz"))
af = apply(genotypes,1,function(x){sum(x+1, na.rm=T)/(2*sum(!is.na(x)))})
genotypes = genotypes[af > 0 & af < 1,]
snps = snps[af > 0 & af < 1,]
n_individuals = ncol(genotypes)
GT = t(genotypes)

#growthrates$id = factor(growthrates$strain, levels = rownames(GT))
K=get.K_ASV(GT)


for(qtl_variance_ratio in QTL_variance_ratio){
  print
  # Permutations
  nullpvals = unlist(parallel::mclapply(1:500, mc.cores = nc, function(p){
    if(p %% 50 == 0){print(p)}
    QTLpos = sample(1:nrow(genotypes),1)
    pheno = simulate_phenotypes(qtl_variance_ratio = qtl_variance_ratio,
                                     QTLpos = QTLpos,
                                     genotypes = genotypes,
                                     heritability = heritability,
                                     phenotypic_variance = phenotypic_variance,
                                     n_loci = n_loci)
    pres = mapping(phenotype = pheno, GT=GT, K, perm=T)
    return(min(pres$pval, na.rm = T))
  }))
  
  save(nullpvals, file = paste0("permutated_minp_",IDOUT, "_", qtl_variance_ratio, ".Rdata"))
  thresholds = data.frame(th =  c(0.01, 0.05, 0.1), p=quantile(nullpvals, probs = c(0.01, 0.05, 0.1)))
  
  
  powersim = do.call(rbind, parallel::mclapply(1:nsimulation, mc.cores = nc, function(n){
    if(n %% 50 == 0){print(n)}
    QTLpos = sample(1:nrow(genotypes),1)
    pheno = simulate_phenotypes(qtl_variance_ratio = qtl_variance_ratio,
                                     QTLpos = QTLpos,
                                     genotypes = genotypes,
                                     heritability = heritability,
                                     phenotypic_variance = phenotypic_variance,
                                     n_loci = n_loci)
    res = mapping(phenotype = pheno, GT=GT, K=K, perm=F)
    res = cbind(snps, res)
    
    #ggplot(res, aes(pos, -log10(pval)))+
    #      geom_point()+facet_wrap(~chrom)+
     #       geom_point(data=res[QTLpos,], aes(pos, -log10(pval)), color='red')+
     #       geom_hline(yintercept = -log10(threshold))+
      #      theme(legend.position = "none")

    statout = do.call(rbind, lapply(1:nrow(thresholds), function(i){
      x = gwasstats(gwas=res, QTLpos=QTLpos, threshold=thresholds$p[i])
      x$alpha = thresholds$th[i]
      x$qtl_variance_ratio = qtl_variance_ratio
      x$nrun = n
      x
    }))
    return(statout)
  }))
  
  
  save(powersim, file = paste0("power_",IDOUT, "_", qtl_variance_ratio, ".Rdata"))
  
}



