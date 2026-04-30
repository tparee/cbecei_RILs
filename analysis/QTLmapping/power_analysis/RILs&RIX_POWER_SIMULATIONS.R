# Load necessary library
library(MASS)
library(sommer)
library(lme4)

calculate.genotype.value = function(gt, effects){
 
  c(t(gt) %*% effects)
}

calculate.scaling.factor = function(gt, effects, aimedVariance){
  ratioSD = sqrt(aimedVariance)/sd(calculate.genotype.value(gt=gt, effects=effects))
  ratioSD
}




# Set simulation parameters
#n_individuals <- 200         # Number of individuals
nsimulation = 500
n_loci <- 1000       # Number of polygenic loci
heritability <- 0.5           # Total heritability
QTL_variance_ratio <- c(0.05, 0.9, 0.12, 0.15, 0.18, 0.25)   # Proportion of phenotypic variance explained by QTL

# Phenotypic variance
phenotypic_variance <- 1


# Simulate genotypes for polygenic loci (-1, 0, 1 coding)
#genotypes <- matrix(sample(c(-1, 0, 1,0.5), 
#                                     n_individuals * n_loci, 
#                                     replace = TRUE), 
#                              ncol = n_individuals, 
#                              nrow = n_loci)

#colnames(genotypes) = paste0("L",1:ncol(genotypes))


# Import genotypes (pruned, LD > 0.9999)
genotypes =  as.matrix(read_csv("~/Downloads/beceiPanels_geno_RILs_pruned0.9999.csv"))
#genotypes = genotypes[,grepl("A_",colnames(genotypes))]
afs = apply(genotypes, 1, function(x){mean(x[x==0 | x==1],na.rm=T)})
snps =  as.data.frame(read_csv("~/Downloads/beceiPanels_snps_pruned0.9999.csv"))
ix = (afs > 0.02 & afs < 0.98)
genotypes = genotypes[ix,]
snps = snps[ix,]
genotypes = (genotypes-0.5)*2
genotypes[is.na(genotypes)]=0
n_individuals = ncol(genotypes)

K = A.mat(t(genotypes))

# Compute the th
minp = lapply(1:nsimulation, function(i){
  
  phenotype = rnorm(n_individuals)
  phenotype  = data.frame(id = colnames(genotypes), idd=colnames(genotypes), value = phenotype)
  
  gwas <- GWAS(value~1,
               random=~vsr(id, Gu=K),
               rcov=~units, M=t(genotypes), gTerm = "u:id",
               data=phenotype)
  
  min(gwas$pvals)
})

minp=unlist(minp)

thresholds = data.frame(th =  c(0.01, 0.05, 0.1), p=quantile(minp, probs = c(0.01, 0.05, 0.1)))


simu_results = NULL
for(qtl_variance_ratio in QTL_variance_ratio){
  
  # Determine variance components
  qtl_variance <- phenotypic_variance * qtl_variance_ratio
  polygenic_variance <- phenotypic_variance * heritability - qtl_variance
  residual_variance <- phenotypic_variance - heritability
  
  res = do.call(rbind, lapply(1:nsimulation, function(ID){
    
    
    # Sample the position of the loci coding the trait
    trait_loci = sample(1:nrow(genotypes), n_loci+1)
    polygenic_background = sort(trait_loci[-1])
    QTLpos = trait_loci[1]
    
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
    phenotype  = data.frame(id = colnames(genotypes), idd=colnames(genotypes), value = phenotype)
    phenotype$cross = ifelse(grepl("A_",phenotype$id), 'A', 'B')
    ##########################################
    ##########################################
    
    
    #K = A.mat(t(genotypes))
    #D = A.mat(t(genotypes))
    
    #mix = mmer(value~1,
    #          random = ~ vsr(id, Gu=K),
    #           rcov = ~units,
    #          data= phenotype)
    
    #mix2 = mmer(value~1,
    #           random = ~ vsr(id, Gu=K) + vsr(idd, Gu=D),
    #           rcov = ~units,
    #           data= phenotype)
    
    
    #summary(mix)$varcomp
    #vpredict(mix, h2 ~ (V1) /  (V1+V2))
    #vpredict(mix2, H2 ~ (V1+V2) /  (V1+V2+V3) )
    
    #gwas <- GWAS(value~cross,
    #             random=~vsr(id, Gu=K),
    #             rcov=~units, M=t(genotypes), gTerm = "u:id",
    #             data=phenotype)
    
    gwas = cbind(snps, p= gwas$pvals)
    #gwas$p[gwas$p==0]=NA
    
    

    #ggplot(gwas, aes(pos, -log10(p)))+
    #  geom_point()+facet_wrap(~chrom)+
    #  geom_point(data=gwas[QTLpos,], aes(pos, -log10(p)), color='red')+
    #  geom_hline(data=thresholds, aes(yintercept = -log10(p), color=as.factor(th)))+
    #  theme(legend.position = "none")
    
    
    
    do.call(rbind,lapply(1:nrow(thresholds), function(i){
      
      th = thresholds[i,]$p
      sign.snps = which(gwas$p < th)
      
      isSignificant = QTLpos %in% sign.snps
      
      if(isSignificant){
        
        highestHit = which.min(gwas$p)
        
        distance_cM = abs(diff(gwas$cm[c(highestHit, QTLpos)]))
        distance_bp = abs(diff(gwas$pos[c(highestHit, QTLpos)]))
        
        QTLsize_cM=diff(range(gwas$cm[sign.snps]))
        QTLsize_bp=diff(range(gwas$pos[sign.snps]))
        
        stats = data.frame(simID=ID, qtl_variance_ratio, heritability,
                           QTLpos,threshold = thresholds[i,]$th,isSignificant,
                           distance_cM,
                           distance_bp,
                           QTLsize_cM,
                           QTLsize_bp)
      }else{
        
        stats = data.frame(simID=ID,  qtl_variance_ratio, heritability,
                           QTLpos,threshold = thresholds[i,]$th, isSignificant,
                           distance_cM=NA,
                           distance_bp=NA,
                           QTLsize_cM=NA,
                           QTLsize_bp=NA)
      }
      
      stats
      
      
    })) 
    
    
  }))
  
  simu_results = rbind(simu_results, res)
  
}


