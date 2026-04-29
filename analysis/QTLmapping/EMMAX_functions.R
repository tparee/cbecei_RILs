
# EMMAX implementation in R
# EMMAX is LMM for GWAS method that improves computational efficiency  
# by estimating variance components once and reusing them across all marker tests.
# Note: This in a approximation because variance components are not reestimated at each marker and is not properly conditionned
# STEPS:
# 0. Pre-compute a genomic relationship (kinship) matrix K from genotype data G 
#    and (no SNP effects), estimate the ratio delta = sigma_e^2 / sigma_g^2. (not this function)
# 1. Spectral transformation: eigen-decompose K and rotate/scale the data in the eigenbasis
#    converting the mixed model with correlated errors into an equivalent model with independent errors
# 2. And then do a fast OLS regression for each marker


# 1) Eigen-decompose K -> eig$values (d) and eig$vectors (U)
# 3) Estimate delta by REML under the null model y ~ X (no SNPs)
# 4) Transform and rescale y and X once: y_t = U' y / sqrt(d + delta), X_t = U' X / sqrt(d + delta)
# 5) For each SNP j do:
#       - impute missing genotypes in SNP j (simple mean imputation)
#       - compute g_t = U' g_j / sqrt(d + delta)
#       - run OLS regression y_t ~ X_t + g_t (no further weights). This gives beta_hat for SNP.
#    This OLS on transformed data is algebraically equivalent to doing a GLS with Var(y)=sigma_g^2(K + delta I).


emmax_assoc <- function(y, G, covar = NULL, kinship, delta, verbose = TRUE){
  # y: vector of phenotypic values of length n
  # G: n x m genotype matrix
  #(!) y and G must match. i.e., genotype G[i,] and phenotype y[i] belong to individual i
  # covar: n x p matrix (including intercept if desired); can be obtained with (model.matrix(~a+b,data)). If NULL, include intercept.
  # kinship: covariance matrix
  # delta: sigma_e^2 / sigma_g^2 (environmental variance / genetic variance)
 
  y <- as.numeric(y)
  n <- length(y)
  if(is.null(covar)) covar <- matrix(1, n, 1) else covar <- as.matrix(covar)
  if(nrow(covar) != n) stop("Covariate matrix (covar) must have same number of rows as y")
  if(nrow(covar) != n) stop("Genotype matrix (G) must have same number of rows as y")
  
  #1) Eigen-decompose K -> eig$values (d) and eig$vectors (U)
  if(verbose) message("Eigen-decomposing kinship matrix...")
  eig <- eigen(kinship, symmetric = TRUE)
  eig$values[eig$values < 0] <- 0 # make sure eigenvalues non-negative (numerical)
  U <- eig$vectors; d <- eig$values
  
  #2) Spectral Transformation: decorrelate data and remove kinship-induced covariance
  # => Precompute transformed y and covariates divided by sqrt(d + delta)
  y_t <- crossprod(U, y) / sqrt(d + delta) #Transformed phenotypes: y_t = U' y / sqrt(d + delta)
  X_t <- crossprod(U, covar) / sqrt(d + delta) # Transformed covariates: X_t = U' X / sqrt(d + delta)
  
  if(verbose) message("OLS completed by:")
  m <- ncol(G) # number of markers
  results = do.call(rbind, lapply(1:m, function(j){
    if(j %% 1000 == 0 & verbose){message( paste0(round(100*j/m, digits = 2), "%") )}
    
    gj <- G[, j] # genotype vector at the marker j
    # Mean impute if some missing
    if(any(is.na(gj))) gj[is.na(gj)] <- mean(gj, na.rm = TRUE)
    g_t <- crossprod(U, gj) / sqrt(d + delta) # spectral transformation
    
    # Ordinary least square regression on tranformed data
    ols = lm(y_t ~ X_t + g_t -1)
    p = summary(ols)$coef["g_t",4]
    beta = summary(ols)$coef["g_t",1]
    se = summary(ols)$coef["g_t",2]
    return(data.frame(pval = p, beta = beta, SE = se))
  }))
    
  return(results)
}


