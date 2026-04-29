
# EMMAX implementation in R (heavily commented)
# File: emmax_r_implementation.R
# Purpose: Provide a clear, well-documented implementation of the EMMAX workflow
#          for genome-wide association testing using a linear mixed model.
#
# High-level algorithm (EMMAX-style):
# 1. Compute a genomic relationship (kinship) matrix K from genotype data G.
# 2. Under the null model (no SNP effects), estimate the ratio delta = sigma_e^2 / sigma_g^2
#    by REML while treating K as known. We use an eigen-decomposition of K to make this fast.
# 3. Using the estimated delta, perform a spectral transformation of data (U' y and U' X,
#    where K = U diag(d) U'), rescale by sqrt(d + delta) so transformed residuals are i.i.d.,
#    and then do a fast ordinary least squares regression for each SNP (this is equivalent to
#    the full mixed-model GLS but much cheaper since delta was estimated once).
#
# Notes on notation used in comments and code:
#  - n: number of individuals (rows of G)
#  - m: number of markers/SNPs (columns of G)
#  - G: genotype matrix coded as 0/1/2 (counts of the reference allele). Missing allowed (NA).
#  - p_j: allele frequency for marker j (frequency of reference allele, between 0 and 1)
#  - Z: centered genotype matrix (G - 2 p_j)
#  - K: realized relationship matrix (n x n)
#  - y: phenotype vector (n x 1)
#  - X: covariate matrix (n x q), include intercept if needed
#  - K eigen-decomposition: K = U diag(d) U', where U (n x n) orthonormal, d vector of eigenvalues
#  - model: y = X b + g + e, with g ~ N(0, sigma_g^2 K) and e ~ N(0, sigma_e^2 I)
#           So Var(y) = sigma_g^2 (K + delta I) where delta = sigma_e^2 / sigma_g^2
#
# --- IMPORTANT PRACTICAL POINTS ---
#  * Many implementations of "kinship" differ only by a constant scaling factor or by whether
#    SNP columns are standardized individually. That changes numerical values and therefore
#    numerical estimates for sigma_g^2 and sigma_e^2, but downstream association p-values
#    are typically invariant to constant rescaling of K (provided delta is estimated accordingly).
#  * If you want exact numerical agreement with packages like 'sommer::A.mat()', select
#    method = "vanraden" (implemented below) which follows VanRaden's Method 1.
#  * If you want each SNP to contribute equally (i.e., standardize each SNP column so it has unit
#    variance under HWE), choose method = "per_snp".
#

# ---- Functions --------

# -----------------------
# REML: estimating delta = sigma_e^2 / sigma_g^2
# -----------------------
# Approach used:
#  - diagonalize K: K = U diag(d) U' (U orthonormal)
#  - transform y and X: y* = U' y, X* = U' X
#  - Under model, y* ~ N(X* b, sigma_g^2 (diag(d) + delta I))
#  - The (restricted) log-likelihood for delta (after profiling out b and sigma_g^2)
#    is a function we can minimize using a univariate optimizer over log(delta).
#
# Implementation notes:
#  - we optimize over log(delta) to ensure positivity and better numerical boundaries
#  - the objective used is the (profiled) REML negative log-likelihood up to a constant:
#      L(delta) = 0.5 * (sum(log(d + delta)) + (n - p) * log(RSS / (n - p)))
#    where RSS = sum( (y* - X* b_hat)' W (y* - X* b_hat) ) with W = diag(1/(d+delta)).
#  - p is the number of fixed-effect parameters (columns of X)
#  - details: the derivation is in many LMM references; this code follows the standard trick
#    used in EMMAX/GEMMA/EMMA implementations.

log_reml_fn <- function(log_delta, y_t, X_t, d){
  # This function returns the profiled REML negative log-likelihood for a given log(delta).
  #  - log_delta: real number (log scale)
  #  - y_t: U' y (vector length n)
  #  - X_t: U' X (n x p)
  #  - d: eigenvalues of K (length n)
  
  delta <- exp(log_delta)
  # weights w_i = 1 / (d_i + delta)
  w <- 1 / (d + delta)
  
  # Compute (X' W X) and solve for b_hat. This is the generalized least squares solution
  # in the transformed (diagonal) space; here we multiply each row of X_t by w (i.e. X_t * w)
  # then compute crossproducts.
  XtW <- crossprod(X_t * w, X_t)
  
  # If XtW is near-singular (e.g., covariates linearly dependent), bail out with huge value
  chol_ok <- TRUE
  ch <- try(chol(XtW), silent = TRUE)
  if(inherits(ch, "try-error")) chol_ok <- FALSE
  if(!chol_ok) return(1e30)
  
  # b_hat under weights
  b_hat <- solve(XtW, crossprod(X_t * w, y_t))
  resid <- y_t - X_t %*% b_hat
  
  # RSS in weighted space: sum(w * resid^2)
  rss <- sum(w * (resid^2))
  
  n <- length(y_t)
  p <- ncol(X_t)
  
  # profiled negative log-likelihood (no constants): 0.5 * (sum(log(d+delta)) + (n - p) * log(rss / (n - p)))
  # we return this for minimization
  ll <- 0.5 * (sum(log(d + delta)) + (n - p) * log(rss / (n - p)))
  return(ll)
}


estimate_delta_reml <- function(y, X, eig, lower = 1e-8, upper = 1e8){
  # y: phenotype vector (n)
  # X: covariate matrix (n x p)
  # eig: list with $values and $vectors for K's eigen-decomposition
  # returns an estimate of delta = sigma_e^2 / sigma_g^2 (a positive number)
  
  d <- eig$values
  U <- eig$vectors
  
  # Project into eigenbasis (U' y and U' X) - these decouple the covariance structure
  y_t <- crossprod(U, y)
  X_t <- crossprod(U, X)
  
  # Optimize over log(delta) to ensure positivity. We use a wide interval by default.
  opt <- optimize(f = function(ld) log_reml_fn(ld, y_t, X_t, d), interval = c(log(lower), log(upper)))
  delta_hat <- exp(opt$minimum)
  return(list(delta = delta_hat, opt = opt))
}




# -----------------------
# compute_kinship
# -----------------------
# Computes a genomic relationship matrix (GRM) / kinship matrix from genotype matrix G.
#  implement VanRaden Method 1: K = ZZ' / (2 sum p_j (1-p_j)),
#    where Z = G - 2 p_j (centered only). This matches sommer::A.mat() and many published uses.
#
# Arguments:
#  - G: numeric matrix n x m of genotypes coded 0,1,2 (rows individuals, cols SNPs). NA allowed.
compute_kinship_vanraden1 <- function(G) {
  G <- as.matrix(G)
  p <- colMeans(G, na.rm = TRUE) / 2
  Z <- sweep(G, 2, 2 * p, "-")
  denom <- 2 * sum(p * (1 - p), na.rm = TRUE)
  A <- tcrossprod(Z) / denom
  A
}



# -----------------------
# Main EMMAX association function
# -----------------------
# Steps implemented:
# 1) Compute or accept kinship K
# 2) Eigen-decompose K -> eig$values (d) and eig$vectors (U)
# 3) Estimate delta by REML under the null model y ~ X (no SNPs)
# 4) Transform and rescale y and X once: y_t = U' y / sqrt(d + delta), X_t = U' X / sqrt(d + delta)
# 5) For each SNP j do:
#       - impute missing genotypes in SNP j (simple mean imputation)
#       - compute g_t = U' g_j / sqrt(d + delta)
#       - run OLS regression y_t ~ X_t + g_t (no further weights). This gives beta_hat for SNP.
#    This OLS on transformed data is algebraically equivalent to doing a GLS with Var(y)=sigma_g^2(K + delta I).
#
# Notes:
#  - Because delta and the transformation are computed once, the per-SNP work is only O(n) for projection
#    and O(p^2) for solving the small linear system; this makes EMMAX very fast for GWAS.
#  - This implementation uses a simple for-loop for clarity. For large m (>100k) you should
#    vectorize or parallelize the SNP loop.
#  - p-values are computed from t-statistics with df = n - rank(X) - 1 (approx). Several implementations
#    use approximate chi-sq statistics; either is fine for large n.

emmax_assoc <- function(y, G, covar = NULL, kinship = NULL, verbose = TRUE){
  # y: vector length n
  # G: n x m genotype matrix (0/1/2)
  # covar: n x p matrix (including intercept if desired). If NULL, include intercept.
  # kinship: optional precomputed kinship matrix. If NULL, compute from G.
  y <- as.numeric(y)
  n <- length(y)
  if(is.null(covar)) covar <- matrix(1, n, 1) else covar <- as.matrix(covar)
  if(nrow(covar) != n) stop("Covariate matrix must have same number of rows as y")
  if(is.null(kinship)){
    if(verbose) message("Computing kinship from genotype matrix...")
    kinship <- compute_kinship_vanraden1(G)
  }
  # eigen decomposition
  if(verbose) message("Eigen-decomposing kinship matrix...")
  eig <- eigen(kinship, symmetric = TRUE)
  # make sure eigenvalues non-negative (numerical)
  eig$values[eig$values < 0] <- 0
  # estimate delta via REML under null (no SNP)
  if(verbose) message("Estimating variance ratio delta via REML (null model)...")
  res_delta <- estimate_delta_reml(y, covar, eig)
  delta <- res_delta$delta
  if(verbose) message(sprintf("Estimated delta = %g", delta))
  # Precompute transformed y and covariates divided by sqrt(d + delta)
  U <- eig$vectors; d <- eig$values
  y_t <- crossprod(U, y) / sqrt(d + delta)
  X_t <- crossprod(U, covar) / sqrt(d + delta)
  # Precompute part for covariates-only projection
  XtX_inv <- solve(crossprod(X_t, X_t))
  P0 <- X_t %*% XtX_inv %*% crossprod(X_t, rep(1,n)) # not used directly; we'll compute residuals per SNP
  # Loop over SNPs and test
  m <- ncol(G)
  p <- ncol(covar)
  snp_beta <- numeric(m); snp_se <- numeric(m); snp_t <- numeric(m); snp_p <- numeric(m); snp_LR <- numeric(m)
  for(j in seq_len(m)){
    gj <- G[, j]
    if(any(is.na(gj))) gj[is.na(gj)] <- mean(gj, na.rm = TRUE)
    g_t <- crossprod(U, gj) / sqrt(d + delta)
    # build design with covariates and SNP
    W <- cbind(X_t, g_t)
    XtW <- crossprod(W, W)
    # check singularity
    if(det(XtW) < 1e-12){
      snp_beta[j] <- NA; snp_se[j] <- NA; snp_t[j] <- NA; snp_p[j] <- NA
      next
    }
    
    
    b_hat <- solve(XtW, crossprod(W, y_t))
    beta_g <- b_hat[p + 1]
    # residuals and sigma2 estimate
    resid <- y_t - W %*% b_hat
    sigma2 <- sum(resid^2) / (n - ncol(W))
    var_b <- sigma2 * solve(XtW)
    se_g <- sqrt(var_b[p + 1, p + 1])
    tstat <- beta_g / se_g
    lr <- 0.5*(tstat)^2
    df <- n - ncol(W)
    pval <- 2 * pt(-abs(tstat), df)
    snp_beta[j] <- beta_g; snp_se[j] <- se_g; snp_t[j] <- tstat; snp_p[j] <- pval; snp_LR[j] <- lr
    
    # #method1
    # m0 = lm(y_t ~ X_t -1)
    # m1 = lm(y_t ~ X_t + g_t -1)
    # summary(m1)
    # anova(m0,m1,test='LRT')
    # 
    # summary(lm(y_t ~ X_t + g_t -1))
    # log_likelihood_m0 <- logLik(m0)
    # log_likelihood_m1 <- logLik(m1)
    # (log_likelihood_m1-log_likelihood_m0)
    # 
    # #method2
    # 0.5*(beta_g/se_g)^2
  }
  
  if(is.null(colnames(G))){colnames(G) = 1:ncol(G)}
  res <- data.frame(marker = colnames(G), beta = snp_beta, se = snp_se, t = snp_t, p = snp_p, LR = snp_LR)
  return(list(results = res, delta = delta, eig = eig))
}



# ----- Example simulation -----
if(FALSE){
  #set.seed(101)
  #n <- 200; m <- 5000
  # simulate genotypes (0/1/2) with MAF between 0.05 and 0.5
  #MAF <- runif(m, 0.05, 0.5)
  #G <- matrix(0, n, m)
  #for(j in seq_len(m)) G[, j] <- rbinom(n, 2, MAF[j])
  
  geno <- read_csv("~/Documents/Documents - MacBook Pro de tom/rockmanlab/becei/RIX_IBD/genotypes/geno_RIXs_pruned0.9999_F&M.csv")
  geno = as.matrix(geno)
  G = t(geno)+1
  m = nrow(geno);n=ncol(geno)
  # simulate kinship and polygenic effect
  K <- compute_kinship_vanraden1(G)
  #u <- t(chol(K + 1e-6 * diag(n))) %*% rnorm(n)
  # covariate: intercept + one covariate
  X <- cbind(1, rnorm(n))
  # pick causal SNPs
  poly <- sort(sample(1:m, 1000))
  qtl = sort(sample(1:m,3))
  beta_snps <- rep(0, m); beta_snps[poly] <- rnorm(1000,sd=0.5)*sample(c(1,-1),replace = T)
  beta_snps[qtl] = c(10,10,-20)
  apply(G[,qtl]/2,2, mean)
  g_effect <- G %*% beta_snps
  y <- X %*% c(1, 0.5) + g_effect 
  y = scale(y) + rnorm(n, sd = 0.5)
  colnames(G) <- paste0("snp", seq_len(m))
  out <- emmax_assoc(y, G, covar = X)
  head(out$results[order(out$results$p), ])
  
  res = out$results
  res$pos = 1:nrow(res)
  res$true_beta_snps = beta_snps
  res$LR <- 0.5*(res$t)^2
  
  pi = 0.01
  res$PP = (pi*exp(res$LR))/(pi*exp(res$LR) + (1-pi))
  
  ggplot(res, aes(pos, PP))+geom_point()+geom_vline(xintercept = qtl)
  ggplot(res, aes(pos, -log10(p)))+geom_point()+geom_vline(xintercept = qtl)
  ggplot(res, aes(pos, abs(beta_snps)))+geom_point()
}





# ----- Compare results to EMMAX implemented in sommer -------
# library(sommer)
# dt = data.frame(pheno=y,covar = X[,2],id = rownames(G))
# test = GWAS(pheno~covar, random=~vsr(id, Gu=K), rcov=~units, M=G, gTerm = "u:id", data=dt)
# res$pval_sommer = test$pval
# ggplot(res, aes(-log10(pval_sommer), -log10(p)))+geom_point()




# whiten_pheno <- function(y, VCOV) {
#   # Input checks
#   if (!is.numeric(y)) stop("y must be numeric")
#   if (!is.matrix(VCOV)) stop("VCOV must be a matrix")
#   if (length(y) != nrow(VCOV)) stop("length of y must match dimensions of VCOV")
#   
#   # Try Cholesky first
#   chol_success <- TRUE
#   L <- tryCatch(chol(VCOV), error = function(e) { chol_success <<- FALSE; NULL })
#   
#   if (chol_success) {
#     # Solve L %*% z = y  →  z = L^{-1} y
#     z <- backsolve(L, y, transpose = FALSE)  # fast triangular solve
#   } else {
#     # Fallback to eigen decomposition
#     eig <- eigen(VCOV, symmetric = TRUE)
#     vals <- eig$values
#     vecs <- eig$vectors
#     
#     # Guard against tiny/negative eigenvalues
#     vals[vals < 1e-10] <- 1e-10
#     
#     # Compute V^{-1/2} = Q * diag(1/sqrt(vals)) * Q^T
#     V_inv_half <- vecs %*% diag(1 / sqrt(vals)) %*% t(vecs)
#     
#     z <- V_inv_half %*% y
#   }
#   
#   return(as.numeric(z))
# }
# 


# -----------------------
# Small utilities & examples
# -----------------------
# Function to compute an approximate narrow-sense heritability (on the data scale) from delta.
# Caution: This is only correct if the K matrix is scaled such that mean(diag(K)) == 1 and
# the model variance decomposition is Var(y) = sigma_g^2 K + sigma_e^2 I. Then:
# h2 = sigma_g^2 / (sigma_g^2 + sigma_e^2) = 1 / (1 + delta)
# (because delta = sigma_e^2 / sigma_g^2). If you don't scale K this way, h2 will be wrong.
#h2_from_delta <- function(delta){
#  return(1 / (1 + delta))
#}



# denoise_wiener <- function(y, V, eps = 1e-10) {
#   # y: numeric vector (length n)
#   # V: n x n symmetric positive semidefinite covariance matrix (measurement error)
#   # eps: small constant to stabilize division when needed
#   
#   # 1. Center y
#   y <- as.numeric(y - mean(y))
#   n <- length(y)
#   
#   # 2. Eigen-decompose noise covariance
#   eig <- eigen(V, symmetric = TRUE)
#   U <- eig$vectors
#   lambda <- pmax(eig$values, 0) # prevent negative eigenvalues due to numerics
#   
#   # 3. Estimate signal variance (scalar)
#   total_var <- sum(y * y) / n
#   noise_var <- sum(lambda) / n
#   sigma_s2 <- max(0, total_var - noise_var)
#   sigma_s2 <- max(sigma_s2, eps)
#   
#   # 4. Compute Wiener weights
#   w <- sigma_s2 / (sigma_s2 + lambda)
#   
#   # 5. Transform y
#   y_t <- U %*% (w * (t(U) %*% y))
#   
#   return(as.numeric(y_t))
# }
