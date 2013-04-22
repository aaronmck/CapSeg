## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetPermStirlingCoefs <- function(W, verbose=FALSE) {
  N <- length(W)
  sum.coefs <- rep(0, N)
  last.coefs <- rep(1, N)
  epsilon <- 5e-04
  C <- Inf
  i <- 1
  
  while (C > epsilon) {
    conf.w <- sample(W, length(W))
    new.coefs <- GetStirlingCoefs(cumsum(conf.w))
    sum.coefs <- sum.coefs + new.coefs
    coefs <- sum.coefs / i 
    i <- i + 1
    C <- abs(1 - sum((coefs / last.coefs), na.rm = TRUE) / N)  
    last.coefs <- coefs
    if (verbose) {
      cat(".")
    }
    if ((i %% 25) == 0) {
      if (verbose) {
        print(C)
      }
    }
  }
  if (verbose) {
    cat("\n")
  }
  return(coefs)
}

GetK0Map <- function(coefs) {
    N <- length(coefs)    
    gamma.grid <- exp(seq(log(5e-04), log(0.03), length.out = 100))
    n.gamma <- length(gamma.grid)
    k.0.map <- matrix(0, ncol = 2, nrow = n.gamma)
    colnames(k.0.map) <- c("k_0", "gamma")
    k.0.map[, 2] <- gamma.grid
    
    dp.prob <- array(0, dim = c(n.gamma, N))
    for (i in 1:n.gamma) {
      for (k in 1:N) {
        dp.prob[i, k] <- exp(CalcDpLoglik(N, k, gamma.grid[i], coefs))
      }
      k.0.map[i, 1] <- sum(dp.prob[i, ] * c(1:N))
    }
    
    return(k.0.map)
}

GetStirlingCoefs <- function(W) {
  ## if W = c(1:N), then this function returns the (unisgned 1st kind of)
  ## Stirling numbers for N, k=c(1:N)
  ## starts to give incorrect results at N = 19 if fft is used
  my.conv <- function(x, y) {
    ## y is len 2
    x <- c(0, x, 0)
    
    N <- length(x) - 2
        res <- rep(0, N + 1)
    
    for (k in 1:(N + 1)) {
      res[k] <- x[k] * y[1] + x[k + 1] * y[2]
    }
    return(res)
  }
  
  N <- length(W)
  n.w <- W
  cres <- 1
  for (i in 1:(N - 1)) {
    cres <- my.conv(cres, rev(c(1, n.w[i])))
  }
  
  return(rev(cres))
}
