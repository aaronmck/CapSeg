## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


DmvNorm <- function (x, mean, sigma, log=FALSE) {
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    
    if (missing(mean)) {
        mean <- rep(0, length = ncol(x))
    }

    if (missing(sigma)) {
        sigma <- diag(ncol(x))
    }

    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }

    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
        stop("sigma must be a symmetric matrix")
    }
    
    if (length(mean) != NROW(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    
    dist.val <- mahalanobis(x, mean, sigma)
    log.det <- sum(log(eigen(sigma, TRUE,
                            only.values=TRUE)[["values"]]))
    logretval <- -(ncol(x) * log(2*pi) + log.det + dist.val) / 2

    if (!log) {
      log.retval <- exp(log.retval)
    }
    
    return(log.retval)
}
  

FastDmvt <- function(x, delta, sigma, invert.sigma=NULL, df=1, log=TRUE) {
#  if (is.null(invert.sigma)) {
#    invert.sigma <- solve(sigma)
#  }
  
  if (df == 0) {
    return(DmvNorm(x, delta, sigma, log=log))
  }

  sigma.ncol <- NCOL(sigma)

 # dist.val <- mahalanobis(x, delta, invert.sigma, inverted=TRUE)
  dist.val <- mahalanobis(x, delta, sigma)
  log.det <- sum(log(eigen(sigma, TRUE, only.values=TRUE)[["values"]]))
  log.retval <- lgamma((sigma.ncol + df) / 2) - (lgamma(df / 2) + 0.5 *
                                                 (log.det + sigma.ncol * logb(pi * df))) -
                                                   0.5 * (df + sigma.ncol) *
                                                     logb(1 + dist.val / df)
  if (!log) {
    log.retval <- exp(log.retval)
  }
  
  return(log.retval)
}


DScaledT <-  function(x, mu, sigma, nu, log=FALSE) {
   log.norm <- log(gamma((nu + 1) /2 ) / (gamma(nu / 2) * sqrt(nu * pi) * sigma))
   log.d <- log(1 + ((x - mu) / sigma)^2 * 1 / nu ) * - (nu + 1) / 2 + log.norm

   if (!log) {
     out <- exp(log.d)
   } else {
     out <- log.d
   }

   return(out)
}

SInvChisqMl <- function(x) {
  res <- GammaMl(1 / x)
  nu <- res[["alpha"]] * 2
  sigma <- sqrt(res[["beta"]] / res[["alpha"]])

  return(list(nu.hat=nu, sigma.hat=sigma))
}

GammaMl <- function(x) {
  obj <- function(k) {
    res <- ((log(k) - digamma(k)) - s)^2
    deriv <- 2*(-s + log(k) - digamma(k)) * (1/k - trigamma(k))
    attr(res, "gradient") <- deriv
    return(res)
  }

  x.len <- length(x)
  s <- log(1/x.len * sum(x)) - 1/x.len * sum(log(x))

  k.init <- 1
  res <- nlm(f=obj, p=k.init)
  k.hat <- res[["estimate"]]
  theta.hat <- sum(x) / (k.hat * x.len)

  return(list(alpha=k.hat, beta=1/theta.hat))
}

RSInvChisq <- function(n, nu, sigma) {
  return(sqrt((sigma^2 * nu) / rchisq(n, nu)))
}
