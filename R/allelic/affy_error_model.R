## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

AffyInitCnBg = function(snp.d, cn.d) {
	# Quickly (and approximately) finds the bg of the cn platform relative to the snp platform.
	n.seg = length(snp.d)
	medians = foreach (i = 1:n.seg, .combine=rbind) %dopar% {
		c(median(colSums(snp.d[[i]])), median(cn.d[[i]]))
	}
	
	snp.medians = medians[,1]
	cn.medians = medians[,2]
	m.snp = median(snp.medians[complete.cases(snp.medians)])
	m.cn = median(cn.medians[complete.cases(cn.medians)])
	return(m.cn - m.snp)
}

AffyInitTheta <- function(h.seg.dat, verbose=FALSE) {

		kHetCovCoefs <- c(-0.05, 1.0 )
		
		## An 'at' value of 0 is equivalent to no attenuation (e.g. WGS)
		theta <- list(sigma.epsilon=0.18, sigma.eta= 0.18, nu=3.5,
				at=0.1, alpha=0, bg=0)
		
		sigma.h <- GetSigmaH(theta[["sigma.epsilon"]], theta[["sigma.eta"]])       
		tmp.hc <- (kHetCovCoefs[1] + theta[["sigma.eta"]] * kHetCovCoefs[2])^2
		tmp.hc <- min(sigma.h^2 - 0.001, tmp.hc)
		theta[["het.cov"]] <- tmp.hc
		## outlier model is uniform on {0-5, 0-5}
		theta[["p.snp.cond.out"]] <- 1 / 25
		
		theta[["cn.bg"]] <- AffyInitCnBg(h.seg.dat[["h.snp.d"]], h.seg.dat[["h.cn.d"]])
		theta[["sigma.eta.scale"]] = 1
		# theta[["sigma.scale.capseg"]] = .25
		# theta[["at.capseg"]]=0.1
		# theta[["nu.capseg"]] = 3.5
		
		return(theta)
}

AffyCalcSnpLogLik <- function(d, delta.tau, out.p, snp.gt.p, theta) {
  return(DoCalcSnpLogLik(d, delta.tau, out.p, snp.gt.p, theta))
}

AffyCalcCnLogLik <- function(d, delta.tau, out.p, theta) { 

  if (length(d) == 0) return(0)
	nu = theta[["nu"]]
	sigma.epsilon <- theta[["sigma.epsilon"]]
	sigma.eta.scale = theta[["sigma.eta.scale"]]
	sigma.eta <- theta[["sigma.eta"]] * sigma.eta.scale
	d <- HTx((d - theta[["cn.bg"]]), sigma.epsilon, sigma.eta)
	jac = HTxVar(d, sigma.epsilon, sigma.eta)
	
	varstab.delta.tau <- HTx(delta.tau, sigma.epsilon, sigma.eta)
	sigma.h = GetSigmaH(sigma.epsilon, sigma.eta)
	
	lik =  log(1-out.p) + d_scaled_t(d, varstab.delta.tau[2], sigma.h, nu, log=TRUE) + log(jac) 
	outlier = rep(log(out.p * (1 / 5) ), length(lik))
	mat = cbind(lik[1,], unlist(outlier))
	
	ll = LogAdd(mat)
	return(ll)
}

AffyGetLL <- function(par, ix, h.seg.dat, out.p, theta) {
  # browser()
	## Platform specific likelihood function for the affy
	## (and arrays in general)

  snp.d <- h.seg.dat[["h.snp.d"]][[ix]]
  cn.d <- h.seg.dat[["h.cn.d"]][[ix]]
  snp.gt.p <- h.seg.dat[["h.snp.gt.p"]][[ix]]

	dist <- abs(par[1])
	t <- par[2]
	
	if ((dist < 0) || (dist > t)) {
		return (Inf)
	}
	delta.tau <- c(dist, t)
	
	log.liks = c( sum(AffyCalcSnpLogLik(snp.d, delta.tau, out.p, snp.gt.p, theta)), 
				sum(AffyCalcCnLogLik(cn.d, delta.tau, out.p, theta)) )
	LL = sum(log.liks[complete.cases(log.liks)])	
	
	
	return(ifelse(is.nan(LL), Inf, -LL))
}

AffyGetLLDEP <- function(par, d, out.p, snp.gt.p, theta) {
  
  ## Platform specific likelihood function for the affy
  ## (and arrays in general)
  dist <- abs(par[1])
  t <- par[2]

  if ((dist < 0) || (dist > t)) {
    return (Inf)
  }

  use.e.mu <- AffyGetMeans(dist, t)
  LL <- sum(AffyCalcSnpLogLik(d, use.e.mu, out.p, snp.gt.p, theta))

  return(ifelse(is.nan(LL), Inf, -LL))
}

AffyGetMeans <- function(d, t) {
  return(c(((t / 2) - (d / 2)), ((t / 2) + (d / 2)), t))
}

AffyGetTau = function(e.mu) {
	# e.mu is a 3 column matrix with t/2 - d/2, t/2 - d/2, and t
	return(e.mu[,3])
}

AffyGetDelta = function(e.mu) {
	return(AffyGetTau(e.mu) - 2 * e.mu[,1])
}

AffyAtten <- function(r, at) {
  return(r * (1 + at) / (1 + at * r))
}

AffyInvAtten <- function(r, at) {
	r / (1 + at - (at * r))
}

AffyDmvFunc <- function(x, mu, sigma, invert.sigma, nu, log=FALSE) {
  return(FastDmvt(x, mu, sigma, invert.sigma=invert.sigma, df=nu, log=log))
}

AffyGetSnpClustLik <- function(d, delta.tau, theta) {
  d <- t(d)
  delta = delta.tau[1]
  tau = delta.tau[2]
  sigma.epsilon <- theta[["sigma.epsilon"]]
  sigma.eta <- theta[["sigma.eta"]]
  d <- HTx((d - theta[["bg"]]), sigma.epsilon, sigma.eta)
  d1.jac <-  HTxVar(d[, 1], sigma.epsilon, sigma.eta)
  d2.jac <-  HTxVar(d[, 2], sigma.epsilon, sigma.eta)
  het.cov <- theta[["het.cov"]]
  unatten.mu1 = (tau - delta) / 2
  unatten.mu2 = (tau + delta) / 2
  unatten.mu3 = tau 
  het.cov <- min(het.cov, (het.cov * unatten.mu1 * unatten.mu2))
  mu <- HTx(AffyAtten(c(unatten.mu1, unatten.mu2, unatten.mu3), theta[["at"]]), sigma.epsilon, sigma.eta)
  names(mu) <- c("atten.mu1", "atten.mu2", "atten.mu3")

  sigma.h <- GetSigmaH(sigma.epsilon, sigma.eta)
  sigma <- matrix(c(sigma.h^2, het.cov, het.cov, sigma.h^2), nrow=2, ncol=2)
  hom.sigma <- matrix(c(sigma.h^2, 0, 0, sigma.h^2), nrow=2, ncol=2)

  out = AffyBuildSnpClustLik(d, mu, sigma, hom.sigma, d1.jac, d2.jac, theta[["nu"]], theta[["p.snp.cond.out"]])
  rownames(out) <- rownames(d)
  return(out)
}

AffyBuildSnpClustLik <- function(d, e.mu, sigma, hom.sigma, d1.jac, d2.jac,
                                 nu, p.snp.cond.out) {
  PartialDmvFunc <- function(mu, sigma, invert.sigma=NULL) {
    return(AffyDmvFunc(d, mu, sigma, invert.sigma, nu))
  }

#  inverted.sigma <- solve(sigma)
#  inverted.hom.sigma <- solve(hom.sigma)
inverted.sigma <- NULL
  inverted.hom.sigma <- NULL
  
  ## FIXME: Can the 4 DmvFuncs be foreach()'d and combined via cbind()?
  clust.lik <- array(0, dim=c(nrow(d), 5))
  clust.lik[, 1] <- PartialDmvFunc(c(0, e.mu[3]), hom.sigma,
                            inverted.hom.sigma)
  clust.lik[, 2] <- PartialDmvFunc(c(e.mu[2], e.mu[1]), sigma,
                            inverted.sigma)
  clust.lik[, 3] <- PartialDmvFunc(c(e.mu[3], 0), hom.sigma,
                            inverted.hom.sigma)
  clust.lik[, 4] <- PartialDmvFunc(c(e.mu[1], e.mu[2]), sigma,
                            inverted.sigma)
  clust.lik[, 1:4] <- clust.lik[, 1:4] * d1.jac * d2.jac
  clust.lik[, 5] <- p.snp.cond.out
  
  return(clust.lik)
}

AffyTauPlatformSpecificInitialization <- function(h.d, seg.info) {
  ## NOTE: seg.info is placeholder
  
  ## platform specific initialization for optim.
  ## This function returns a three part data frame, with the 
  ## following values for tau:
  ## ( 0.25 median of the column sums, 0.75 median, 1.0 * median)  
  h.mu.t <- sapply(sapply(h.d, colSums), median)
  h.mu.t[h.mu.t < 0] <- 0
  
  i.mu.t <- h.mu.t
  half <- i.mu.t / 2
  quart <- i.mu.t / 4

  return(cbind(quart, (i.mu.t - quart), i.mu.t))
}

AffyPlatformSpecificOptimization <- function(idx, delta.tau, h.seg.dat, out.p, theta, verbose=FALSE) {
	## platform specific optimization - the optimization
	## we do for the affy (and arrays in general)
	
	return(GridstartH1OptMeans(idx, h.seg.dat, delta.tau, out.p, theta, verbose=verbose))
}

AffyThetaOpt <- function(h.seg.dat, delta.tau, out.p, theta, verbose=FALSE) {
	if (verbose) print("Optimizing Affy Theta")	
	theta[["at"]] <- HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "at", list("lower"= 0.0, "upper"=0.2), "%", 1e-3, probe.types=c("snp"), verbose=verbose)
	theta[["sigma.eta"]] <- HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "sigma.eta", list("lower"= 0.05, "upper"=1.0), ".", 1e-3, probe.types=c("snp", "cn"), verbose=verbose)
	theta[["sigma.epsilon"]] <- HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "sigma.epsilon", list("lower"= 0.001, "upper"=0.5), "*", 1e-3, probe.types=c("snp", "cn"), verbose=verbose)
	
	max.het.cov <- GetMaxHetCov(theta)
	if (max.het.cov > 0) {
		theta[["het.cov"]] = HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "het.cov", list("lower"= 0, "upper"=max.het.cov), "~", 1e-4, probe.types=c("snp"), verbose=verbose)
	} else {
		theta[["het.cov"]] <- 0
	}
	
	theta[["nu"]] <- HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "nu", list("lower"= 2, "upper"=15), "&", 1e-1, probe.types=c("snp"), verbose=verbose)
	theta[["bg"]] <- HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "bg", list("lower"= -0.1, "upper"=0.1), "@", 1e-4, probe.types=c("snp"), verbose=verbose)
	theta[["cn.bg"]] <- HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "cn.bg", list("lower"= -0.1, "upper"=1), "()", 1e-4, probe.types=c("cn"), verbose=verbose)
	theta[["sigma.eta.scale"]] <- HThetaOpt(h.seg.dat, delta.tau, out.p, theta, "sigma.eta.scale", list("lower"= 0, "upper"=30), "$", 1e-4, probe.types=c("cn"), verbose=verbose)
	
	return(theta)
}


AffyThetaOptDEP <- function(h.d, h.snp.gt.p, h.e.mu, h.snp.clust.p, out.p, theta, verbose=FALSE) {
	
	theta[["at"]] <- HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta, "at", list("lower"= 0.0, "upper"=0.2), "%", 1e-3, verbose=verbose)
	theta[["sigma.eta"]] <- HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta, "sigma.eta", list("lower"= 0.05, "upper"=1.0), ".", 1e-3, verbose=verbose)
	theta[["sigma.epsilon"]] <- HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta, "sigma.epsilon", list("lower"= 0.001, "upper"=0.5), "*", 1e-3, verbose=verbose)
	
	max.het.cov <- GetMaxHetCov(theta)
	if (max.het.cov > 0) {
		theta[["het.cov"]] = HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta, "het.cov", list("lower"= 0, "upper"=max.het.cov), "~", 1e-4, verbose=verbose)
	} else {
		theta[["het.cov"]] <- 0
	}
	
	theta[["nu"]] <- HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta, "nu", list("lower"= 2, "upper"=15), "&", 1e-1, verbose=verbose)
	theta[["bg"]] <- HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta, "bg", list("lower"= -0.1, "upper"=0.1), "@", 1e-4, verbose=verbose)
	
	return(theta)
}

GetMaxHetCov <- function(theta) {
  sigma.h <- GetSigmaH(theta[["sigma.epsilon"]], theta[["sigma.eta"]])
  return(theta[["sigma.eta"]]^2 / 2)
}

GridstartH1OptMeans <- function(idx, h.seg.dat, delta.tau, out.p, theta, n.grid=5, verbose=FALSE) {
  # browser()

	# mode.t <- e.mu[idx, 3]
	mode.t <- delta.tau[idx, 2]
	mode.d.vec <- rep(NA, n.grid + 1)
	ll.vec <- rep(NA, n.grid + 1)
	# mode.d.vec[1] <- abs(e.mu[idx, 2] - e.mu[idx,1])
	mode.d.vec[1] <- delta.tau[idx, 1]
	mode.d.vec[c(2:(n.grid + 1))] <- seq(0, mode.t, length=n.grid)
	for (i in 1:length(mode.d.vec)) {
    # i <- 1
		cur.par <- c(mode.d.vec[i], mode.t)
    ll.vec[i] <- AffyGetLL(cur.par, idx, h.seg.dat, out.p, theta )
	}
	mode.d.vec <- mode.d.vec[is.finite(ll.vec)]
	ll.vec <- ll.vec[is.finite(ll.vec)]
	
	start.ix <- which.min(ll.vec)
	if (length(start.ix) == 0)  {
		## no valid eval points
		ll <- -Inf
		while (is.finite(ll) == FALSE) {
			cur.par[2] <- runif(1, 0, 6)
			cur.par[1] <- runif(1, 0, cur.par[2])
			if (verbose) {
				cat("?")
			}
			# ll <- AffyGetLL(cur.par, h.seg.dat[["h.snp.d"]][[idx]], h.seg.dat[["h.cn.d"]][[idx]], out.p, h.seg.dat[["h.snp.gt.p"]][[idx]], theta )
      ll <- AffyGetLL(cur.par, idx, h.seg.dat, out.p, theta )
		}
	} else {
		cur.par <- c(mode.d.vec[start.ix], mode.t)
		if ((verbose) && (start.ix != 1)) {
			cat(start.ix)
		}
	}
  return(AffyDeltaTauOptum(cur.par, idx, h.seg.dat, out.p, theta))
}


GridstartH1OptMeansDEP <- function(d, e.mu, out.p, snp.gt.p, theta,
                                n.grid=5, verbose=FALSE) {
  
  mode.t <- e.mu[3]
  mode.d.vec <- rep(NA, n.grid + 1)
  ll.vec <- rep(NA, n.grid + 1)
  mode.d.vec[1] <- abs(e.mu[2] - e.mu[1])
  mode.d.vec[c(2:(n.grid + 1))] <- seq(0, mode.t, length=n.grid)
  for (i in 1:length(mode.d.vec)) {
    cur.par <- c(mode.d.vec[i], mode.t)
    ll.vec[i] <- GetLL(cur.par, d, out.p, snp.gt.p, theta )
  }
  mode.d.vec <- mode.d.vec[is.finite(ll.vec)]
  ll.vec <- ll.vec[is.finite(ll.vec)]

  start.ix <- which.min(ll.vec)
  if (length(start.ix) == 0)  {
    ## no valid eval points
    ll <- -Inf
    while (is.finite(ll) == FALSE) {
      cur.par[2] <- runif(1, 0, 6)
      cur.par[1] <- runif(1, 0, cur.par[2])
      if (verbose) {
        cat("?")
      }
      ll <- GetLL(cur.par, d, out.p, snp.gt.p, theta)
    }
  } else {
    cur.par <- c(mode.d.vec[start.ix], mode.t)
    if ((verbose) && (start.ix != 1)) {
      cat(start.ix)
    }
  }
  
  return(DeltaTauOptum(cur.par, d, out.p, snp.gt.p, theta))
}

AffyDeltaTauOptum <- function(cur.par, ix, h.seg.dat, out.p, theta) {
   snp.d <- h.seg.dat[["h.snp.d"]][[ix]] 
   cn.d <- h.seg.dat[["h.cn.d"]][[ix]]
   snp.gt.p <- h.seg.dat[["h.snp.gt.p"]][[ix]]

	res <- optim(cur.par, AffyGetLL, ix=ix, h.seg.dat=h.seg.dat, out.p=out.p, theta=theta )
	opt.d <- abs(res[["par"]][1])
	opt.t <- res[["par"]][2]
	
	return(c(opt.d, opt.t))
}


DeltaTauOptumDEP <- function(cur.par, d, out.p, snp.gt.p, theta) {

	res <- optim(cur.par, GetLL, d=d, out.p=out.p, snp.gt.p=snp.gt.p, theta=theta )
   	opt.d <- abs(res[["par"]][1])
   	opt.t <- res[["par"]][2]
        
   	return(AffyGetMeans(opt.d, opt.t))
}

