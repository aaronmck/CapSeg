## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

# move into affy specific -- fix me
## FIXME: what's the deal w/ the above coment? this is called
## from non-affy specific location and I don't see another declaration

SegPostSdExtreme <- function(snp.d, cn.d, delta.tau, out.p, h.snp.gt.p, theta) {
	CobCalcLL <- function(par) {
		dist <- par[1]
		total <- par[2]
		e.mu <-  c(((total / 2) - (dist / 2)), ((total / 2) + (dist / 2)), total)
		
		logliks = c(
				sum(CalcSnpLogLik(snp.d, e.mu, out.p, h.snp.gt.p, theta)),
				sum(CalcCnLogLik(cn.d, e.mu, out.p, theta))
#				sum(CalcCaptureLogLik(capseg.d, e.mu, out.p, theta))
				)
		logliks = logliks[!is.nan(logliks)]
		res = -sum( logliks)
		return(res)
	}
	h1.mode <- GetMeans(delta.tau[1], delta.tau[2])
	m.delta <- abs(h1.mode[1] - h1.mode[2])
	m.tau <- h1.mode[3]
	hess.mat <- hessian(CobCalcLL, x=c(m.delta, m.tau))
	ihm <- try(solve(hess.mat), silent=TRUE)
	if (inherits(ihm, "try-error")) {
		sd.mat <- matrix(NA, nrow=2, ncol=2)
	} else {
		sd.mat <- sqrt(abs(ihm))
	}
	se <- c(sd.mat[1,1], sd.mat[2,2]) 
	names(se) <- c("delta", "tau")
	return(se)
}


SegPostSd <- function(d, e.mu, out.p, het.prob, theta) {
  CobCalcLL <- function(par) {
    dist <- par[1]
    total <- par[2]
    e.mu <-  c(((total / 2) - (dist / 2)), ((total / 2) + (dist / 2)), total)
    res <- -sum(CalcSnpLogLik(d, e.mu, out.p, het.prob, theta))

    return(res)
  }

  h1.mode <- e.mu
  m.delta <- abs(h1.mode[1] - h1.mode[2])
  m.tau <- h1.mode[3]

  hess.mat <- hessian(CobCalcLL, x=c(m.delta, m.tau))
  ihm <- try(solve(hess.mat), silent=TRUE)
  if (inherits(ihm, "try-error")) {
    sd.mat <- matrix(NA, nrow=2, ncol=2)
  } else {
    sd.mat <- sqrt(abs(ihm))
  }

  se <- c(sd.mat[1,1], sd.mat[2,2]) 
  names(se) <- c("delta", "tau")

  return(se)
}

HThetaOptExtreme <- function(h.seg.dat, delta.tau, out.p, theta, parname, limits, symbol, opttol, probe.types, verbose=FALSE) {

	LL <- function(par) {
		theta[parname] <- par
		
		delta.tau.rows <- nrow(delta.tau)
		seg.loglik <- rep(0, delta.tau.rows)
		for (s in 1:delta.tau.rows) {
			
			logliks = c( sum(CalcSnpLogLik(h.seg.dat[["h.snp.d"]][[s]], delta.tau[s, ], out.p, h.seg.dat[["h.snp.gt.p"]][[s]], theta)), 
							sum( CalcCnLogLik(h.seg.dat[["h.cn.d"]][[s]], delta.tau[s, ], out.p, theta)), 
							sum( CalcCaptureLogLik(h.seg.dat[["h.capseg.d"]][[s]], delta.tau[s, 2], out.p, theta ))
							)
			seg.loglik[s] = sum(logliks[c("snp", "cn", "cap") %in% probe.types])
			
		}
		cur.loglik <- sum(seg.loglik[complete.cases(seg.loglik)])
		
		if (verbose) {
			if (!is.finite(cur.loglik)) {
				print(paste(parname, ": Non-finite log-liklihood!", sep=""))
			} else {
				cat(symbol)
			}
		}
		
		return(-cur.loglik)
	}
	
	res <- optimize(LL, lower=limits[["lower"]], upper=limits[["upper"]], tol=opttol, maximum=FALSE)
	
	return(res[["minimum"]])
}


HThetaOpt <- function(h.d, h.snp.gt.p, h.e.mu, out.p, theta, parname, limits, symbol, opttol, verbose=FALSE) {
	LL <- function(par) {
		theta[parname] <- par
		
		h.e.mu.rows <- nrow(h.e.mu)
		seg.loglik <- rep(0, h.e.mu.rows)
		
		for (s in 1:h.e.mu.rows) {
			snp.loglik <- CalcSnpLogLik(h.d[[s]], h.e.mu[s, ], out.p, h.snp.gt.p[[s]], theta)
			seg.loglik[s] <- sum(snp.loglik)
		}
		cur.loglik <- sum(seg.loglik)
		
		if (verbose) {
			if (!is.finite(cur.loglik)) {
				print(paste(parname, ": Non-finite log-liklihood!", sep=""))
			} else {
				cat(symbol)
			}
		}
		
		return(-cur.loglik)
	}
	
	res <- optimize(LL, lower=limits[["lower"]], upper=limits[["upper"]], tol=opttol, maximum=FALSE)
	
	return(res[["minimum"]])
}

