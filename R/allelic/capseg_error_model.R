CalcAllelicSegSEMs <- function(h.seg.dat, f.hat, f.p.H0, tau, Theta) 
{
   # h.seg.dat <- cap.res[["as.res"]][["h.seg.dat"]]
   # f.hat <- F[, "f.hat"]
   # f.p.H0 <- F[, "p.H0"]

   n.seg = length(h.seg.dat[["h.capseg.d"]])

   cols = c("sigma1", "sigma2", "sigma3")
   result = matrix( NA, nrow=n.seg, ncol=length(cols) )
   colnames(result)=cols

   ab = CalcExpectedAlphaBeta(h.seg.dat, f.hat, Theta)
   for(i in 1:n.seg )
   {
      # i = 1; 
      d = h.seg.dat[["h.capseg.d"]][[i]]
      if (length(d) < 2 ) {
         result[i,] = c(NA, NA, sigma.mu)
         next
      } 

      sigma <- tau[i] * Theta[["sigma.scale.capseg"]]
      mu.hat = tau[i]

## TODO: refactor ?
      loglik <- function(x) sapply(x, function(n) sum(dnorm(d, mean=n, sd=sigma, log=T)))
      
      C = -hessian(loglik, c(tau[i]), "Richardson")[1]
      sigma.mu = C^(-1/2) 

      if (is.na(f.hat[i])) {
         result[i,] = c(NA, NA, sigma.mu)
         next
      }
      
      ## Then calc e.mu[1] and e.mu[2] using e.mu[3], f, and 1-f, and calc their error bars. 
      ## Gaussian approximation
      n <- 10000
      n.H0 <- round(n * f.p.H0[i])
      H0.f.sample <- rep(0.5, n.H0)
      H1.f.sample <- rbeta(n - n.H0, ab$alpha[i], ab$beta[i])
      f.sample <- c(H0.f.sample, H1.f.sample)[sample.int(n, n)] # mix up sample
      tau.sample <- rnorm(n, mean=mu.hat, sd=sigma.mu)
      sample.minor <- f.sample * tau.sample 
      sample.major <- (1 - f.sample) * tau.sample
#      mu.minor <- mean(sample.minor)
#      mu.major <- mean(sample.major)
      sigma.minor <- sd(sample.minor)
      sigma.major <- sd(sample.major)
      
      result[i,] = c(sigma.minor, sigma.major, sigma.mu)
   }
   rownames(result) <- NULL
   return(result)
}



CalcExpectedAlphaBeta <- function(h.seg.dat, f.hat, Theta) { 

   foreach (i=1:length(h.seg.dat[["gh.wes.allele.d"]]), .combine=rbind) %dopar% {

      d = h.seg.dat[["gh.wes.allele.d"]][[i]]
      
      alt.phase.prob = CapturePhaseProb( alt=d["alt",], ref= d["ref",], f.hat[i], Theta)[,1]
      e.alpha.seg = 1 + sum(d["alt",] * alt.phase.prob) + sum(d["ref",] * (1 - alt.phase.prob))
      e.beta.seg = 1 + sum(d["ref",] * alt.phase.prob) + sum(d["alt",] * (1 - alt.phase.prob))
      data.frame(alpha=e.alpha.seg, beta=e.beta.seg, stringsAsFactors=F)
   }
   
}

CalcCaptureFErrorBars <- function(h.seg.dat, f.hat, conf=.95) {
   ab = CalcExpectedAlphaBeta(h.seg.dat, f.hat)
   eb = foreach (i=1:length(h.seg.dat[["gh.wes.allele.d"]]), .combine=rbind) %dopar% {
      lb = qbeta((1 - conf) / 2, ab[i,"alpha"], ab[i, "beta"])
      ub = qbeta((1 + conf) / 2, ab[i,"alpha"], ab[i, "beta"])
      data.frame(lb=lb, f.hat=f.hat[i], ub=ub)
   }
   return(eb)
}

CapturePhaseProb <- function(alt, ref, f, Theta) 
{
   mat = CaptureAllelicLoglikMat( alt, ref, f, Theta )
   probs = exp(mat - LogAdd( mat ))
   return(probs)
}

CaptureAllelicLoglikMat = function( alt, ref, f, Theta )
{
   cov = alt + ref
   # rho = effective sample size for beta
   # linear model for log rho vs. f_skew, learned from normal samples.
   f_skew = Theta[["f_skew"]]
   rho = exp(Theta[["m"]] * (f_skew/2) +  Theta[["b"]])
   out.p = Theta[["allelic_out.p"]]

#   A1 = f_skew*f*rho
#   B1 = (1/f_skew)*(1-f)*rho
   A1 = f_skew*f*rho
   B1 = (1-f_skew*f)*rho
   lik1 = log(1/2) + log(1-out.p) + d_beta_binom(alt, A1, B1, cov, log=TRUE)

#   A2 = f_skew*(1-f)*rho
#   B2 = (1/f_skew)*f*rho
   A2 = (1 - f_skew^-1 * f) * rho
   B2 = f_skew^-1 * f * rho
   lik2 =  log(1/2) + log(1-out.p) + d_beta_binom(alt, A2, B2, cov, log=TRUE)

## ideal model for unbiased data (WGS?)
## This is equivalent to above with lim rho -> Inf, f_skew=1 
   # lik1 = log(1/2) + log(1-out.p) + dbinom(alt, cov, f, log=TRUE)
   # lik2 =  log(1/2) + log(1-out.p) + dbinom(alt, cov, 1-f, log=TRUE)
   
   outlier = log(out.p )  ## uniform outlier model
   mat = cbind(lik1, lik2, outlier)

   return(mat)
}

CalcCaptureSegAllelicLogLik <- function(alt, ref, f, Theta )
{
   mat = CaptureAllelicLoglikMat( alt, ref, f, Theta )

   ll = LogAdd(mat)
   nan.idx = is.nan(ll)
   ll[nan.idx] <- -Inf

   return(sum(ll))
}


CalcCaptureASLogLik <- function(h.seg.dat, f, Theta) 
{
   n.segs <- length(h.seg.dat[["gh.wes.allele.d"]])
   ll <- sapply(1:n.segs, function(i) {
         d <- h.seg.dat[["gh.wes.allele.d"]][[i]]
         alt = d["alt",]
         ref = d["ref",]
         CalcCaptureSegAllelicLogLik(alt, ref, f[i], Theta)
   })      
   return(sum(ll))
}


CalcCaptureSegTauLogLik <- function(d, tau, Theta) 
{
   if (length(d) == 0 ) return(NaN)
   # Gaussian Model
   atten.mu3 = Atten(tau, Theta[["at.capseg"]])
   sigma <- tau * Theta[["sigma.scale.capseg"]]
   out.p = Theta[["CR.out.p"]]

   lik = log(1 - out.p) + dnorm(d, atten.mu3, sigma, log=TRUE)
   outlier = rep(log(out.p * ( 1 / 5 ) ), length(lik))  ## outlier model is uniform on [0-5]?
   
   # Scaled T-Distribution Model
   # lik = log(1 - out.p) + DScaledT(d, e.mu[3], Theta[["sigma.scale.capseg"]], Theta[["nu.capseg"]], log=TRUE)
   # outlier = rep(log(out.p * (1 / 5)), length(lik))
   
   mat = cbind(lik[1,], unlist(outlier))
   ll = sum(LogAdd(mat))
   if(is.nan(ll)) { ll = -Inf }
  
   return(ll)
}


CalcCaptureCNLogLik <- function(h.seg.dat, tau, Theta) {
   n.segs <- length(h.seg.dat[["h.capseg.d"]])
   ll <- sapply(1:n.segs, function(i) {
         d <- h.seg.dat[["h.capseg.d"]][[i]]
         CalcCaptureSegTauLogLik(d, tau[i], Theta)
   })      
   return(sum(ll))
}


