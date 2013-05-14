## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


HscrSegFitExtreme <- function(h.seg.dat, theta, eps=1e-5, out.p=1e-3, min.iter=1, max.iter=10, force.diploid=FALSE, verbose=FALSE) {
   
   if (force.diploid) {
      min.iter = 1
   }
   
   if (verbose) {
      print(paste("OUT_P is", out.p))
   }
   smms = LoadCached({
            smms <- SegMeansMStepsExtreme(h.seg.dat, out.p, theta, eps, force.diploid, min.iter, max.iter, verbose=verbose)
            saveRDS (smms, file = file.path(RESULTS.DIR, "smms.res.rds"))
            smms }, overwrite=TRUE, res.fn = file.path(RESULTS.DIR, "smms.res.rds"), mod.name="SegMeansMStepsExtreme")
   
   theta <- smms[["theta"]]
   h.snp.clust.p <- smms[["h.snp.clust.p"]]
   delta.tau <- smms[["delta.tau"]]
   
   wes.f = OptimizeCaptureF(h.seg.dat, out.p=0 , verbose=verbose)
   het.phase.log.p = lapply(1:length(h.seg.dat[["gh.wes.allele.d"]]), function(i) 
            CapturePhaseProb( alt=h.seg.dat[["gh.wes.allele.d"]][[i]]["alt",], 
               ref= h.seg.dat[["gh.wes.allele.d"]][[i]]["ref",], wes.f[i, "f.hat"], out.p=0))
   
   cap.e.mu = CalcCaptureMu(h.seg.dat, wes.f[, "f.hat"])
   
   delta.tau.sd <- BuildHMuSdExtreme(h.seg.dat, delta.tau, out.p, theta, verbose=verbose)
   seg.log.ev <- BuildSegLogEvExtreme(h.seg.dat, theta, out.p, delta.tau, verbose=verbose)
   seg.expected.phase <- BuildSegExpectedPhase(length(h.seg.dat[["h.snp.d"]]), h.snp.clust.p)
   
   return(list(snp.clust.p=h.snp.clust.p, h.snp.gt.p=h.seg.dat[["h.snp.gt.p"]],
               delta.tau=delta.tau,
               delta.tau.sd=delta.tau.sd, sigma.h=smms[["sigma.h"]], loglik=smms[["loglik"]],
               seg.log.ev=seg.log.ev, seg.expected.phase=seg.expected.phase,
               theta=theta, wes.f=wes.f, het.phase.log.p=het.phase.log.p, cap.e.mu=cap.e.mu))
}


SegMeansMStepsExtreme <- function(h.seg.dat, out.p, theta, eps, force.diploid, min.iter, max.iter, verbose=FALSE) {

   # h.seg.dat <- iams.res[["as.res"]][["h.seg.dat"]]
   # theta <- iams.res[["em.fit"]][["theta"]]
   # eps <- iams.res[["use.eps"]]
   # min.iter <- 1
   # max.iter <- 10

   n.segs <- length(h.seg.dat[[1]])
   if (verbose) {
      print(paste("h.d size =", n.segs))
   }
   h.ngt.t <- array(NA, dim=c(n.segs, 5))
   delta.tau <- InitDeltaAndTau(h.seg.dat, theta, seg.info, force.diploid)
   
   loglik <- -Inf  
   iter <- 1
   
   while (TRUE) {
      ##  M-steps for seg-means
      delta.tau <- OptimizeDeltaTauArray(delta.tau, h.seg.dat, out.p, theta, force.diploid, verbose=verbose)
      
      h.snp.clust.p <- ClustProbsESteps(h.seg.dat[["h.snp.d"]], delta.tau, h.seg.dat[["h.snp.gt.p"]], out.p, theta)
      for (i in seq(n.segs)) {
         h.ngt.t[i, ] <- colSums(h.snp.clust.p[[i]])
      }
      
      p.het <- CalculatePHet(h.ngt.t)
      
      theta <- HscrSegFitThetaOptExtreme(iter, min.iter, p.het, h.seg.dat, delta.tau, out.p, theta, verbose=verbose)
      
      cur.loglik <- GetCurLogLikExtreme(h.seg.dat, delta.tau, out.p, theta)
      
      cond <- abs(cur.loglik - loglik) / abs(cur.loglik) 
      loglik <- cur.loglik
      sigma.h <- GetSigmaH(theta[["sigma.epsilon"]], theta[["sigma.eta"]])
      
      if (verbose) {
         print(round(cur.loglik, 4))
         PrintTheta(theta, cond)
      }
      
      
      if ((iter > min.iter) && ((cond < eps) || (iter >= max.iter))) {
         break
      }
      
      iter <- iter + 1
   }
   
   return(list(loglik=cur.loglik, theta=theta, delta.tau=delta.tau, sigma.h=sigma.h, h.snp.clust.p=h.snp.clust.p))
}


   
InitDeltaAndTauExtreme <- function(h.seg.dat, theta, seg.info, force.diploid) {
   
   ## just return delta and tau
   
    n.segs <- length(h.seg.dat[[1]])
   if (force.diploid) {
      h.e.mu <- array(NA, dim=c(n.segs, 3))
      for (i in seq(n.segs)) {
         h.e.mu[i, ] <- c(1, 1, 2)
      }
   } else {
      init <- function(x) {
         d <- foreach(i=c("h.snp.d", "h.cn.d", "h.capseg.d"), .combine=c) %do% {colSums(h.seg.dat[[i]][[x]])}
         out = median(d)
         out = ifelse(is.na(out), 2, out)
         return(median(out))   
      }
      h.mu.t <- unlist(lapply(1:n.segs, init))
      
      h.mu.t[h.mu.t < 0] <- 0
      i.mu.t <- h.mu.t
      quart <- i.mu.t / 4
      # h.e.mu <- cbind(quart, i.mu.t - quart, i.mu.t)
      
      delta = 2 * (i.mu.t - quart) - i.mu.t
      tau = i.mu.t
      dt = cbind(delta, tau)
   }
   return(dt)
}
