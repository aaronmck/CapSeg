## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CaptureHscrSegFit <- function(h.seg.dat, tol=1e-5, verbose=verbose)
{
   # Initialize Theta_CN
   # while
   #    for i in segs:
   #       optimize tau_i given Theta_CN
   #    optimize Theta given tau
   #    break when cond < 1e-4
   #   
   # Initialize Theta_AS
   # while
   #    for i in segs:
   #       optimize F_i given Theta_AS
   #    optimize Theta_AS given F
   #    break when cond < 1e-4
   #
   # for i in segs:
   #    optimize F_i | tau_i
   #    compute posterior on F_i
   #    compute posterior on tau_i
   #    compute posterior on F_i * tau_i and (1 - F_i * tau_i)
   # return F, tau, F*tau, (1-F) * tau, and all their error bars

   # h.seg.dat <- cap.dat[["as.res"]][["h.seg.dat"]]
   # tol=1e-5;

   CR.res = CaptureCNModelFit( h.seg.dat, verbose=verbose )
   tau = CR.res[["tau"]]
   Theta_CR = CR.res[["Theta"]]

   AS.res = CaptureASModelFit( h.seg.dat, verbose=verbose )
   wes.f = AS.res[["wes.f"]]
   Theta_AS = AS.res[["Theta"]]
   het.phase.prob = AS.res[["het.phase.prob"]]


   Theta = c(Theta_CR, Theta_AS)

   mu.minor = wes.f[,"f.hat"] * tau
   mu.major = (1-wes.f[,"f.hat"]) * tau

   seg.sems <- CalcAllelicSegSEMs(h.seg.dat, wes.f[,"f.hat"], wes.f[, "p.H0"], tau, Theta)

   mus = cbind( "mu1"=mu.minor, "mu2"=mu.major, tau, seg.sems )
   
##? Do we need the delta / tau stuff?
## Yes - because the plotting code expects it.
   delta.tau <- data.frame(delta=mus[,"mu2"] - mus[,"mu1"], tau=tau)
   delta.tau.sd <- data.frame(delta=mus[,"sigma1"]^2 + mus[,"sigma2"]^2, tau=mus[,"sigma3"])

   out <- list(delta.tau=delta.tau, delta.tau.sd=delta.tau.sd, Theta=Theta, wes.f=wes.f, cap.e.mu=mus, "het.phase.prob"=het.phase.prob )

   return(out)
}

CaptureASModelFit = function( h.seg.dat, verbose=FALSE )
{
   min.iter=1
   max.iter=10
   eps = 1e-4

   n.segs <- length(h.seg.dat[[1]])

   Theta <- InitCaptureASTheta(verbose=verbose)
   wes.f <- InitCaptureF( h.seg.dat, Theta ) 

   loglik <- -Inf
   iter <- 1
   while(TRUE) 
   {
      Theta <- CaptureASThetaOpt(h.seg.dat, wes.f, Theta, verbose=verbose)
      wes.f <- OptimizeCaptureF(h.seg.dat, wes.f, Theta, tol=1e-5, verbose=verbose)

      cur.loglik <- CalcCaptureASLogLik(h.seg.dat, wes.f[,"f.hat"], Theta)
      cond <- abs(cur.loglik - loglik) / abs(cur.loglik) 
      loglik <- cur.loglik

      if (verbose) {
#         PrintTheta(Theta, cond, sigma.h=FALSE)
         print(paste("f_skew = ", round(Theta[["f_skew"]],5), ", LogLik:", round(loglik, 4), ", cond = ", round(cond,5),  sep="" ))
      }
      if ((iter > min.iter) && ((cond < eps) || (iter >= max.iter))) {
         break
      }
      iter <- iter + 1
   }      

   het.phase.prob = list()
   for( i in 1:n.segs )
   {
      alt = h.seg.dat[["gh.wes.allele.d"]][[i]]["alt",]
      ref = h.seg.dat[["gh.wes.allele.d"]][[i]]["ref",]
      f.hat = wes.f[i,"f.hat"]
      het.phase.prob[[i]] = CapturePhaseProb( alt, ref, f.hat, Theta)
   }

   return( list("Theta"=Theta, "wes.f"=wes.f, "het.phase.prob"=het.phase.prob) )
}


CaptureASThetaOpt <- function(h.seg.dat, wes.f, Theta, verbose=FALSE) {
   LL = function(par)
   {
      Theta[["f_skew"]] = par
      if (verbose) cat("-")
      return( CalcCaptureASLogLik(h.seg.dat, wes.f[,"f.hat"], Theta) )
   }      

   if (verbose) print("Optimizing Capture allelic Theta")   

   res <- optimize(LL, lower=0.6, upper=1, tol=1e-4, maximum=TRUE )
   if(verbose) cat("\n")
   Theta[["f_skew"]] = res[["maximum"]]

   return(Theta)
}



CaptureCNModelFit = function( h.seg.dat, verbose=FALSE )
{
   min.iter=1
   max.iter=10
   eps= 1e-4

   Theta <- InitCaptureCRTheta(verbose=verbose)
   tau <- CaptureInitTau(h.seg.dat) 
   n.segs <- length(h.seg.dat[["h.capseg.d"]])
   loglik <- -Inf
   iter <- 1

   while(TRUE) 
   {  
      tau <- CaptureOptimizeTau(h.seg.dat, tau, Theta, tol=1e-4, verbose=verbose)
      Theta <- CaptureCNThetaOpt(h.seg.dat, tau, Theta, verbose=verbose)
      cur.loglik <- CalcCaptureCNLogLik(h.seg.dat, tau, Theta)

      cond <- abs(cur.loglik - loglik) / abs(cur.loglik) 
      loglik <- cur.loglik

      if (verbose) {
         print(paste("LogLik:", round(loglik, 4)))
         print(paste("sigma.scale.capseg = ", round(Theta[["sigma.scale.capseg"]],5), ", LogLik:", round(loglik, 4), ", cond = ", round(cond,5),  sep="" ))

      }
      if ((iter > min.iter) && ((cond < eps) || (iter >= max.iter))) {
         break
      }
      
      iter <- iter + 1
   }      

   return( list("Theta"=Theta, "tau"=tau) )
}



CaptureCNThetaOpt <- function(h.seg.dat, tau, Theta, verbose=FALSE) {
   LL = function(par)
   {
      Theta[["sigma.scale.capseg"]] = par
      if(verbose) cat("~")
      return( CalcCaptureCNLogLik(h.seg.dat, tau, Theta) )
   }      

   if (verbose) print("Optimizing Capture CR Theta")   

   res <- optimize(LL, lower=0, upper=20, tol=1e-4, maximum=TRUE)
   if(verbose) cat("\n")
   Theta[["sigma.scale.capseg"]] = res[["maximum"]]

   return(Theta)
}



## Is this used??
CapsegThetaOptExtreme <- function(h.seg.dat, delta.tau, out.p, Theta, verbose=FALSE) {
   
   if (verbose) print("Optimizing Capseg Theta")   
   Theta[["at.capseg"]] <- HThetaOptExtreme(h.seg.dat, delta.tau, out.p, Theta, "at.capseg", list("lower"= 0, "upper"=5), "+", 1e-4, probe.types=c("cap"), verbose=verbose)
   Theta[["sigma.scale.capseg"]] <- HThetaOptExtreme(h.seg.dat, delta.tau, out.p, Theta, "sigma.scale.capseg", list("lower"= 0, "upper"=20), "~", 1e-4, probe.types=c("cap"), verbose=verbose)
   # Theta[["nu.capseg"]] <- HThetaOptExtreme(h.seg.dat, h.e.mu, out.p, Theta, "nu.capseg", list("lower"= 1, "upper"=25), "-", 1e-1, probe.types=c("cap"), verbose=verbose)
   
   return(Theta)
}

InitCaptureCRTheta <- function(verbose=FALSE) 
{
   if (verbose) print("Initializing Capture-only Theta")
   Theta <- list()
   Theta[["at.capseg"]] <- 0
   Theta[["sigma.scale.capseg"]] <- .15    ## this is the only free param
   Theta[["CR.out.p"]] = 0.05

   return(Theta)
}

InitCaptureASTheta = function(verbose=FALSE)
{
   Theta = list()
# linear model of log(rho) = m*f_skew/2 + b

# fit with out.p = 0.001
   Theta[["b"]] = -33.67232
   Theta[["m"]] = 82.77464  

# fit with out.p = 0.005
#   Theta[["b"]] = -45.86149
#   Theta[["m"]] = 110.51106 

   Theta[["f_skew"]] = 0.48 * 2     ## this is the only free param
   Theta[["allelic_out.p"]] = 0.005

   return(Theta)
}



CapturePlatformSpecificOptimizationExtreme <- function(idx, delta.tau, h.seg.dat, out.p, Theta, verbose=FALSE) {
   ## platform specific optimization - the optimization
   ## we do for the affy (and arrays in general)
   
   return(GridstartH1OptMeansExtreme(idx, h.seg.dat, delta.tau, out.p, Theta, verbose=verbose))
}


F_opt_func <- function( Par, alt, ref, Theta, dom )
{
   f = Par
#   if (f < dom[1] ) { f=dom[1] }
#   if( f > dom[2] ) { f=dom[2] }

   if (f < dom[1] | f > dom[2] ) { return(-1e99) }

   ll =  CalcCaptureSegAllelicLogLik(alt, ref, f, Theta ) 
   if(ll== -Inf) { ll = -1e100}
   return( ll )
}

# For nlm-based optimizer
neg_F_opt_func <- function( Par, alt, ref, Theta, dom )
{
   return( -F_opt_func( Par, alt, ref, Theta, dom ) )
}


InitCaptureF = function( h.seg.dat, Theta, verbose=FALSE ) {
   ## Modify Theta[["allelic_out.p"]] to low tolerance for outliers.  For provisional fitting
   Theta[["allelic_out.p"]] = 0
   n.segs = length(h.seg.dat[["gh.wes.allele.d"]])

   fdom = c(0,0.5)
   wes.f= matrix( c(NA, 0.5, 0.5), nrow=n.segs, ncol=3, byrow=TRUE)
   f.hat = rep(NA, n.segs)
   colnames(wes.f) = c("f.hat", "p.H0", "p.H1")

   if (verbose) print(paste("Initializing f.hat for ", n.segs, " segs:", sep="") )

   for( s in 1:n.segs )
   {
      alt <- h.seg.dat[["gh.wes.allele.d"]][[s]]["alt", ]
      ref <- h.seg.dat[["gh.wes.allele.d"]][[s]]["ref", ]

      f_grid = seq(fdom[1], fdom[2], length=10)
      opt_f = array(NA, dim=c( length(f_grid), 2))
      for (i in seq_along(f_grid))
      {
# faster
         opt = try( nlm(f=neg_F_opt_func, p=f_grid[i], alt=alt, ref=ref, Theta=Theta, dom=fdom) )
         if( class(opt) != "try-error" )
         {
            opt_f[i, ] = c(opt$estimate, -opt$minimum)
         } 
         else
         {
            opt <- optim(par=f_grid[i], fn=F_opt_func, method="L-BFGS-B", lower=0, upper=0.5, control=list(fnscale=-1), alt=alt, ref=ref, Theta=Theta, dom=fdom)
            opt_f[i,] = c( opt[["par"]], opt[["value"]] )
         }
      }
      if (verbose) cat("@")
#      f.hat[s] = f_grid[which.max( t(opt_f[,2]))]
      m.ix = which.max( t(opt_f[,2]))
      f.hat[s] = opt_f[ m.ix, 1 ]
   }
   if(verbose) cat("\n")
   wes.f[,"f.hat"] = f.hat

   return(wes.f)
}



OptimizeCaptureF <- function(h.seg.dat, wes.f, Theta, tol=1e-5, verbose=FALSE) 
{
   n.segs = length(h.seg.dat[["gh.wes.allele.d"]])
   prev_f = wes.f[,"f.hat"]
   res = matrix( NA, nrow=n.segs, ncol=3)
   fdom = c(0,0.5)

#   res = foreach(r=1:n.segs, .combine=rbind) %dopar% 
   for(r in 1:n.segs)
   {
      # r = 54
      n.hets = ncol(h.seg.dat[["gh.wes.allele.d"]][[r]])
      if( (n.hets > 0) & (sum(h.seg.dat[["gh.wes.allele.d"]][[r]]) > 0) ) 
      {
         alt <- h.seg.dat[["gh.wes.allele.d"]][[r]]["alt", ]
         ref <- h.seg.dat[["gh.wes.allele.d"]][[r]]["ref", ]

#         H1.opt.res <- optimize(F_opt_func, interval=c(0, 0.5), tol=tol, maximum=TRUE)
#          H1.f.hat= H1.opt.res[["maximum"]]
#         H1.opt.f.loglik = H1.opt.res[["objective"]]

         opt = nlm(f=neg_F_opt_func, p=prev_f[r], alt=alt, ref=ref, Theta=Theta, dom=fdom )
         H1.f.hat = opt$estimate
         H1.opt.f.loglik = opt$minimum

#         opt <- optim(par=prev_f[r], fn=F_opt_func, method="L-BFGS-B", lower=0, upper=0.5, control=list(fnscale=-1), alt=alt, ref=ref, Theta=Theta, dom=fdom)
#         H1.f.hat = opt[["par"]]
#         H1.opt.f.loglik = opt[["value"]]

         # Let H0 be for F=0.5, and H1 be for F=f
         H0.f.hat = 0.5
         ev.H0 = CalcCaptureSegAllelicLogLik(alt, ref, f=H0.f.hat, Theta )
         ev.H1 <- GridEstimateOfFNormalizingConstant(alt, ref, Theta, dx=1e-2)

         ev.sum = LogAdd(c(ev.H0, ev.H1))
         log.p.H0 <- ev.H0 - ev.sum
         log.p.H1 <- ev.H1 - ev.sum

         if( n.hets > 10 )  ## model averaging
         {
            f.hat <- exp(log.p.H0) * H0.f.hat + exp(log.p.H1) * H1.f.hat
         }
         else              ## Winner take all
         {
            f.hat = ifelse(log.p.H0 > log.p.H1, H0.f.hat, H1.f.hat )
         }

         res[r,] = c(f.hat=f.hat, p.H0=exp(log.p.H0), p.H1=exp(log.p.H1))

         if( FALSE & verbose) 
         {
            cat("Highest Evidence F Model: ")
            if (ev.H0 > ev.H1) cat("H0\n") else cat("H1\n")
            print( paste("H1 fhat = ", round(H1.f.hat,5), sep=""))
            print( paste("p.H0 = ", round(exp(log.p.H0),5), sep=""))
            print( paste("H0 loglik = ", round(ev.H0,5), sep=""))
            print( paste("H1 opt f loglik = ", round(H1.opt.f.loglik,5), sep=""))
         }
         if( verbose) 
         {  
            if (ev.H0 > ev.H1) cat("0") else cat("1")
         }         

      }  else {
         res[r,] = c(f.hat=NA, p.H0=NA, p.H1=NA)
      }
   }

  if( verbose) 
  {  
     cat("\n")
  }         

   rownames(res) <- NULL
   colnames(res) = colnames(wes.f)

   return(res)
}

GridEstimateOfFNormalizingConstant <- function(alt, ref, Theta, dx=1e-2) {
   prior <- 2
   likelihood <- function(x) sapply(x, function(x) CalcCaptureSegAllelicLogLik(alt, ref, f=x, Theta ) )
   post.numerator <- function(x) sapply(x, function(x) log(prior) + likelihood(x) )
   y <- post.numerator(seq(0, 0.5, by=dx))
   norm.const <- LogAdd(y + log(dx))
   return(norm.const)
}



CaptureOptimizeTau <- function(h.seg.dat, tau, Theta, tol=1e-4, verbose=FALSE) 
{
   
   if (verbose) print("Optimizing Tau: ")
   lb <- 0
   LL <- function(par, d, Theta){
      
      par <- max(lb, par)
      tau = par
      CalcCaptureSegTauLogLik(d, tau, Theta)
   }
 
   n.segs = length(h.seg.dat[["h.capseg.d"]])
   opt.tau = rep(NA, n.segs)
   # for( i in 1:n.segs ) {
   #    d <- h.seg.dat[["h.capseg.d"]][[i]]
   #    if (length(as.vector(d)) == 0) return(NA)
#      if (length(as.vector(d)) == 0) { opt.tau[i] = NA; next }
   #    opt.tau[i] <- optim(par=tau[i], fn=LL, method="L-BFGS-B", lower=0, upper=Inf, control=list(factr=tol, fnscale=-1), d=d, Theta=Theta)[["par"]]
   #    if (verbose) cat(paste(i, sep=""), ".")
   # }

   opt.tau <- foreach(i=1:n.segs, .combine=c) %dopar% {
   # opt.tau <- sapply(1:n.segs, function(i) {
      d <- h.seg.dat[["h.capseg.d"]][[i]]
      if (length(as.vector(d)) == 0) { return(NA) }
      opt.tau.i <- optim(par=tau[i], fn=LL, method="L-BFGS-B", lower=lb, upper=Inf, control=list(factr=tol, fnscale=-1), d=d, Theta=Theta)[["par"]]
      if (verbose) cat(paste(i, sep=""), ".")
      return(opt.tau.i)
   }

   if( verbose ){ cat("\n") }

   return(opt.tau)
}
# i <- 73
# d <- h.seg.dat[["h.capseg.d"]][[i]]
# curve( sapply(x, function(x) LL(x, d, Theta)), from=0, to=4)
# LL(-4.44089209850063e-16, d, Theta)


CaptureInitTau <- function(h.seg.dat) 
{
   n.segs = length(h.seg.dat[["h.capseg.d"]])
   tau = rep(NA, n.segs)
   for(i in 1:n.segs )
   {
      d <- h.seg.dat[["h.capseg.d"]][[i]]
      if( length(as.vector(d) ) > 0) { tau[i] = median(d) }
   }

   return(tau)
}

