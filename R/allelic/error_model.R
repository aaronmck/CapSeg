## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

IsAffyPlatform <- function(platform) {
  return(platform %in% c("SNP_6.0", "SNP_250K_STY"))
}
  
SetPlatformSpecificFuncs <- function(platform) {
  ## There are a variety of functions which vary depending on the
  ## platform. For now, dynamically assign those based on platform
  
  ## FIXME: This is pretty dangerous, set up a better way of handling
  ## this (which doesn't involve branching everywhere or passing platform
  ## around constantly)

  if (IsAffyPlatform(platform)) {
    return(SetAffyPlatformSpecificFuncs())
  } else if (platform == "WES") {
    return(SetCapturePlatformSpecificFuncs())
  } else {
    stop("Unsupported platform: ", platform)
  }
}


SetSegPlatformSpecificFuncs <- function() {
  InitTheta <<- SegInitTheta
  ThetaOpt <<- SegThetaOpt
  Atten <<- SegAtten
  TauPlatformSpecificInitialization <<- SegTauPlatformSpecificInitialization
  PlatformSpecificOptimization <<- SegPlatformSpecificOptimization
  GetMeans <<- SegGetMeans
  GetLL <<- SegGetLL
  CalcSnpLogLik <<- SegCalcSnpLogLik
  GetSnpClustLik <<- SegGetSnpClustLik
  DmvFunc <<- SegDmvFunc

  return(TRUE)
}

SetCapturePlatformSpecificFuncs <- function() {
  Atten <<- AffyAtten
  InvAtten <<- AffyInvAtten
  GetMeans <<- AffyGetMeans
  CalcSnpLogLik <<- AffyCalcSnpLogLik
  CalcCnLogLik <<- AffyCalcCnLogLik
  GetSnpClustLik <<- AffyGetSnpClustLik
  DmvFunc <<- AffyDmvFunc

  return(TRUE)
}

# SetAffyPlatformSpecificFuncsDEP <- function() {
#   InitTheta <<- AffyInitTheta
#   InitThetaExtreme <<- AffyInitThetaExtreme
#   ThetaOpt <<- AffyThetaOpt
#   ThetaOptExtreme <<- AffyThetaOptExtreme
#   InitCnBg <<- AffyInitCnBg
#   Atten <<- AffyAtten
#   InvAtten <<- AffyInvAtten
#   TauPlatformSpecificInitialization <<- AffyTauPlatformSpecificInitialization
#   PlatformSpecificOptimization <<- AffyPlatformSpecificOptimization
#   PlatformSpecificOptimizationExtreme <<- AffyPlatformSpecificOptimizationExtreme
#   GetMeans <<- AffyGetMeans
#   GetTau <<- AffyGetTau
#   GetDelta <<- AffyGetDelta
#   GetLL <<- AffyGetLL
#   DeltaTauOptumExtreme <<- AffyDeltaTauOptumExtreme
#   CalcSnpLogLik <<- AffyCalcSnpLogLik
#   CalcCnLogLik <<- AffyCalcCnLogLik
#   GetSnpClustLik <<- AffyGetSnpClustLik
#   DmvFunc <<- AffyDmvFunc

#   return(TRUE)
# }

# SetAffyAndCaptureFuncs <- function() {
#   Atten <<- AffyAtten
#   InvAtten <<- AffyInvAtten
#   PlatformSpecificOptimization <<- AffyPlatformSpecificOptimization
#   GetMeans <<- AffyGetMeans
#   CalcSnpLogLik <<- AffyCalcSnpLogLik
#   CalcCnLogLik <<- AffyCalcCnLogLik
#   GetSnpClustLik <<- AffyGetSnpClustLik
#   DmvFunc <<- AffyDmvFunc

#   return(TRUE)
# }
  
DFunc <- function(x, mu, sigma, nu) {
  return(d_scaled_t(x, mu, sigma, nu))
}

ExomeDFunc = function(x, tau, theta) {
	return(dnorm(x, mean=tau, sd=tau * theta[["sigma.scale.capseg"]]))
}


DoCalcSnpLogLik <- function(d, delta.tau, out.p, snp.gt.p, theta) {
  ## This function calculates the SNP log likelihood given the
  ## following parameters:
  ## d - the input data, containing a row for each allel (A & B)
  ## e.mu - the error model
  if (length(d) == 0 ) return(0)
  clust.lik <- AffyGetSnpClustLik(d, delta.tau, theta)

  snp.gt.p <- snp.gt.p * (1 - out.p)
  cp <- cbind(snp.gt.p, out.p) 
  
  return(log(rowSums(cp * clust.lik)))
}
