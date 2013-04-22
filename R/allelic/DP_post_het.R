## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CalcDpLoglik <- function(N, K, alpha, coefs) {
  ## actually returns the log of the normalized probability 
  ## Escobar and West 1995 for DP loglik.  No correction for segment sizes!
  ## Carter and Getz 2008 with seglen correction.
  return(log(coefs[K]) + K * log(alpha) - log(EvalMPolynom(alpha, coefs)))
}

EvalMPolynom <- function(alpha, coefs) {
  pows <- c(1:length(coefs))
  return(sum(coefs * alpha^pows))
}
