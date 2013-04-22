## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


HTx <- function(x, sigma.epsilon, sigma.eta, alpha=0) {
 # if (is.matrix(x)) {
 #   out <- .Call("HTxMatrix", x, sigma.epsilon, sigma.eta, alpha, PACKAGE="HAPSEG")
 # } else {
 #   out <- .Call("HTxVector", x, sigma.epsilon, sigma.eta, alpha, PACKAGE="HAPSEG")
 # }
 # return(out)

  b = (sqrt( exp(sigma.eta^2) - 1) ) / sigma.epsilon
  a = -alpha * b
  asinh(a+b*x)

}

HTxVar <- function(x, sigma.epsilon, sigma.eta, alpha=0) {
#  return(.Call("HTxVar", x, sigma.epsilon, sigma.eta, alpha, PACKAGE="HAPSEG"))
   b = (sqrt( exp(sigma.eta^2) - 1) ) / sigma.epsilon
   a = -alpha * b

  (b/(sqrt(1+(a+b*x)^2)) )


}

HTxInv <- function(x, sigma.epsilon, sigma.eta, alpha=0) {
  b <- (sqrt(exp(sigma.eta^2) - 1)) / sigma.epsilon
  a <- -alpha * b
  
  return((sinh(x) -a) / b)
}

GetSigmaH <- function(sigma.epsilon, sigma.eta, alpha=0) {
  return(sigma.epsilon * HTxVar(alpha, sigma.epsilon, sigma.eta))
}

