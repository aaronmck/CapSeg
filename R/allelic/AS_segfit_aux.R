## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

EStep <- function(d, delta.tau, snp.gt.p, out.p, theta) {
   clust.lik <- AffyGetSnpClustLik(d, delta.tau, theta)

   ## use priors for SNP allele frequency data
   snp.gt.p <- snp.gt.p * (1 - out.p)
   cp <- cbind(snp.gt.p, out.p ) 

   clust.lik <- clust.lik * cp
   snp.clust.p <- clust.lik / rowSums(clust.lik)

   if (any(is.na(snp.clust.p))) {
     stop("NA values found in snp.clust.p")
   }
   
   return(snp.clust.p)
}
