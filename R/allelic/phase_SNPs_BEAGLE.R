## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

PhaseSnps <- function(h.snp.clust.p, h.snp.gt.p, h.snp.annot, platform, tmp.dir, plate.name, phased.bgl.dir, verbose=FALSE) {
  ## kMinSnps was originally passed in but was passed in as a constant
  ## value so redefining it here
  kMinSnps <- 100
  
  n.seg <- length(h.snp.clust.p)

  new.h.snp.gt.p <- list()
  unfixed.h.snp.clust.p <- list()
  h.switch.ix <- list()

  ## FIXME: This could be foreached if we wanted to ...
  for (s in seq(n.seg)) {
    if (verbose) {
      print(paste("Phasing segment: ", s, "/", n.seg, sep=""))
    }
    n.snp <- nrow(h.snp.clust.p[[s]])
    new.h.snp.gt.p[[s]] <- h.snp.gt.p[[s]] 
    unfixed.h.snp.clust.p[[s]] <- h.snp.clust.p[[s]]
    
    allele.tab <- array(h.snp.clust.p[[s]], dim=c(1, n.snp, 5))
	  dnames = list()
	  dnames[[2]] = rownames(h.snp.gt.p[[s]])
	  dimnames(allele.tab) <- dnames

    chr <- h.snp.annot[[s]][["chr"]]
    dbSNP.annot <- h.snp.annot[[s]][["dbSNP"]]

    ix <- rep(FALSE, nrow(dbSNP.annot))
    ix[which(rowSums(is.na(dbSNP.annot)) == 0)] <- TRUE

    new.h.snp.gt.p[[s]][!ix, ] <- h.snp.gt.p[[s]][!ix, ]
    unfixed.h.snp.clust.p[[s]][!ix, ] <- h.snp.clust.p[[s]][!ix, ]

    n.annot.snp <- sum(ix)
    if (n.annot.snp < kMinSnps) {
      next
    }
    
    dbSNP.annot <- dbSNP.annot[ix, ]
    cols <- colnames(dbSNP.annot)
    dbSNP.annot <- cbind(dbSNP.annot, h.snp.annot[[s]][["pos"]][ix])
    colnames(dbSNP.annot) <- c(cols, "pos")
    
    allele.tab <- allele.tab[, ix, , drop=FALSE]
    phase.res <- CallBeagle(chr, allele.tab, dbSNP.annot,
                            platform, tmp.dir, plate.name,
                            phased.bgl.dir, verbose=verbose)
    
    if (phase.res[["fail"]] == FALSE) {
      phased.snp.clust.p <- phase.res[["phase.snp.clust.p"]][1, , ]
      unfixed.h.snp.clust.p[[s]][ix, c(1:4)] <- phased.snp.clust.p
      unfixed.h.snp.clust.p[[s]][ix, 5] <- 0
      
      ## add some uncertainty to phased results
      phased.snp.clust.p <- phased.snp.clust.p + 1e-2
      phased.snp.clust.p <- phased.snp.clust.p / rowSums(phased.snp.clust.p)

      res <- HmmResolveSwitchErrors(phased.snp.clust.p, allele.tab[1, , c(1:4)],
                                    verbose=verbose)
      
      new.h.snp.gt.p[[s]][ix, ] <- res[["phased.snp.clust.p"]]
      
      h.switch.ix[[s]] <- res[["switch.ix"]]
    }
  }
  if (verbose) {
    cat("\n")
  }
  
  return(list(h.snp.gt.p=new.h.snp.gt.p, unfixed.h.snp.clust.p=unfixed.h.snp.clust.p, h.switch.ix=h.switch.ix))
}

HmmResolveSwitchErrors <- function(phased.snp.clust.p, snp.clust.p,
                                   verbose=FALSE) {
  res <- BestViterbiPath(phased.snp.clust.p, snp.clust.p)
  snp.path <- res[["best.path"]]

  ix <- which(snp.path == 2) 
  tmp <- phased.snp.clust.p
  phased.snp.clust.p[ix, 2] <- tmp[ix, 4]
  phased.snp.clust.p[ix, 4] <- tmp[ix, 2]

  switch.ix <- which(diff(snp.path) != 0)

  n.switch <- length(switch.ix) 
  n.snp <- length(snp.path)

  if (verbose) {
    print(paste(n.switch, " switch errors detected over ", n.snp,
                " markers.", sep=""))
  }
  
  return(list(phased.snp.clust.p=phased.snp.clust.p,
              switch.ix=switch.ix))
}

BestViterbiPath <- function(phased.snp.clust.p, snp.clust.p) {
  kStartP <- c(0.5, 0.5)
  kTransP <- c(0.9, 0.1)

  rev.phased.snp.clust.p <- phased.snp.clust.p
  rev.phased.snp.clust.p[, 2] <- phased.snp.clust.p[, 4]
  rev.phased.snp.clust.p[, 4] <- phased.snp.clust.p[, 2]

  ph1.snp.ll <- LogAdd(log(snp.clust.p) + log(phased.snp.clust.p))
  ph2.snp.ll <- LogAdd(log(snp.clust.p) + log(rev.phased.snp.clust.p))
  prob.snp <- cbind(ph1.snp.ll, ph2.snp.ll)
  log.prob.snp <- prob.snp - LogAdd(prob.snp) 
  
  n.snp <- nrow(log.prob.snp)   
  log.state.prob <- matrix(NA, nrow=2, ncol=n.snp)
  
  which.state <- matrix(NA, nrow=2, ncol=n.snp)
  
  log.state.prob[1, 1] <- log.prob.snp[1, 1] + log(kStartP[1])
  log.state.prob[2, 1] <- log.prob.snp[1, 2] + log(kStartP[2])
    
  for (i in 2:n.snp) {
    prev.state <- c(log(kTransP[1]) + log.prob.snp[i-1, 1], 
                    log(kTransP[2]) + log.prob.snp[i-1, 2])
    
    which.state[1, i] <- which.max(prev.state)
    log.state.prob[1, i] <- log.prob.snp[i, 1] + max(prev.state)
    prev.state <- c(log(kTransP[2]) + log.prob.snp[i-1, 1], 
                    log(kTransP[1]) + log.prob.snp[i-1, 2])
    
    which.state[2, i] <- which.max(prev.state)
    log.state.prob[2, i] <- log.prob.snp[i, 2] + max(prev.state)
  }
  
  best.path <- rep(NA, n.snp)
  best.path[1] <- which.max(log.state.prob[, 1])
  
  best.path[n.snp] <- which.max(log.state.prob[, n.snp ])
  state <- best.path[n.snp]
  
  for (i in 1:(n.snp - 1)) {
    ix <- n.snp - i
    best.path[ix] <- which.state[state, ix] 
    state <- best.path[ix]
  }
  
  return(list(best.path=best.path, which.state=which.state,
              log.state.prob=log.state.prob))
}

