PhaseSnps.chunked <- function(h.snp.clust.p, h.snp.gt.p, h.snp.annot,
                      platform, tmp.dir, plate.name, phased.bgl.dir, verbose=FALSE) {
  kMinSnps <- 100

  new.h.snp.gt.p <- list()
  unfixed.h.snp.clust.p <- list()
  h.switch.ix <- list()

  n.segs <- length(h.snp.clust.p)
  for (s in seq_len(n.segs)) {
    if (verbose) {
      print(paste("Phasing segment: ", s, "/", n.segs, sep=""))
    }
    cur.snp.clust.p <- h.snp.clust.p[[s]]
    cur.snp.gt.p <- h.snp.gt.p[[s]]
    cur.snp.annot <- h.snp.annot[[s]]      
    cur.chr <- cur.snp.annot[["chr"]]
    cur.snp.pos <- cur.snp.annot[["pos"]]
    cur.dbSNP <- cur.snp.annot[["dbSNP"]]

    ix <- rep(FALSE, nrow(cur.dbSNP))
    ix[which(rowSums(is.na(cur.dbSNP)) == 0)] <- TRUE
    n.annot.snp <- sum(ix)
    if (n.annot.snp < kMinSnps) {
      new.h.snp.gt.p[[s]] <- h.snp.gt.p[[s]]
      unfixed.h.snp.clust.p[[s]] <- h.snp.clust.p[[s]]
      next
    }
    
    phase.results <- ProcessSegPhasing(cur.snp.clust.p, cur.snp.gt.p,
                                       cur.chr, cur.snp.pos,
                                       cur.dbSNP, platform, tmp.dir,
                                       plate.name, phased.bgl.dir,
                                       ix, kMinSnps, verbose=verbose)
    
    new.h.snp.gt.p[[s]] <- phase.results[["h.snp.gt.p"]]
    unfixed.h.snp.clust.p[[s]] <- phase.results[["unfixed.h.snp.clust.p"]]
    h.switch.ix[[s]] <- phase.results[["h.switch.ix"]]
  }
  if (verbose) {
    cat("\n")
  }
  
  return(list(h.snp.gt.p=new.h.snp.gt.p,
              unfixed.h.snp.clust.p=unfixed.h.snp.clust.p,
              h.switch.ix=h.switch.ix))
}

ProcessSegPhasing <- function(h.snp.clust.p, h.snp.gt.p, chr,
                              snp.annot.pos, dbSNP.annot,
                              platform, tmp.dir, plate.name,
                              phased.bgl.dir, ix, min.snps,
                              verbose=FALSE) {
  new.h.snp.gt.p <- h.snp.gt.p 
  unfixed.h.snp.clust.p <- h.snp.clust.p

  chunk.factor <- GetChunkFactor(snp.annot.pos, min.snps)
  chunked.snp.clust.p.rows <- split(rownames(h.snp.clust.p), chunk.factor)
  chunked.snp.gt.p.rows <- split(rownames(h.snp.gt.p), chunk.factor)
  chunked.snp.pos <- split(snp.annot.pos, chunk.factor)
  chunked.dbSNP.rows <- split(rownames(dbSNP.annot), chunk.factor)
  chunked.ix <- split(ix, chunk.factor)
  
  beagle.failure <- FALSE
  phased.snp.clust.p <- NULL
  allele.tab <- NULL
  ## FIXME: factor out this for loop into a separate function - also
  ## foreach() it for potential parallelization
  for (i in seq_along(levels(chunk.factor))) {
    snp.clust.p.chunk <- h.snp.clust.p[chunked.snp.clust.p.rows[[i]], ]
    snp.gt.p.chunk <- h.snp.gt.p[chunked.snp.gt.p.rows[[i]], ]
    snp.pos.chunk <- chunked.snp.pos[[i]]
    snp.dbSNP.chunk <- dbSNP.annot[chunked.dbSNP.rows[[i]], ]
    cur.ix <- chunked.ix[[i]]
    
    cur.allele.tab <- array(snp.clust.p.chunk, dim=c(1, nrow(snp.clust.p.chunk), 5))
    dimnames(cur.allele.tab)[[2]] <- rownames(snp.gt.p.chunk)
    
    snp.dbSNP.chunk <- snp.dbSNP.chunk[cur.ix, ]
    cols <- colnames(snp.dbSNP.chunk)
    snp.dbSNP.chunk <- cbind(snp.dbSNP.chunk, snp.pos.chunk[cur.ix])
    colnames(snp.dbSNP.chunk) <- c(cols, "pos")
    
    cur.allele.tab <- cur.allele.tab[, cur.ix, , drop=FALSE]
    phase.res <- CallBeagle(chr, cur.allele.tab, snp.dbSNP.chunk,
                            platform, tmp.dir, plate.name,
                            phased.bgl.dir, verbose=verbose)
    if (phase.res[["fail"]]) {
      beagle.failure <- TRUE
    }
    ## FIXME: The raw rbinds should be replaced by a pre-allocated
    ## matrix. I doubt it matters *that* much, but ....
    if (is.null(phased.snp.clust.p)) {
      phased.snp.clust.p <- phase.res[["phase.snp.clust.p"]][1, , ]
    } else {
      phased.snp.clust.p <- rbind(phased.snp.clust.p,
                                  phase.res[["phase.snp.clust.p"]][1, , ])
    }
    if (is.null(allele.tab)) {
      allele.tab <- cur.allele.tab[1, , ]
    } else {
      allele.tab <- rbind(allele.tab, cur.allele.tab[1, , ])
    }
  }
  
  if (beagle.failure) {
    h.switch.ix <- numeric()
  } else {
    unfixed.h.snp.clust.p[ix, c(1:4)] <- phased.snp.clust.p
    unfixed.h.snp.clust.p[ix, 5] <- 0
    
    ## add some uncertainty to phased results
    phased.snp.clust.p <- phased.snp.clust.p + 1e-2
    phased.snp.clust.p <- phased.snp.clust.p / rowSums(phased.snp.clust.p)
    
    res <- HmmResolveSwitchErrors(phased.snp.clust.p, allele.tab[, c(1:4)],
                                  verbose=verbose)
    
    new.h.snp.gt.p[ix, ] <- res[["phased.snp.clust.p"]]
    
    h.switch.ix <- res[["switch.ix"]]
  }
  
  return(list(h.snp.gt.p=new.h.snp.gt.p,
              unfixed.h.snp.clust.p=unfixed.h.snp.clust.p,
              h.switch.ix=h.switch.ix))
}  

GetChunkFactor <- function(seg.pos, min.snps) {
  ## The idea is to split the segment into chunk.size blocks.
  ## This won't perfectly do this - rather it assumes a fairly
  ## even distribution. Should be tweaked, but it works
  ## close enough for now
  kBeagleThreshold <- 5000000

  seg.size <- GetSegmentSize(seg.pos)
  ## Don't use integer division as it'll not necessarily ceiling
  ## the operation - which we want to guarantee at least 1 chunk
  n.chunks <- ceiling(seg.size / kBeagleThreshold)
  if (n.chunks > 1) {
    out.factor <- cut(seg.pos, breaks=n.chunks)
  } else {
    out.factor <- as.factor(rep(1, length(seg.pos)))
  }

  ## Make sure that every chunk has at least min.snps snps in
  ## it
  factor.table <- table(out.factor)
  ## first make sure that the first chunk is ok - keep combining the
  ## first segment with the next one until either it meets the minimum
  ## or we're down to a single chunk
  while (factor.table[1] < min.snps) {
    if (length(factor.table) == 1) {
      break
    }
    levels(out.factor)[2] <- levels(out.factor)[1]
    factor.table <- table(out.factor)
  }
  ## Now verify the rest of the chunks using the opposite logic. Move
  ## through the levels and every time one does not meet the minimum, combine
  ## it with the previously known valid level.
  last.valid <- 1
  factor.table <- table(out.factor)
  while (last.valid < length(factor.table)) {
    next.level <- last.valid + 1
    if (factor.table[next.level] < min.snps) {
      levels(out.factor)[next.level] <- levels(out.factor)[last.valid]
      factor.table <- table(out.factor)
      next
    } else {
      last.valid <- next.level
    }
  }
  
  return(out.factor)
}
  
GetSegmentSize <- function(seg.pos) {
  return (seg.pos[length(seg.pos)] - seg.pos[1])
}

