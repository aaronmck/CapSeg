JoinSmallSegsCapture = function(res, min.probes, verbose=FALSE )
{
   h.seg.dat <- res[["as.res"]][["h.seg.dat"]]
# h.capseg.d
# gh.wes.allele.d
# h.capseg.annot
# h.peobe.annot

   seg.dat <- res[["seg.dat"]]
   
   dat.types <- names(h.seg.dat)
   annot.types <- grep("annot", names(h.seg.dat), value=T)
   d.types <- setdiff(dat.types, annot.types)
   
   
   l = sapply(dat.types, function(n) length(h.seg.dat[[n]])) 
   if( !all(l[1] == l[2:length(l)])) stop ("There is different segmentations between Capseg probes / SNPs / annotations.")

   n.seg <- l[1]
   seg.chrs <- sapply(seq(n.seg), function(i) h.seg.dat[[annot.types[1]]][[i]][["chr"]])
   
   seg.n.probes <- foreach(x=d.types, .combine=cbind) %do% {
      unlist(lapply(h.seg.dat[[x]], ncol))
      }; colnames(seg.n.probes) <- d.types
   
   suff.probes.on.chr <- sapply(1:nrow(seg.n.probes), function(i) {
         chr.ix <- which(seg.chrs == seg.chrs[i])
         if (sum(seg.n.probes[chr.ix,"h.capseg.d"]) < min.probes ) {
            return(FALSE) # There aren't enough probes on the chr.  Force false.
         } else return(TRUE)
      })
   # ix <- which(rowSums(seg.n.probes) < min.probes | rowSums(seg.n.probes[,c("h.snp.d", "h.cn.d")]) < min.probes)
   ix <- which(seg.n.probes[,"h.capseg.d"] < min.probes & suff.probes.on.chr)
   
   IterateMergeSegs <- function(ix.1, ix.2) {

      out <- lapply(d.types, function(n) {
            d <- h.seg.dat[[n]]
            annot <- h.seg.dat[[gsub("\\.d$", ".annot", n)]]
            MergeTwoSegsExtreme(ix.1, ix.2, d, annot)
            }) 
      names(out) <- d.types
      return(out)
   }

   while (length(ix) > 0) {
      
      ix <- ix[1]
      if (verbose) {
         cat(paste(ix, ":", sep=""))
         # print(paste("left:", seg.chrs[ix-1], "center:", seg.chrs[ix], "right:", seg.chrs[ix+1]) )
      }
      if (ix == 1) {
         merge.res <- IterateMergeSegs(1, 2)
         seg.dat <- MergeSegTab(1, 2, seg.dat)
      } else if (ix == n.seg) {
         merge.res <- IterateMergeSegs(n.seg - 1, n.seg)
         seg.dat <- MergeSegTab(n.seg - 1, n.seg, seg.dat)
      } else {
         idxs = GetMergeableSegIndicesExtreme(seg.chrs, ix, verbose=verbose, h.seg.dat[d.types])
         merge.res <- IterateMergeSegs(idxs[1], idxs[2])
         seg.dat = MergeSegTab(idxs[1], idxs[2], seg.dat)
      }
      if (verbose) {
         cat(", ")
      }

      for (n in names(merge.res) ) { 
         # n = "h.snp.d"
         annot.n <- gsub("\\.d$", ".annot", n)
         h.seg.dat[[n]] <- merge.res[[n]][["h.d"]]
         h.seg.dat[[annot.n]] <- merge.res[[n]][["h.probe.annot"]]
      }

      l = sapply(dat.types, function(n) length(h.seg.dat[[n]])) 
      if( !all(l[1] == l[2:length(l)])) stop ("Small segments were joined and now probe segmentations are off.")
      
      n.seg <- length(h.seg.dat[[d.types[1]]])
      seg.chrs <- sapply(seq(n.seg), function(i) h.seg.dat[[annot.types[1]]][[i]][["chr"]])
      seg.n.probes <- foreach(x=d.types, .combine=cbind) %do% {
         unlist(lapply(h.seg.dat[[x]], ncol))
      }; colnames(seg.n.probes) <- d.types


      suff.probes.on.chr <- sapply(1:nrow(seg.n.probes), function(i) {
         chr.ix <- which(seg.chrs == seg.chrs[i])
         if (sum(seg.n.probes[chr.ix,"h.capseg.d"]) < min.probes ) {
            return(FALSE) # There aren't enough probes on the chr.  Force false.
         } else return(TRUE)
      })
      # ix <- which(rowSums(seg.n.probes[,c("h.snp.d", "h.cn.d")]) < min.probes | suff.probes.on.chr)
      ix <- which((seg.n.probes[,"h.capseg.d"] < min.probes) & suff.probes.on.chr)
   }

   if (verbose) {
      cat("\n")
   }
   
   res[["as.res"]][["h.seg.dat"]] <- h.seg.dat
   res[["seg.dat"]] <- seg.dat
   
   return(res)
}




GetSegLogEvCapture <- function(idx, h.seg.dat, theta, f.tau=NULL, f.approx.method="grid", verbose=FALSE) {
   # idx <- 498
   # f.tau <- f.tau[idx, ]

   GetTauSegLL <- function(x, d.cn, Theta) {
      tau <- x[1]
      return(CalcCaptureSegTauLogLik(d.cn, tau, Theta) )
   }

   GetAllelicSegLL <- function(x, d.allele, Theta) {
      alt <- d.allele["alt", ]
      ref <- d.allele["ref", ]
      f <- x[1]
      if (f < 0) f <- 1e-20 # This can happen when f is near zero and the hessian function drops a grid stupid
      res <- CalcCaptureSegAllelicLogLik(alt, ref, f, Theta )
      # print(paste("f=", f, "ll=", res))
      return( res)
   }
   
   if (is.null(f.tau)) {
      single.h.seg.dat <- lapply(h.seg.dat, "[", idx)
      f <- CaptureASModelFit ( single.h.seg.dat, verbose=FALSE )[["wes.f"]][, "f.hat"]
      tau <- CaptureCNModelFit ( single.h.seg.dat, verbose=FALSE )[["tau"]]
      f.tau <- cbind(f, tau)
   }

   f <- f.tau[1]
   tau <- f.tau[2]

   # Tau Evidence Calculations   
   if ( length(h.seg.dat[["h.capseg.d"]][[idx]]) != 0 && !is.na(tau)) {
      hess.mat.tau <- hessian(GetTauSegLL, c(tau), "Richardson", d.cn=h.seg.dat[["h.capseg.d"]][[idx]], Theta=theta) 
      curv.tau <- abs(det(hess.mat.tau / (2 * pi)))
      ll.tau <- GetTauSegLL(c(tau), h.seg.dat[["h.capseg.d"]][[idx]], Theta=theta); 
      log.ev.tau <- ll.tau + (log((curv.tau)^(-1 / 2))) - log(5) # prior is 1/5

   } else {
      ll.tau <- -Inf
      log.ev.tau <- 0   
   }

   # F Evidence Calculations
   if (length(h.seg.dat[["gh.wes.allele.d"]][[idx]]) != 0 && !is.na(f) ) {
      ll.f <- GetAllelicSegLL(c(f), h.seg.dat[["gh.wes.allele.d"]][[idx]], Theta=theta); 
      if (grepl("laplace", f.approx.method, ignore.case=TRUE)) {
         # Laplace approximation to the normalizing constant
         hess.mat.f <- hessian(GetAllelicSegLL, c(f), "Richardson", d.allele=h.seg.dat[["gh.wes.allele.d"]][[idx]], Theta=theta) 
         curv.f <- abs(det(hess.mat.f / (2 * pi)))
         log.ev.f <- ll.f + (log((curv.f)^(-1 / 2))) - log(0.5) # Prior is 2 because 0 <= F <= 0.5   

      } else if (grepl("grid", f.approx.method, ignore.case=TRUE)) {
         # Grid approximation to the normalizing constant
         log.ev.f <- GridEstimateOfFNormalizingConstant(alt=h.seg.dat[["gh.wes.allele.d"]][[idx]]["alt", ], ref=h.seg.dat[["gh.wes.allele.d"]][[idx]]["ref", ], Theta=theta, dx=1e-2)
      } else stop(paste("F.approx.method:", f.approx.method, "not supported/recognized."))

   } else {
      ll.f <- -Inf
      log.ev.f <- 0   
   }
   
   return(list(log.ev=log.ev.tau + log.ev.f, ll=ll.tau + ll.f, log.ev.tau=log.ev.tau, log.ev.f=log.ev.f, ll.tau=ll.tau, ll.f=ll.f))
   
}

CalculateH1EvidenceCapture <- function(ix.1, ix.2, h.seg.dat, theta, f.tau, verbose=FALSE) {      
   
   mrg.capseg.d <- cbind(h.seg.dat[["h.capseg.d"]][[ix.1]], h.seg.dat[["h.capseg.d"]][[ix.2]])
   mrg.allele.d <- cbind(h.seg.dat[["gh.wes.allele.d"]][[ix.1]], h.seg.dat[["gh.wes.allele.d"]][[ix.2]])
   
   mrg.h.seg.dat <- h.seg.dat
   n.segs <- length(h.seg.dat[[1]] )
   # Merge the data for segs ix.1 and ix.2 into ix.1 to get H1 Evidence.  
   # Note that this isn't comprehensive merge and that there will be other lists (eg, h.snp.annot) that will be off, but these are not returned from this function
   if (ix.1 == 1) {
      mrg.h.seg.dat[["h.capseg.d"]][[1]] <- mrg.capseg.d
      mrg.h.seg.dat[["h.capseg.d"]][2:(n.segs-1)] <- h.seg.dat[["h.capseg.d"]][3:n.segs]
      mrg.h.seg.dat[["h.capseg.d"]][[n.segs]] <- NULL

      mrg.h.seg.dat[["gh.wes.allele.d"]][[1]] <- mrg.allele.d
      mrg.h.seg.dat[["gh.wes.allele.d"]][2:(n.segs-1)] <- h.seg.dat[["gh.wes.allele.d"]][3:n.segs]
      mrg.h.seg.dat[["gh.wes.allele.d"]][[n.segs]] <- NULL

   } else if (ix.1 == n.segs - 1) {
      mrg.h.seg.dat[["h.capseg.d"]][[n.segs-1]] <- mrg.capseg.d
      mrg.h.seg.dat[["h.capseg.d"]][[n.segs]] <- NULL

      mrg.h.seg.dat[["gh.wes.allele.d"]][[n.segs-1]] <- mrg.allele.d
      mrg.h.seg.dat[["gh.wes.allele.d"]][[n.segs]] <- NULL
   } else {
      mrg.h.seg.dat[["h.capseg.d"]][[ix.1]] <- mrg.capseg.d
      mrg.h.seg.dat[["h.capseg.d"]][ix.2:(n.segs-1)] <- h.seg.dat[["h.capseg.d"]][(ix.2+1):(n.segs)]
      mrg.h.seg.dat[["h.capseg.d"]][[n.segs]] <- NULL

      mrg.h.seg.dat[["gh.wes.allele.d"]][[ix.1]] <- mrg.allele.d
      mrg.h.seg.dat[["gh.wes.allele.d"]][ix.2:(n.segs-1)] <- h.seg.dat[["gh.wes.allele.d"]][(ix.2+1):(n.segs)]
      mrg.h.seg.dat[["gh.wes.allele.d"]][[n.segs]] <- NULL
   }

   res <- GetSegLogEvCapture(ix.1, mrg.h.seg.dat, theta, f.tau=NULL, verbose=verbose)
   return(list(tot.log.ev=res[["log.ev"]], tau=res[["log.ev.tau"]], f=res[["log.ev.f"]]))
}

CalculateH0EvidenceCapture <- function(ix.1, ix.2, h.seg.dat, theta, f.tau, verbose=FALSE) {

#   ev1 <- GetSegLogEvCapture(ix.1, h.seg.dat, theta, f.tau[ix.1,], verbose=verbose)
#   ev2 <- GetSegLogEvCapture(ix.2, h.seg.dat, theta, f.tau[ix.2,], verbose=verbose)
   
   ev1 <- GetSegLogEvCapture(ix.1, h.seg.dat, theta, f.tau=NULL, verbose=verbose)
   ev2 <- GetSegLogEvCapture(ix.2, h.seg.dat, theta, f.tau=NULL, verbose=verbose)

   log.ev1 <- ev1[["log.ev"]]
   log.ev2 <- ev2[["log.ev"]]
   
   return(log.ev1 + log.ev2)
}

CalcSegMergeProbCapture <- function(ix.1, ix.2, h.seg.dat, theta, f.tau, verbose=FALSE) {
   ## H0:  two segs are seperate
   ## H1:  two segs are same

   h0.log.ev <- CalculateH0EvidenceCapture(ix.1, ix.2, h.seg.dat, theta, f.tau, verbose=verbose)
   h1.log.ev <- CalculateH1EvidenceCapture(ix.1, ix.2, h.seg.dat, theta, f.tau, verbose=verbose)[["tot.log.ev"]]
   h1.probs <- CalculateH1Probs( h0.log.ev, h1.log.ev )
   
   return(data.frame(merge.prob=h1.probs[["h1.prob"]], log10.merge.prob=h1.probs[["log10.h1.prob"]], h0.log.ev=h0.log.ev, h1.log.ev=h1.log.ev))

}

JoinCloseSegsCapture <- function(h.seg.dat, capture.em.fit, merge.thresh=0.5, verbose=FALSE) 
{
   # h.seg.dat <- tot.res[["as.res"]][["h.seg.dat"]]
   # capture.em.fit <- tot.res[["capture.em.fit"]]
   # merge.thresh <- seg.merge.thresh
   # verbose=TRUE

   theta <- capture.em.fit[["Theta"]]
   f.tau <- cbind(f=capture.em.fit[["wes.f"]][,"f.hat"], tau=capture.em.fit[["delta.tau"]][["tau"]])
   if (verbose) {
      print("Merging close segments on capture platform...")
   }
   
   n.seg <- length(h.seg.dat$h.capseg.d)
   print(paste("n.seg:", n.seg))
   col.names <- c('merge.prob', 'log10.merge.prob', 'h0.log.ev.cap', 'h1.log.ev.cap')
   merge.prob <- matrix(NA, nrow=(n.seg - 1), ncol=length(col.names)) 
   colnames(merge.prob) <- col.names
   merged.loci <- data.frame()
   
   seg.chrs <- sapply(h.seg.dat[["h.capseg.annot"]], "[[", "chr")
   d.types <- c("h.capseg.d")
   seg.n.probes <- cbind(foreach(dt=d.types, .combine=cbind) %do% sapply(h.seg.dat[[dt]], ncol))
   colnames(seg.n.probes) <- gsub("(^h\\.)|(\\.d$)", "", d.types)
   
   while (TRUE) {

      na.ix <- which(is.na(merge.prob[, "merge.prob"]))
      delta.chr.ix <- (diff(seg.chrs) != 0)
      
      for (i in 1:length(na.ix)) {
         # i <- 1

         ## never merge across different chrs
         if (delta.chr.ix[na.ix[i]] == FALSE) {
            res <- CalcSegMergeProbCapture(na.ix[i], na.ix[i] + 1, h.seg.dat, theta, f.tau, verbose=verbose)
            if ( !is.finite(res[1, "merge.prob"]) ) {
               merge.prob[na.ix[i], ] <- c(-1, NaN, res[,"h0.log.ev"], res[, "h1.log.ev"] )
               if (verbose) {
                  print(paste( "CalcSegMergeProb error: ", na.ix[i], sep="" ))
               }
            } else {
               merge.prob[na.ix[i], ] <- as.matrix(res[1,])
            }
            if (verbose) cat(".")
         } else {
            merge.prob[na.ix[i], ] <- c(-1, rep(NaN, length(col.names)-1) ) 
         }
      }
      if (verbose) {
         cat("\n")
      }
      ## merge pair with highest prob or stop if none
      if (max(merge.prob[, "merge.prob"]) < merge.thresh) {
         break
      }
      
      m.ix <- which.max(merge.prob[, "merge.prob"]) 
      if (verbose) {
         cat(paste("merging #", m.ix, ", prob= ", signif(merge.prob[m.ix, "merge.prob"], 4), ", log10 prob= ", signif(merge.prob[m.ix, "log10.merge.prob"], 4),  ",  h0.log.ev.cap= ", signif(merge.prob[m.ix, "h0.log.ev.cap"], 4), ",  h1.log.ev.cap= ", signif(merge.prob[m.ix, "h1.log.ev.cap"], 4), "\n", sep=""))
      }

      annot <- h.seg.dat[["h.capseg.annot"]][[m.ix]]
      chr <- h.seg.dat[["h.capseg.annot"]][[m.ix]][["chr"]]
      chr <- ifelse(!all(complete.cases(chr)), max(chr), chr[1])
      pos <- h.seg.dat[["h.capseg.annot"]][[m.ix]][["pos"]][length(h.seg.dat[["h.capseg.annot"]][[m.ix]][["pos"]])]
      merged.loci <- rbind(merged.loci, c(chr, pos, merge.prob[m.ix, ]) )
      
      merge.cap.res <- MergeTwoSegsExtreme(m.ix, m.ix + 1, h.seg.dat[["h.capseg.d"]], h.seg.dat[["h.capseg.annot"]])
      h.seg.dat[["h.capseg.d"]] <- merge.cap.res[["h.d"]]
      h.seg.dat[["h.capseg.annot"]] <- merge.cap.res[["h.probe.annot"]]   

      merge.gh.res <- MergeTwoSegsExtreme(m.ix, m.ix + 1, h.seg.dat[["gh.wes.allele.d"]], h.seg.dat[["gh.wes.allele.annot"]])
      h.seg.dat[["gh.wes.allele.d"]] <- merge.gh.res[["h.d"]]
      h.seg.dat[["gh.wes.allele.annot"]] <- merge.gh.res[["h.probe.annot"]]   
      
      
      ## update merging probs
      nmp <- matrix(NA, nrow=0, ncol=ncol(merge.prob))
      
      if (m.ix > 1) {
         nmp <- merge.prob[c(1:(m.ix - 1)), , drop=FALSE]
         nmp[m.ix - 1, ] <- NA
      }
      
      if (m.ix < nrow(merge.prob)) {
         nmp <- rbind(nmp, merge.prob[c((m.ix + 1):nrow(merge.prob)), ])
         nmp[m.ix, ] <- NA
      }
      
      merge.prob <- nmp
      n.seg <- length(h.seg.dat$h.capseg.d)
      seg.chrs <- sapply(h.seg.dat[["h.capseg.annot"]], "[[", "chr")
      seg.n.probes <- cbind(foreach(dt=d.types, .combine=cbind) %do% sapply(h.seg.dat[[dt]], ncol))
      colnames(seg.n.probes) <- gsub("(^h\\.)|(\\.d$)", "", d.types)
      
   } ## End while(TRUE)   
   
   if (length(merged.loci) > 0) {
      colnames(merged.loci) <- c("chr", "pos", "prob", "log10_prob", "h0.log.ev.cap", "h1.log.ev.cap")
      colnames(merge.prob) <- c("merge.prob", "log10.merge.prob", "h0.log.ev.cap", "h1.log.ev.cap")
   }
   
   list(h.seg.dat=h.seg.dat, merged.loci=merged.loci, final.merge.prob=merge.prob)
}

