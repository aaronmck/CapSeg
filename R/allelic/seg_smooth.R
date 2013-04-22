## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

JoinSmallSegs <- function(res, min.probes, verbose=FALSE) {
	
	# res = jss.res ; min.probes = min.seg.size
	 
	h.seg.dat <- res[["as.res"]][["h.seg.dat"]]
	seg.dat <- res[["seg.dat"]]
	
	dat.types <- names(h.seg.dat)
	annot.types <- grep("annot", names(h.seg.dat), value=T)
	d.types <- setdiff(dat.types, c(annot.types, "h.snp.gt.p"))
	
	
	l = sapply(dat.types, function(n) length(h.seg.dat[[n]])) 
	if( !all(l[1] == l[2:length(l)])) stop ("There is different segmentations between SNP, CN, or Capseg probes.")

	n.seg <- l[1]
	seg.chrs <- sapply(seq(n.seg), function(i) h.seg.dat[[annot.types[1]]][[i]][["chr"]])
	
	seg.n.probes <- foreach(x=d.types, .combine=cbind) %do% {
		unlist(lapply(h.seg.dat[[x]], ncol))
		}; colnames(seg.n.probes) <- d.types
	
	suff.probes.on.chr <- sapply(1:nrow(seg.n.probes), function(i) {
			chr.ix <- which(seg.chrs == seg.chrs[i])
			if (sum(seg.n.probes[chr.ix,c("h.snp.d", "h.cn.d")]) < min.probes ) {
				return(FALSE) # There aren't enough probes on the chr.  Force false.
			} else return(TRUE)
		})
	# ix <- which(rowSums(seg.n.probes) < min.probes | rowSums(seg.n.probes[,c("h.snp.d", "h.cn.d")]) < min.probes)
	ix <- which(rowSums(seg.n.probes[,c("h.snp.d", "h.cn.d")]) < min.probes & suff.probes.on.chr)
	
	IterateMergeSegs <- function(ix.1, ix.2) {

		out <- lapply(d.types, function(n) {
				d <- h.seg.dat[[n]]
				annot <- h.seg.dat[[gsub("\\.d$", ".annot", n)]]
				if (n == "h.snp.d"){
					h.snp.gt.p <- h.seg.dat[['h.snp.gt.p']]
					MergeTwoSegsExtreme(ix.1, ix.2, d, annot, h.snp.gt.p)	
				} else {
					MergeTwoSegsExtreme(ix.1, ix.2, d, annot)
				}
				}) 
		names(out) <- d.types
		return(out)
	}

	while (length(ix) > 0) {
		# browser()
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
			idxs = GetMergeableSegIndicesExtreme(seg.chrs, ix, h.seg.dat[d.types])
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
			if (n == "h.snp.d") {
				h.seg.dat[["h.snp.gt.p"]] <- merge.res[[n]][["h.snp.gt.p"]]
			}
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
			if (sum(seg.n.probes[chr.ix,c("h.snp.d", "h.cn.d")]) < min.probes ) {
				return(FALSE) # There aren't enough probes on the chr.  Force false.
			} else return(TRUE)
		})
		# ix <- which(rowSums(seg.n.probes[,c("h.snp.d", "h.cn.d")]) < min.probes | suff.probes.on.chr)
		ix <- which((rowSums(seg.n.probes[,c("h.snp.d", "h.cn.d")]) < min.probes) & suff.probes.on.chr)
	}

	if (verbose) {
		cat("\n")
	}
	
	res[["as.res"]][["h.seg.dat"]] <- h.seg.dat
	res[["seg.dat"]] <- seg.dat
	
	return(res)
	
	 
}


MergeSegTab = function(ix.1, ix.2, seg.tab) {

	a.cols = c("Chromosome", "Start.bp", "End.bp")
	if (!all(a.cols %in% colnames(seg.tab))) stop(paste("seg.tab does not have correct column names. I need: ", paste(a.cols, collapse=", "), " You have: ", colnames(seg.tab)))
	unmergeable.cols = setdiff(names(seg.tab), a.cols)
	if (length(unmergeable.cols) > 0 ) warning(paste("Dropping columns: ", paste(unmergeable.cols, collapse=", ")))
	if (ix.2 - ix.1 != 1) stop (paste("Cant merge non-adjacent segs: ", ix.1, ix.2))
	if (seg.tab$Chromosome[ix.1] != seg.tab$Chromosome[ix.2]) stop (paste("Cant merge segs on different chromosomes: ", seg.tab$Chromosome[ix.1], seg.tab$Chromosome[ix.2]))
	
	if (nrow(seg.tab) == 2) {
		out = data.frame(Chromosome = seg.tab$Chromosome[ix.1], Start.bp = seg.tab$Start.bp[ix.1], End.bp = seg.tab$End.bp[ix.2])
	} else if (ix.1 == 1) { 
		mrg = data.frame(Chromosome = seg.tab$Chromosome[ix.1], Start.bp = seg.tab$Start.bp[ix.1], End.bp = seg.tab$End.bp[ix.2])
		base = seg.tab[(ix.2 + 1):nrow(seg.tab), a.cols]
		out = rbind(mrg, base)
	} else if (ix.2 == nrow(seg.tab)) {
		head = seg.tab[1:( ix.1 - 1), a.cols]
		mrg = data.frame(Chromosome = seg.tab$Chromosome[ix.1], Start.bp = seg.tab$Start.bp[ix.1], End.bp = seg.tab$End.bp[ix.2])
		out = rbind(head, mrg)
	} else {
		head = seg.tab[1:( ix.1 - 1), a.cols]
		mrg = data.frame(Chromosome = seg.tab$Chromosome[ix.1], Start.bp = seg.tab$Start.bp[ix.1], End.bp = seg.tab$End.bp[ix.2])
		base = seg.tab[(ix.2 + 1):nrow(seg.tab), a.cols]
		out = rbind(head, mrg, base)
	}
	return(out)
}


GetMergeableSegIndicesExtreme = function(seg.chrs, ix, ...) {

	dat = list(...)[[1]]
	# dat = list(h.seg.dat[d.types])[[1]]
	if( sum(unlist(lapply(dat, function(x) ncol(x[[ix]]) ))) > 0 ) {
		smean <- median(unlist(lapply(dat, function(x) colSums(x[[ix]]) )))
		m1 <- median(unlist(lapply(dat, function(x) colSums(x[[ix-1]]) )))
		m2 <- median(unlist(lapply(dat, function(x) colSums(x[[ix+1]]) )))
	} else { # force a left merge
		smean = 2
		m1 = 2
		m2 =1
	}
	if (is.na(m2)) m2 <- m1 + smean + 1 # Don't merge to an empty segment first.
	
	d1 <- c(abs(smean - m1), abs(smean - m2))
	
	if ((seg.chrs[ix - 1] == seg.chrs[ix]) && (seg.chrs[ix] == seg.chrs[ix + 1])) {
		## seg before
		if (d1[1] < d1[2]) {
			ix.1 <- ix - 1
			ix.2 <- ix
			if (verbose) {
				cat("1")
			}
		} else {
			ix.1 <- ix
			ix.2 <- ix + 1
			if (verbose) {
				cat("2")
			}
		}
	} else {
		if (seg.chrs[ix - 1] == seg.chrs[ix]) {
			ix.1 <- ix - 1
			ix.2 <- ix
			if (verbose) {
				cat("c1")
			}
		}
		if (seg.chrs[ix] == seg.chrs[ix + 1]) {
			ix.1 <- ix
			ix.2 <- ix + 1
			if (verbose) {
				cat("c2")
			}
		}
	}
	return(c(ix.1, ix.2))
}

GetMergeableSegIndices = function(seg.chrs, ix, h.d, h.probe.annot) {
	smean <- median(colSums(h.d[[ix]]))
	m1 <- median(colSums(h.d[[ix - 1]]))
	m2 <- median(colSums(h.d[[ix + 1]]))
	
	d1 <- c(abs(smean - m1), abs(smean - m2))
	
	if ((seg.chrs[ix - 1] == seg.chrs[ix]) && (seg.chrs[ix] == seg.chrs[ix + 1])) {
		## seg before
		if (d1[1] < d1[2]) {
			ix.1 <- ix - 1
			ix.2 <- ix
			if (verbose) {
				cat("1")
			}
		} else {
			ix.1 <- ix
			ix.2 <- ix + 1
			if (verbose) {
				cat("2")
			}
		}
	} else {
		if (seg.chrs[ix - 1] == seg.chrs[ix]) {
			ix.1 <- ix - 1
			ix.2 <- ix
			if (verbose) {
				cat("c1")
			}
		}
		if (seg.chrs[ix] == seg.chrs[ix + 1]) {
			ix.1 <- ix
			ix.2 <- ix + 1
			if (verbose) {
				cat("c2")
			}
		}
	}
	return(c(ix.1, ix.2))
}

GetInteriorMergeResExtreme <- function(seg.chrs, ix, h.d, h.probe.annot, verbose=verbose, ...) {
#	seg.chrs = seg.chrs.capseg
#	ix = 4
#	h.d = h.capseg.d
#	h.probe.annot = h.capseg.annot
	
	## compare medians to before/after
	idxs = GetMergeableSegIndices(seg.chrs, ix, h.d, h.probe.annot) 
	return(MergeTwoSegsExtreme(idxs[1], idxs[2], h.d, h.probe.annot, ...))
}



MergeTwoSegsExtreme <- function(ix.1, ix.2, h.d, h.probe.annot, ...) {
	## ... is meta data that need not be shared between all h.* data
	## For example, snp probes contain h.snp.gt.p information, that CN probes
	## do not have.
	## All metadata must be lists, with each corresponding to a segments and
	## have probes on the rows since they will be merged with rbind.
	if (h.probe.annot[[ix.1]][["chr"]] != h.probe.annot[[ix.2]][["chr"]]) {
		stop("Attempt to merge segs from different chromosomes!")
	}
	meta = list(...)
	if (length(meta) > 0) {
		meta.names = as.list(substitute({...})[-1])
	}					
	n.seg <- length(h.d)
	
	new.h.d <- vector(mode="list", length=n.seg - 1)
	new.h.probe.annot <- vector(mode="list", length=n.seg - 1)
	if (length(meta) > 0 ) {
		new.meta = lapply(1:length(meta), function(i) vector(mode="list", length=n.seg - 1))
	}
	
	mrg.h.d <- cbind(h.d[[ix.1]], h.d[[ix.2]])
	if (length(meta) > 0 ) {
		mrg.meta = lapply(meta, function(elem) rbind(elem[[ix.1]], elem[[ix.2]]))
	}
	
	mrg.h.probe.annot <- h.probe.annot[[ix.1]]
	for (name in names(mrg.h.probe.annot)) { 
		if (name == "chr") {
			mrg.h.probe.annot[[name]] = h.probe.annot[[ix.1]][[name]]
		} else {
			if (is.matrix(mrg.h.probe.annot[[name]]) | is.data.frame(mrg.h.probe.annot[[name]])) {
				mrg.h.probe.annot[[name]] = rbind(h.probe.annot[[ix.1]][[name]], h.probe.annot[[ix.2]][[name]])
			} else {
				mrg.h.probe.annot[[name]] = c(h.probe.annot[[ix.1]][[name]], h.probe.annot[[ix.2]][[name]])
			}
		}
	}
	
	
	
	if (ix.1 > 1) {
		ix <- c(1:(ix.1 - 1))
		new.h.d[ix] <- h.d[ix]
		if (length(meta) > 0 ) {
			new.meta = lapply(1:length(new.meta), function(i) {
						new.meta[[i]][ix] <- meta[[i]][ix]
						return(new.meta[[i]])
					})
		}
		new.h.probe.annot[ix] <- h.probe.annot[ix]
	}
	
	new.ix <- ix.1
	new.h.d[[new.ix]] <- mrg.h.d
	if (length(meta) > 0 ) {
		new.meta = lapply(1:length(new.meta), function(i) {
					new.meta[[i]][[new.ix]] = mrg.meta[[i]]
					return(new.meta[[i]])
				})
#							new.h.snp.gt.p[[new.ix]] <- mrg.h.snp.gt.p
	}
	new.h.probe.annot[[new.ix]] <- mrg.h.probe.annot
	
	if (ix.2 < n.seg) {
		ix <- c((ix.2 + 1):n.seg) 
		new.h.d[ix - 1] <- h.d[ix] 
		if (length(meta) > 0 ) {
			new.meta = lapply(1:length(new.meta), function(i) {
						new.meta[[i]][ix - 1] = meta[[i]][ix]
						return(new.meta[[i]])
					})
#							new.h.snp.gt.p[ix - 1] <- h.snp.gt.p[ix]
		}
		new.h.probe.annot[ix - 1] <- h.probe.annot[ix] 
	}
	if (length(meta) > 0 ) {
		names(new.meta) <- meta.names
		out = c(list(h.d=new.h.d, h.probe.annot=new.h.probe.annot), new.meta)
	} else {
		out = list(h.d=new.h.d, h.probe.annot=new.h.probe.annot)
	}
	return(out)
}


JoinCloseSegs <- function(h.seg.dat, theta, force.diploid, out.p, merge.thresh=0.5, verbose=FALSE) {
	
	h.seg.dat <- iams.res[["as.res"]][["h.seg.dat"]]
	theta <- iams.res[["array.em.fit"]][["theta"]]

	if (verbose) {
		print("Merging close segments...")
	}
	
	n.seg <- length(h.seg.dat[[1]])
	merge.prob <- matrix(NA, nrow=(n.seg - 1), ncol=8) 
	colnames(merge.prob) <- c('merge.prob', 'log10.merge.prob', 'h0.log.ev.cn', 'h0.log.ev.snp', 'h1.log.ev.cn', 'h1.log.ev.snp', 'merge.prob.cn', 'merge.prob.snp')
	merged.loci <- data.frame()
	
	seg.chrs <- sapply(h.seg.dat[["h.snp.annot"]], "[[", "chr")
	# seg.n.snps <- sapply(h.d[["snp"]], ncol)
	d.types <- c("h.snp.d", "h.cn.d")
	seg.n.probes <- foreach(dt=d.types, .combine=cbind) %do% sapply(h.seg.dat[[dt]], ncol)
	colnames(seg.n.probes) <- gsub("(^h\\.)|(\\.d$)", "", d.types)
	
	WhileFunc <- function() {
		while (TRUE) {

		na.ix <- which(is.na(merge.prob[, "merge.prob"]))
		delta.chr.ix <- (diff(seg.chrs) != 0)
		
		# for (i in 1:length(na.ix)) {
		for (i in 244:246) {
			# i <- 245

			## never merge across different chrs
			if (delta.chr.ix[na.ix[i]] == FALSE) {
				# res <- CalcSegMergeProb(na.ix[i], na.ix[i] + 1, h.snp.d, h.snp.gt.p, theta, force.diploid, out.p, verbose=verbose)
				# res <- CalcSegMergeProbExtreme(na.ix[i], na.ix[i] + 1, h.d, h.snp.gt.p, theta, force.diploid, out.p, verbose=verbose)
				res <- CalcSegMergeProbExtreme(na.ix[i], na.ix[i] + 1, h.seg.dat, theta, force.diploid, out.p, verbose=verbose)
				if (is.finite(res[1, "merge.prob"]) == FALSE) {
					merge.prob[na.ix[i], ] <- c(-1, NaN, res[,"h0.log.ev.cn"], res[, "h0.log.ev.snp"], res[, "h1.log.ev.cn"], 
							res[, "h1.log.ev.snp"], res[, "merge.prob.cn"], res[, "merge.prob.snp"])
					if (verbose) {
						print(paste( "CalcSegMergeProb error: ", na.ix[i], sep="" ))
					}
				} else {
					merge.prob[na.ix[i], ] <- as.matrix(res[1,])
				}
				if (verbose) {
					cat(".")
				}
			} else {
				merge.prob[na.ix[i], ] <- c(-1, NaN, NaN, NaN, NaN, NaN, NaN, NaN) 
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
			cat(paste("merging #", m.ix, ", prob= ", round(merge.prob[m.ix, 1], 4), ", log10 prob= ", round(merge.prob[m.ix, 2], 4),  
							",  h0.log.ev.cn= ", round(merge.prob[m.ix, 3], 4), ",  h0.log.ev.snp= ", round(merge.prob[m.ix, 4], 4),
							",  h1.log.ev.cn= ", round(merge.prob[m.ix, 5], 4), ",  h1.log.ev.snp= ", round(merge.prob[m.ix, 6], 4),
							",  merge.prob.cn= ", round(merge.prob[m.ix, 7], 4), ",  merge.prob.snp= ", round(merge.prob[m.ix, 8], 4), "\n", sep=""))

		}
		annot <- h.seg.dat[["h.snp.annot"]][[m.ix]]
		chr <- c(h.seg.dat[["h.snp.annot"]][[m.ix]][["chr"]], h.seg.dat[["h.cn.annot"]][[m.ix]][["chr"]])
		chr <- ifelse(!all(complete.cases(chr)), max(chr), chr[1])
		snp.pos <- h.seg.dat[["h.snp.annot"]][[m.ix]][["pos"]][length(h.seg.dat[["h.snp.annot"]][[m.ix]][["pos"]])]
		cn.pos <- h.seg.dat[["h.cn.annot"]][[m.ix]][["pos"]][length(h.seg.dat[["h.cn.annot"]][[m.ix]][["pos"]])]
		pos <- max(snp.pos, cn.pos)
		merged.loci <- rbind(merged.loci, c(chr, pos, merge.prob[m.ix, ]) )

		merge.snp.res <- MergeTwoSegsExtreme(m.ix, m.ix + 1, h.seg.dat[["h.snp.d"]], h.seg.dat[["h.snp.annot"]], h.seg.dat[["h.snp.gt.p"]])
		# h.d[["snp"]] <- merge.snp.res[["h.d"]]
		h.seg.dat[["h.snp.d"]] <- merge.snp.res[["h.d"]]
		# h.snp.gt.p <- merge.snp.res[["h.snp.gt.p"]]
		h.seg.dat[["h.snp.gt.p"]] <- merge.snp.res[["h.snp.gt.p"]]
		# h.probe.annot[['snp']] <- merge.snp.res[["h.probe.annot"]]
		h.seg.dat[["h.snp.annot"]] <- merge.snp.res[["h.probe.annot"]]
		
		merge.cn.res <- MergeTwoSegsExtreme(m.ix, m.ix + 1, h.seg.dat[["h.cn.d"]], h.seg.dat[["h.cn.annot"]])
		# h.d[["cn"]] <- merge.cn.res[["h.d"]]
		h.seg.dat[["h.cn.d"]] <- merge.cn.res[["h.d"]]
		# h.probe.annot[['cn']] <- merge.cn.res[["h.probe.annot"]]
		h.seg.dat[["h.cn.annot"]] <- merge.cn.res[["h.probe.annot"]]
		
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
		n.seg <- length(h.seg.dat[[1]])
		seg.chrs <- sapply(h.seg.dat[["h.snp.annot"]], "[[", "chr")
		seg.n.probes <- foreach(dt=d.types, .combine=cbind) %do% sapply(h.seg.dat[[dt]], ncol)
		colnames(seg.n.probes) <- gsub("(^h\\.)|(\\.d$)", "", d.types)
		
		} ## End while(TRUE)	
		return(list(h.seg.dat=h.seg.dat, merge.prob=merge.prob, merged.loci=merged.loci))
	}

	wl.res <- LoadCached({
		wl.res <- WhileFunc()
		saveRDS(wl.res, file.path(results.dir, "jcs.wlres.rds"))	
		wl.res}, 
		cached=T, res.fn=file.path(results.dir, "jcs.wlres.rds"), mod.name="WhileLoopJoinCloseSegs") 
	
	
	h.seg.dat <- wl.res$h.seg.dat
	merged.loci <- wl.res$merged.loci
	merge.prob <- wl.res$merge.prob
	
	if (length(merged.loci) > 0) {
		colnames(merged.loci) <- c("chr", "pos", "prob", "log10_prob", "h0.log.ev.cn", "h0.log.ev.snp", "h1.log.ev.cn", "h1.log.ev.snp" )
		colnames(merge.prob) <- c("merge.prob", "log10.merge.prob", "h0.log.ev.cn", "h0.log.ev.snp", "h1.log.ev.cn", "h1.log.ev.snp", "merge.prob.cn", "merge.prob.snp" )
	}
	
	# list(h.d=list(snp=h.d[["snp"]], cn=h.d[["cn"]]), h.snp.gt.p=h.snp.gt.p, h.probe.annot=list(snp=h.probe.annot[['snp']], cn=h.probe.annot[['cn']]),
	# 		merged.loci=merged.loci, final.merge.prob=merge.prob)
	list(h.seg.dat=h.seg.dat, merged.loci=merged.loci, final.merge.prob=merge.prob)
}


CalcSegMergeProbExtreme <- function(ix.1, ix.2, h.seg.dat, theta, force.diploid, out.p, verbose=FALSE) {
	## H0:  two segs are seperate
	## H1:  two segs are same
	
	if (force.diploid == TRUE) {
		delta.tau <- c(1, 2)
	} else{
		delta.tau <- NA
	}
	h0.res <- CalculateH0EvidenceExtreme(ix.1, ix.2, h.seg.dat, theta, out.p, delta.tau, verbose=verbose)
	h0.log.ev.cn <- h0.res[["cn"]]
	h0.log.ev.snp <- h0.res[["snp"]]
	
	h1.res <- CalculateH1EvidenceExtreme(ix.1, ix.2, h.seg.dat, theta, out.p, delta.tau, verbose=verbose)
	h1.log.ev.cn <- h1.res[["cn"]]
	h1.log.ev.snp <- h1.res[["snp"]]
	
	h0.ev.vec <- c(h0.log.ev.cn, h0.log.ev.snp)
	h0.ev.vec <- h0.ev.vec[complete.cases(h0.ev.vec)]
	h1.ev.vec <- c(h1.log.ev.cn, h1.log.ev.snp)
	h1.ev.vec <- h1.ev.vec[complete.cases(h1.ev.vec)]

	h1.probs <- CalculateH1Probs( sum(h0.ev.vec), sum(h1.ev.vec) )
	h1.probs.cn <- CalculateH1Probs(h0.log.ev.cn , h1.log.ev.cn )
	h1.probs.snp <- CalculateH1Probs(h0.log.ev.snp , h1.log.ev.snp )

	fact.prob <- c(h1.probs.cn[["log10.h1.prob"]], h1.probs.snp[["log10.h1.prob"]])
	fact.prob <- fact.prob[complete.cases(fact.prob)]
	comb.prob <- h1.probs[["log10.h1.prob"]]
	if (abs(sum(fact.prob) - comb.prob) > 1) {
		warning(paste("Factored probs differ by", signif(abs(h1.probs.cn[["log10.h1.prob"]]+h1.probs.snp[["log10.h1.prob"]] - h1.probs[["log10.h1.prob"]])), 2) )
	}
	
	# return(list(merge.prob=h1.probs[["h1.prob"]], log10.merge.prob=h1.probs[["log10.h1.prob"]], merge.prob.cn = h1.probs.cn[["h1.prob"]], merge.prob.snp=h1.probs.snp[["h1.prob"]], h0.log.ev.cn=h0.log.ev.cn, h0.log.ev.snp=h0.log.ev.snp, h1.log.ev.cn=h1.log.ev.cn, h1.log.ev.snp=h1.log.ev.snp))

	return(data.frame(merge.prob=h1.probs[["h1.prob"]], log10.merge.prob=h1.probs[["log10.h1.prob"]], h0.log.ev.cn=h0.log.ev.cn, h0.log.ev.snp=h0.log.ev.snp, h1.log.ev.cn=h1.log.ev.cn, h1.log.ev.snp=h1.log.ev.snp, merge.prob.cn = h1.probs.cn[["h1.prob"]], merge.prob.snp=h1.probs.snp[["h1.prob"]]))

}

CalcSegMergeProb <- function(ix.1, ix.2, h.d, h.snp.gt.p, theta,
                             force.diploid, out.p, verbose=FALSE) {
  ## H0:  two segs are seperate
  ## H1:  two segs are same
  if (force.diploid == TRUE) {
    e.mu <- c(1, 1, 2)
  } else{
    e.mu <- NA
  }
  h0.log.ev <- CalculateH0Evidence(ix.1, ix.2, h.d, h.snp.gt.p,
                                   theta, out.p, e.mu, verbose=verbose)
  h1.log.ev <- CalculateH1Evidence(ix.1, ix.2, h.d, h.snp.gt.p,
                                   theta, out.p, e.mu, verbose=verbose)
  h1.probs <- CalculateH1Probs(h0.log.ev, h1.log.ev)

  return(list(merge.prob=h1.probs[["h1.prob"]],
              log10.merge.prob=h1.probs[["log10.h1.prob"]]))
}

CalculateH1EvidenceExtreme <- function(ix.1, ix.2, h.seg.dat, theta, out.p, delta.tau, verbose=FALSE) {
	   
	# browser()
	mrg.snp.d <- cbind(h.seg.dat[["h.snp.d"]][[ix.1]], h.seg.dat[["h.snp.d"]][[ix.2]])
 	mrg.cn.d <- cbind(h.seg.dat[["h.cn.d"]][[ix.1]], h.seg.dat[["h.cn.d"]][[ix.2]])
	mrg.het.prob <- rbind(h.seg.dat[["h.snp.gt.p"]][[ix.1]], h.seg.dat[["h.snp.gt.p"]][[ix.2]])
	
	mrg.h.seg.dat <- h.seg.dat
	n.segs <- length(h.seg.dat[[1]] )
	# Merge the data for segs ix.1 and ix.2 into ix.1 to get H1 Evidence.  
	# Note that this isn't comprehensive merge and that there will be other lists (eg, h.snp.annot) that will be off, but these are not returned from this function
	if (ix.1 == 1) {
		mrg.h.seg.dat[["h.snp.d"]][[1]] <- mrg.snp.d
		mrg.h.seg.dat[["h.cn.d"]][[1]] <- mrg.cn.d
		mrg.h.seg.dat[["h.snp.gt.p"]][[1]] <- mrg.het.prob

		mrg.h.seg.dat[["h.snp.d"]][2:(n.segs-1)] <- h.seg.dat[["h.snp.d"]][3:n.segs]
		mrg.h.seg.dat[["h.cn.d"]][2:(n.segs-1)] <- h.seg.dat[["h.cn.d"]][3:n.segs]
		mrg.h.seg.dat[["h.snp.gt.p"]][2:(n.segs-1)] <- h.seg.dat[["h.snp.gt.p"]][3:n.segs]

		mrg.h.seg.dat[["h.snp.d"]][[n.segs]] <- NULL
		mrg.h.seg.dat[["h.cn.d"]][[n.segs]] <- NULL
		mrg.h.seg.dat[["h.snp.gt.p"]][[n.segs]] <- NULL

	} else if (ix.1 == n.segs - 1) {
		mrg.h.seg.dat[["h.snp.d"]][[n.segs-1]] <- mrg.snp.d
		mrg.h.seg.dat[["h.cn.d"]][[n.segs-1]] <- mrg.cn.d
		mrg.h.seg.dat[["h.snp.gt.p"]][[n.segs-1]] <- mrg.het.prob

		mrg.h.seg.dat[["h.snp.d"]][[n.segs]] <- NULL
		mrg.h.seg.dat[["h.cn.d"]][[n.segs]] <- NULL
		mrg.h.seg.dat[["h.snp.gt.p"]][[n.segs]] <- NULL
	} else {
		mrg.h.seg.dat[["h.snp.d"]][[ix.1]] <- mrg.snp.d
		mrg.h.seg.dat[["h.cn.d"]][[ix.1]] <- mrg.cn.d
		mrg.h.seg.dat[["h.snp.gt.p"]][[ix.1]] <- mrg.het.prob

		mrg.h.seg.dat[["h.snp.d"]][ix.2:(n.segs-1)] <- h.seg.dat[["h.snp.d"]][(ix.2+1):(n.segs)]
		mrg.h.seg.dat[["h.cn.d"]][ix.2:(n.segs-1)] <- h.seg.dat[["h.cn.d"]][(ix.2+1):(n.segs)]
		mrg.h.seg.dat[["h.snp.gt.p"]][ix.2:(n.segs-1)] <- h.seg.dat[["h.snp.gt.p"]][(ix.2+1):(n.segs)]

		mrg.h.seg.dat[["h.snp.d"]][[n.segs]] <- NULL
		mrg.h.seg.dat[["h.cn.d"]][[n.segs]] <- NULL
		mrg.h.seg.dat[["h.snp.gt.p"]][[n.segs]] <- NULL
	}
	
	res <- GetSegLogEvExtreme(ix.1, mrg.h.seg.dat, theta, out.p, delta.tau, verbose=verbose)
	log.ev.cn <- res[["log.ev.cn"]]
	log.ev.snp <- res[["log.ev.snp"]]
	return(list(cn=log.ev.cn, snp=log.ev.snp))
}

CalculateH1Evidence <- function(ix.1, ix.2, h.d, h.snp.gt.p, theta, out.p, e.mu, verbose=FALSE) {
  
  mrg.d <- cbind(h.d[[ix.1]], h.d[[ix.2]])
  mrg.het.prob <- rbind(h.snp.gt.p[[ix.1]], h.snp.gt.p[[ix.2]])

  return(GetSegLogEv(mrg.d, mrg.het.prob, theta, out.p, e.mu, verbose=verbose)[1])
}


CalculateH0EvidenceExtreme <- function(ix.1, ix.2, h.seg.dat, theta, out.p, delta.tau, verbose=FALSE) {
	# GetLogEvN <- function(ix, verbose=FALSE) {
	# 	# h.d.i = list(snp=h.d[["snp"]][[ix]], cn=h.d[["cn"]][[ix]])
	# 	GetSegLogEvExtreme(ix, h.seg.dat, theta, out.p, delta.tau, verbose=verbose)
	# }
	
	# ev1 <- GetLogEvN(ix.1, verbose=verbose)
	# ev2 <- GetLogEvN(ix.2, verbose=verbose)
	
	ev1 <- GetSegLogEvExtreme(ix.1, h.seg.dat, theta, out.p, delta.tau, verbose=verbose)
	ev2 <- GetSegLogEvExtreme(ix.2, h.seg.dat, theta, out.p, delta.tau, verbose=verbose)
	
	log.ev1 <- ev1[["log.ev"]]
	log.ev2 <- ev2[["log.ev"]]
	log.ev1.cn <- ev1[["log.ev.cn"]]
	log.ev2.cn <- ev2[["log.ev.cn"]]
	log.ev1.snp <- ev1[["log.ev.snp"]]
	log.ev2.snp <- ev2[["log.ev.snp"]]
	
	# return(log.ev1 + log.ev2)
	return(list(cn=log.ev1.cn + log.ev2.cn, snp = log.ev1.snp + log.ev2.snp))
}


CalculateH0Evidence <- function(ix.1, ix.2, h.d, h.snp.gt.p,
                                theta, out.p, e.mu, verbose=FALSE) {
  GetLogEvN <- function(ix, verbose=FALSE) {
    GetSegLogEv(h.d[[ix]], h.snp.gt.p[[ix]], theta, out.p, e.mu,
                verbose=verbose)
  }

  log.ev1 <- GetLogEvN(ix.1, verbose=verbose)[1]
  log.ev2 <- GetLogEvN(ix.2, verbose=verbose)[1]

  return(log.ev1 + log.ev2)
}

CalculateH1Probs <- function(h0.log.ev, h1.log.ev) {
  ve <- c(h0.log.ev, h1.log.ev)
  mx <- max(ve)
  nf <- sum(exp(ve - mx))
  dens <- exp(ve - mx) / nf
  h1.prob <- dens[2]

  log.dens <- (ve - mx) - log(nf)
  log10.dens <- log.dens * log(exp(1), 10)
  log10.h1.prob <- log10.dens[2]

  return(list(h1.prob=h1.prob, log10.h1.prob=log10.h1.prob))
}




GetSegLogEvExtreme <- function(idx, h.seg.dat, theta, out.p, delta.tau=NA, verbose=FALSE) {
	
	het.prob <- h.seg.dat[["h.snp.gt.p"]][[idx]]
	GetSnpSegLL <- function(x, d, het.prob, theta, out.p) {
		dist <- x[1]
		total <- x[2]
		return(sum(CalcSnpLogLik(d, delta.tau=x, out.p, het.prob, theta)))
	}
	
	GetCnSegLL = function(x, d, theta, out.p ) {
		dist <- x[1]
		total <- x[2]
		return(sum(CalcCnLogLik(d, delta.tau=x, out.p, theta) ))
	}
	

	
	if (is.na(delta.tau)) {
		init.tau <- sapply(seq_along(h.seg.dat[["h.snp.d"]]), function(i) {
				snp.d <- colSums(h.seg.dat[["h.snp.d"]][[i]])
				cn.d <- h.seg.dat[["h.cn.d"]][[i]]
				capseg.d <- h.seg.dat[["h.capseg.d"]][[i]]
				median(c(snp.d, cn.d, capseg.d ) ) 
			}
		)
		quart <- init.tau / 4   
		init.e.mu <- c(quart, init.tau - quart, init.tau)
		init.delta <- 2 * (init.tau - quart) - init.tau
		
		delta.tau <- GridstartH1OptMeansExtreme(idx, h.seg.dat, delta.tau=cbind(init.delta, init.tau), out.p, theta, verbose=verbose)
	}
	
	tau <- delta.tau[2]
	dist <- delta.tau[1]
	
	hess.mat.snp <- hessian(GetSnpSegLL, c(dist, tau), "Richardson", d=h.seg.dat[["h.snp.d"]][[idx]], het.prob=het.prob, theta=theta, out.p=out.p) 
	curv.snp <- abs(det(hess.mat.snp / (2 * pi)))
	ll.snp <- GetSnpSegLL(c(dist, tau), h.seg.dat[["h.snp.d"]][[idx]], het.prob, theta, out.p); 
	log.ev.snp <- ll.snp + (log((curv.snp)^(-1 / 2))) - (log((5 / 2)^2))
	
	hess.mat.cn <- hessian(GetCnSegLL, c(dist, tau), "Richardson", d=h.seg.dat[["h.cn.d"]][[idx]], theta=theta, out.p=out.p) 
	curv.cn <- abs(det(hess.mat.cn / (2 * pi)))
	ll.cn <- GetCnSegLL(c(dist, tau), h.seg.dat[["h.cn.d"]][[idx]], theta, out.p); 
	log.ev.cn <- ll.cn + (log((curv.cn)^(-1 / 2))) - log(5)
	
	if (all(c("snp", "cn") %in% PROBE.TYPES)) {
		return(list(log.ev=log.ev.snp + log.ev.cn, ll=ll.snp + ll.cn, log.ev.snp=log.ev.snp, log.ev.cn=log.ev.cn, ll.snp=ll.snp, ll.cn=ll.cn))
	} else if ("snp" %in% PROBE.TYPES) {
		return(list(log.ev=log.ev.snp, ll=ll.snp, log.ev.snp=log.ev.snp, log.ev.cn=log.ev.cn, ll.snp=ll.snp, ll.cn=ll.cn))
	} else if ("cn" %in% PROBE.TYPES) {
		return(list(log.ev=log.ev.cn, ll=ll.cn, log.ev.snp=log.ev.snp, log.ev.cn=log.ev.cn, ll.snp=ll.snp, ll.cn=ll.cn))
	} else {
		stop ("cn or snp wasn't found in PROBE.TYPES when trying to GetSegLogEvExtreme")
	}
}

GetSegLogEv <- function(d, het.prob, theta, out.p, e.mu=NA, verbose=FALSE) {
  GetSegLL <- function(x, d, het.prob, theta, out.p) {
    dist <- x[1]
    total <- x[2]
    e.mu <-  c((total / 2 - dist / 2), (total / 2 + dist / 2), total)

    return(sum(CalcSnpLogLik(d, e.mu, out.p, het.prob, theta)))
  }

  if (is.na(e.mu)) {
    init.tau <- median(colSums(d)) 
    quart <- init.tau / 4   
    init.e.mu <- c(quart, init.tau - quart, init.tau)
    e.mu <- GridstartH1OptMeans(d, init.e.mu, out.p, het.prob, theta,
                                verbose=verbose) 
  }
  
  tau <- e.mu[3]
  dist <- abs(e.mu[2] - e.mu[1])

  hess.mat <- hessian(GetSegLL, c(dist, tau), "Richardson", d=d,
                      het.prob=het.prob, theta=theta, out.p=out.p) 
  curv <- abs(det(hess.mat / (2 * pi)))
  ll <- GetSegLL(c(dist, tau), d, het.prob, theta, out.p)
  log.ev <- ll + (log((curv)^(-1 / 2))) - (log((5 / 2)^2))

  return(c(log.ev, ll, curv))
}
