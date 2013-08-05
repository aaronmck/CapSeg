## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


MergeSegs = function(seg1, seg2) {    
  # seg1 = array.seg.dat[[1]]
  # seg2 = capseg.seg.dat[[1]]
  	chr = unique(c(seg1[,"Chromosome"], seg2[,"Chromosome"]))
  
  	out = foreach(x=chr, .combine=rbind) %dopar% {
  		# x= 1
  		bp = sort(c(seg1[seg1[, "Chromosome"] == x, "Start.bp"], seg2[seg2[, "Chromosome"] == x, "Start.bp"], seg1[seg1[, "Chromosome"] == x, "End.bp"], seg2[seg2[, "Chromosome"] == x, "End.bp"]))
  		starts = bp[1:(length(bp)-1)]
  		ends = bp[2:length(bp)]-1
  		ends[length(ends)] <- ends[length(ends)] + 1
  		return(data.frame(Chromosome = rep(x, length(bp)-1), Start.bp=starts, End.bp=ends, stringsAsFactors=F))
   	}

}




ConstructAsSegMat <- function(h.d, h.snp.annot, em.fit) {
	n.seg <- length(h.snp.annot)
	seg.chrs <- rep(NA, n.seg)
	for (i in 1:n.seg) {
		seg.chrs[i] <- h.snp.annot[[i]][["chr"]]
	}
	seg.n.snps <- sapply(h.d, ncol)
	
	cols <- c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", 
			"A1.Seg.CN", "A2.Seg.CN", "tCN.Seg.sd", "AS.Seg.sd" )
	as.seg.mat <- matrix(nrow=0, ncol=length(cols))
	colnames(as.seg.mat) <- cols
	
	for (i in 1:n.seg) {
		start.bp <- h.snp.annot[[i]][["pos"]][1]
		end.bp <- h.snp.annot[[i]][["pos"]][seg.n.snps[i]]
		
		as.seg.mat <- rbind(as.seg.mat, c(seg.chrs[i], 
						start.bp,
						end.bp, 
						seg.n.snps[i],
						end.bp - start.bp,
						em.fit[["e.mu"]][i, 1], 
						em.fit[["e.mu"]][i, 2], 
						em.fit[["mu.post.sd"]][i, "tau"], 
						em.fit[["mu.post.sd"]][i, "delta"]))
	}
	
	cols <- colnames(as.seg.mat)
	as.seg.mat <- as.data.frame(as.seg.mat)
	
	as.seg.mat[as.seg.mat[,1] == "X", 1] <- 23
	as.seg.mat[as.seg.mat[,1] == "Y", 1] <- 24
	
	as.seg.mat <- as.matrix(as.seg.mat)
	as.seg.mat <- matrix(as.numeric(as.seg.mat), ncol=length(cols))
	colnames(as.seg.mat) = cols
	
	return(as.seg.mat)
}

PostCalibrateAsDat <- function(dat, snp.tx, adj.atten=FALSE, verbose=FALSE) {
	InvAtten <- function(r, at) {
		r / (1 + at - (at * r))
	}
	
	## align snp.tx to dat
	dat <- dat[intersect(rownames(dat), rownames(snp.tx)), ]
	snp.tx <- snp.tx[rownames(dat), ]
	
	n.snp <- nrow(snp.tx)
	tx.snp.dat <- array(NA, dim=c(n.snp, 2))
	
	## affine tx
	tx.snp.dat[, 1] <- dat[, 1] + snp.tx[, 2] * dat[, 2]
	tx.snp.dat[, 2] <- dat[, 1] * snp.tx[, 1] + dat[, 2] *
			(1 + apply(snp.tx[, c(1:2)], 1, prod))
	
	## re-center and scale
	tx.snp.dat <- (tx.snp.dat + snp.tx[, c(4, 3)]) * snp.tx[, c(6, 5)] 
	
	snp.tx[is.na(snp.tx)] <- 0
	
	## Attenuation adjustment 
	if (adj.atten) {
		a.at.hat <- snp.tx[, 7]
		b.at.hat <- snp.tx[, 8]
		
		tx.snp.dat[, 1] <- InvAtten(tx.snp.dat[,1], a.at.hat)
		tx.snp.dat[, 2] <- InvAtten(tx.snp.dat[,2], b.at.hat)
	}
	
	return(tx.snp.dat)
}

GetAlleleSegData <- function(snp.d, cn.d, capseg.d, germline.hets, snp.annot, glad.mat, snp.freqs,
		use.pop, use.normal, normal,
		bad.snps, dbSNP.annot, impute.gt,
		mn.calls.fn, mn.sample, verbose=FALSE) {
	
	## First process SNP probes
	
	## remove SNPs with NAs
	sna <- apply(is.na(snp.d), 1, sum)
	snp.d <- snp.d[sna == 0, ]
	
	snp.d <- snp.d[intersect(rownames(snp.annot), rownames(snp.d)), ]
	snp.d <- cbind(snp.annot[rownames(snp.d), ], snp.d)
	
	## remove 'Bad' SNPs
	if (!is.na(bad.snps)) {
		bad <- which(rownames(snp.d) %in% bad.snps) 
		snp.d <- snp.d[-bad, ]
		if (verbose) {
			print(paste("Removing ", sum(bad), " \'bad\' SNPs", sep=""))
		}
	}
	
	## check for matched-normal sample
	found.matched.normal <- FALSE
	mn.calls <- NULL
	
	## read matched-norm GT calls
	## These days mn.sample should only be passed in as a string or NULL,
	## but historically it was NA or even "NA" so catch these cases just
	## in case old scripts are used
	if (!normal && use.normal && (!is.null(mn.sample)) &&
			(!is.na(mn.sample)) && (mn.sample != "NA")) {
		if (is.null(mn.calls.fn)) {
			stop("Matched normal sample specified, but no calls file provided")
		}
		if (verbose) {
			print("reading column data")
		}
		mn.calls <- ReadCol(mn.calls.fn, mn.sample, save.rownames=TRUE)
	}
	
	if (verbose) {
		print(summary(mn.calls))
	}
	
	if (use.normal && (!is.null(mn.calls))) {
		found.matched.normal <- TRUE
		
		mn.calls <- mn.calls[rownames(snp.d), ]
		snp.gt.p <- matrix(0, nrow=4, ncol=length(mn.calls), byrow=TRUE)
		snp.gt.p[, mn.calls == 1] <- c(0.015, 0.985, 0.015, 0.985) / 2
		snp.gt.p[, mn.calls == 0] <- c(0.005, 0.005, 0.985, 0.005) 
		snp.gt.p[, mn.calls == 2] <- c(0.985, 0.005, 0.005, 0.005) 
		snp.gt.p[, mn.calls == -1] <- c(0.39, 0.11, 0.39, 0.11)
		
		snp.gt.p <- t(snp.gt.p)
		
		if (verbose) {
			n.het <- sum(mn.calls == 1)
			print(paste("Normal ", round(n.het / length(mn.calls) * 100, 2), "% Het",
							sep=""))
		}
	} else {
		if (is.na(use.pop) || use.pop == "NA") {
			if (verbose) {
				print("Using Null pop.")
			}
			snp.gt.p <- matrix(c(0.39, 0.11, 0.39, 0.11), nrow=nrow(snp.d), ncol=4,
					byrow=TRUE)
		} else {
			if (verbose) {
				print(paste( "Using population allele-frequency data for: ",
								use.pop, sep=""))
			}
			## lookup population allele-frequencies
			pop.cols <- c(paste(use.pop, "_A",sep=""),  paste(use.pop, "_B", sep=""))
			
			## match SNP names
			snp.freqs <- snp.freqs[intersect(rownames(snp.d), rownames(snp.freqs)), ]
			
			## pseduo-count
			freqs <- snp.freqs[, pop.cols] + 0.01
			freqs <- freqs / rowSums(freqs)
			
			aa.probs <- freqs[, 1]^2
			bb.probs <- freqs[, 2]^2
			het.probs <- 1 - (aa.probs + bb.probs)
			
			snp.gt.p <- rbind(bb.probs, het.probs / 2, aa.probs, het.probs / 2) 
			colnames(snp.gt.p) <- rownames(snp.d)
			rownames(snp.gt.p)[c(1, 3)] <- c("BB_probs", "AA_probs")
			
			na.ix <- apply(is.na(snp.gt.p), 2, sum) > 0
			
			snp.gt.p[, na.ix] <- apply(snp.gt.p, 1, mean, na.rm=TRUE)
			snp.gt.p <- t(snp.gt.p)  
		}
	}
	
	snp.d[, 1] <- as.character(snp.d[, 1])
	
	
	## Then process CN probes
	
	## remove CNs with NAs
	
	sna <- apply(is.na(cn.d), 1, sum)
	cn.d <- cn.d[sna == 0, , drop=F]
	
	cn.d <- cn.d[intersect(rownames(snp.annot), rownames(cn.d)), , drop=F]
	cn.d <- cbind(snp.annot[rownames(cn.d), ], cn.d)
	cn.d[["Chromosome"]] <- as.character(cn.d[["Chromosome"]])
	
	## Using all probe evidence, split into segments.
	## Now, small segments will only be dropped if the sum of all probes types in 
	## a given segment is less than the minimum allowed number of probes.
	
	h.seg.dat <- GetAsSegs(glad.mat, snp.d, cn.d, capseg.d, germline.hets, snp.gt.p, dbSNP.annot, impute.gt, verbose=verbose)
	
	return(list(h.seg.dat=h.seg.dat, mn.sample=mn.sample, found.matched.normal=found.matched.normal))
}

GetArrayAlleleSegData <- function(snp.d, cn.d, snp.annot, glad.mat, snp.freqs,
		use.pop, use.normal, normal,
		bad.snps, dbSNP.annot, impute.gt,
		mn.calls.fn, mn.sample, verbose=FALSE) {
	
	## First process SNP probes
	
	## remove SNPs with NAs
	sna <- apply(is.na(snp.d), 1, sum)
	snp.d <- snp.d[sna == 0, ]
	
	snp.d <- snp.d[intersect(rownames(snp.annot), rownames(snp.d)), ]
	snp.d <- cbind(snp.annot[rownames(snp.d), ], snp.d)
	
	## remove 'Bad' SNPs
	if (!is.na(bad.snps)) {
		bad <- which(rownames(snp.d) %in% bad.snps) 
		snp.d <- snp.d[-bad, ]
		if (verbose) {
			print(paste("Removing ", sum(bad), " \'bad\' SNPs", sep=""))
		}
	}
	
	## check for matched-normal sample
	found.matched.normal <- FALSE
	mn.calls <- NULL
	
	## read matched-norm GT calls
	## These days mn.sample should only be passed in as a string or NULL,
	## but historically it was NA or even "NA" so catch these cases just
	## in case old scripts are used
	if (!normal && use.normal && (!is.null(mn.sample)) &&
			(!is.na(mn.sample)) && (mn.sample != "NA")) {
		if (is.null(mn.calls.fn)) {
			stop("Matched normal sample specified, but no calls file provided")
		}
		if (verbose) {
			print("reading column data")
		}
		mn.calls <- ReadCol(mn.calls.fn, mn.sample, save.rownames=TRUE)
	}
	
	if (verbose) {
		print(summary(mn.calls))
	}
	
	if (use.normal && (!is.null(mn.calls))) {
		found.matched.normal <- TRUE
		
		mn.calls <- mn.calls[rownames(snp.d), ]
		snp.gt.p <- matrix(0, nrow=4, ncol=length(mn.calls), byrow=TRUE)
		snp.gt.p[, mn.calls == 1] <- c(0.015, 0.985, 0.015, 0.985) / 2
		snp.gt.p[, mn.calls == 0] <- c(0.005, 0.005, 0.985, 0.005) 
		snp.gt.p[, mn.calls == 2] <- c(0.985, 0.005, 0.005, 0.005) 
		snp.gt.p[, mn.calls == -1] <- c(0.39, 0.11, 0.39, 0.11)
		
		snp.gt.p <- t(snp.gt.p)
		
		if (verbose) {
			n.het <- sum(mn.calls == 1)
			print(paste("Normal ", round(n.het / length(mn.calls) * 100, 2), "% Het",
							sep=""))
		}
	} else {
		if (is.na(use.pop) || use.pop == "NA") {
			if (verbose) {
				print("Using Null pop.")
			}
			snp.gt.p <- matrix(c(0.39, 0.11, 0.39, 0.11), nrow=nrow(snp.d), ncol=4,
					byrow=TRUE)
		} else {
			if (verbose) {
				print(paste( "Using population allele-frequency data for: ",
								use.pop, sep=""))
			}
			## lookup population allele-frequencies
			pop.cols <- c(paste(use.pop, "_A",sep=""),  paste(use.pop, "_B", sep=""))
			
			## match SNP names
			snp.freqs <- snp.freqs[intersect(rownames(snp.d), rownames(snp.freqs)), ]
			
			## pseduo-count
			freqs <- snp.freqs[, pop.cols] + 0.01
			freqs <- freqs / rowSums(freqs)
			
			aa.probs <- freqs[, 1]^2
			bb.probs <- freqs[, 2]^2
			het.probs <- 1 - (aa.probs + bb.probs)
			
			snp.gt.p <- rbind(bb.probs, het.probs / 2, aa.probs, het.probs / 2) 
			colnames(snp.gt.p) <- rownames(snp.d)
			rownames(snp.gt.p)[c(1, 3)] <- c("BB_probs", "AA_probs")
			
			na.ix <- apply(is.na(snp.gt.p), 2, sum) > 0
			
			snp.gt.p[, na.ix] <- apply(snp.gt.p, 1, mean, na.rm=TRUE)
			snp.gt.p <- t(snp.gt.p)  
		}
	}
	
	snp.d[, 1] <- as.character(snp.d[, 1])
	
	
	## Then process CN probes
	
	## remove CNs with NAs
	
	sna <- apply(is.na(cn.d), 1, sum)
	cn.d <- cn.d[sna == 0, , drop=F]
	
	cn.d <- cn.d[intersect(rownames(snp.annot), rownames(cn.d)), , drop=F]
	cn.d <- cbind(snp.annot[rownames(cn.d), ], cn.d)
	cn.d[["Chromosome"]] <- as.character(cn.d[["Chromosome"]])
	
	## Using all probe evidence, split into segments.
	## Now, small segments will only be dropped if the sum of all probes types in 
	## a given segment is less than the minimum allowed number of probes.
	
	h.seg.dat <- GetArrayAsSegs(glad.mat, snp.d, cn.d, snp.gt.p, dbSNP.annot, impute.gt, verbose=verbose)
	
	return(list(h.seg.dat=h.seg.dat, mn.sample=mn.sample, found.matched.normal=found.matched.normal))
}


GetAsSegs <- function(glad.mat, snp.d, cn.d, capseg.d, germline.hets, snp.gt.p, dbSNP.annot, impute.gt, verbose=FALSE) {
	GetProbeIx <- function(i, probe.dat) {
		(probe.dat[, 1] == glad.mat[i, "Chromosome"]) &
					((probe.dat[, 2] >= glad.mat[i, "Start.bp"]) &
					(probe.dat[, 2] <= glad.mat[i, "End.bp"]))
	}
	
	if (verbose) {
		print(paste("GetAsSegs .... processing", nrow(glad.mat), "segments"))
	}
	##  check number of total probes in each segment
	n.seg <- nrow(glad.mat)
	n.snp.probes <- foreach(i = 1:n.seg, .combine=c) %dopar% {
		sum(GetProbeIx(i, snp.d)) 
	}
	n.cn.probes <- foreach(i = 1:n.seg, .combine=c) %dopar% {
		sum(GetProbeIx(i, cn.d))
	}
	n.capseg.probes <- foreach(i = 1:n.seg, .combine=c) %dopar% {
		sum(GetProbeIx(i, capseg.d))
	}
		
	# too.small.ix = n.snp.probes + n.cn.probes + n.capseg.probes > 2
	# glad.mat <- glad.mat[, ]
	
	if (verbose) {
		print(paste(nrow(glad.mat), " segs with > 2 SNPs and > 2 CN probes", sep=""))
	}
	n.seg <- nrow(glad.mat)
	
	if (verbose) {
		print(paste("get_AS_segs_ n.seg = ",n.seg))
	}
	
	
	snp.res = foreach(i = 1:n.seg) %dopar% {
		snp.ix <- GetProbeIx(i, snp.d)
		raw.snp.seg <- snp.d[snp.ix, c(3, 4), drop=FALSE]
		if (verbose) {
			print(paste("Processing Segment: ", i, "SNP size = ", dim(raw.snp.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"]))
		}
		
		h.snp.gt.p <- snp.gt.p[snp.ix, , drop=FALSE]
		rownames(h.snp.gt.p) <- rownames(snp.d)[snp.ix]
		
		h.snp.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=snp.d[snp.ix, 2]) 
		
		if (impute.gt) {
			dbSNP.ix <- match(rownames(snp.d)[snp.ix], rownames(dbSNP.annot))
			h.snp.annot[["dbSNP"]] <- dbSNP.annot[dbSNP.ix, , drop=FALSE] 
		}
		
		list(h.snp.d=t(raw.snp.seg), h.snp.gt.p=h.snp.gt.p, h.snp.annot=h.snp.annot)
	}
	h.snp.d = lapply(snp.res, "[[", "h.snp.d")
	h.snp.gt.p = lapply(snp.res, "[[", "h.snp.gt.p")
	h.snp.annot = lapply(snp.res, "[[", "h.snp.annot")

	
	# CN Probes
	# for (i in 1:n.seg) {
	cn.res = foreach(i = 1:n.seg) %dopar% {
		cn.ix <- GetProbeIx(i, cn.d)
		raw.cn.seg <- cn.d[cn.ix, "Intensity", drop=FALSE]
		if (verbose) {
			print(paste("Processing Segment: ", i, "CN size = ", dim(raw.cn.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"]))
		}
		h.cn.d <- t(raw.cn.seg)
		h.cn.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=cn.d[cn.ix, 2])
		list(h.cn.d=h.cn.d, h.cn.annot=h.cn.annot)
	}
	h.cn.d = lapply(cn.res, "[[", "h.cn.d")
	h.cn.annot = lapply(cn.res, "[[", "h.cn.annot")
	
	# Capseg Probes
	# for (i in 1:n.seg) {
	capseg.res = foreach(i = 1:n.seg) %dopar% {
		capseg.ix <- GetProbeIx(i, capseg.d)
		raw.capseg.seg <- capseg.d[capseg.ix, "Intensity", drop=FALSE]
		if (verbose) {
			print(paste("Processing Segment: ", i, "Capseg size = ", dim(raw.capseg.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"]))
		}
		h.capseg.d <- t(raw.capseg.seg)
		h.capseg.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=capseg.d[capseg.ix, 2])
		list(h.capseg.d=h.capseg.d, h.capseg.annot=h.capseg.annot)
	}
	h.capseg.d = lapply(capseg.res, "[[", "h.capseg.d")
	h.capseg.annot = lapply(capseg.res, "[[", "h.capseg.annot")
	
	# Germline Het Allele Counts
	GetHetIx <- function(i, het.dat) {
		(het.dat$Chromosome == glad.mat[i, "Chromosome"]) & (het.dat$Start_position >= glad.mat[i, "Start.bp"]) & (het.dat$Start_position <= glad.mat[i, "End.bp"])
	}
	gh.res = foreach(i = 1:n.seg) %dopar% {
		gh.ix <- GetHetIx(i, germline.hets)
		raw.gh.seg <- germline.hets[gh.ix, c("i_t_ref_count", "i_t_alt_count"), drop=FALSE]
		rownames(raw.gh.seg) = apply(germline.hets[gh.ix, c("Chromosome", "Start_position", "Hugo_Symbol"), drop=FALSE], 1, function(x) paste(trim(x), collapse="_"))
		
		if (verbose) {
			print(paste("Processing Segment: ", i, "Germline Het size = ", dim(raw.gh.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"], ifelse(dim(raw.gh.seg)[1]==0, ": No HETS", "")))
		}
		gh.wes.allele.d <- t(raw.gh.seg)
		gh.wes.allele.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=germline.hets$Start_position[gh.ix], hugo.symbol=germline.hets$Hugo_Symbol[gh.ix], 
				ref.allele=germline.hets$Reference_Allele[gh.ix], alt.allele=germline.hets$Tumor_Seq_Allele1[gh.ix], dbSNP=germline.hets$dbSNP[gh.ix])
		list(gh.wes.allele.d=gh.wes.allele.d, gh.wes.allele.annot=gh.wes.allele.annot)
	}
	gh.wes.allele.d = lapply(gh.res, "[[", "gh.wes.allele.d")
	gh.wes.allele.annot = lapply(gh.res, "[[", "gh.wes.allele.annot")
	
	return(list(h.snp.d=h.snp.d, h.cn.d=h.cn.d, h.capseg.d=h.capseg.d, h.snp.gt.p=h.snp.gt.p, h.snp.annot=h.snp.annot,
				h.cn.annot=h.cn.annot, h.capseg.annot=h.capseg.annot,  gh.wes.allele.d=gh.wes.allele.d, gh.wes.allele.annot=gh.wes.allele.annot))
}

GetCaptureAsSegs <- function(glad.mat, capseg.d, germline.hets, verbose=FALSE) {
	GetProbeIx <- function(i, probe.dat) {
		(probe.dat[, 1] == glad.mat[i, "Chromosome"]) &
					((probe.dat[, 2] >= glad.mat[i, "Start.bp"]) &
					(probe.dat[, 2] <= glad.mat[i, "End.bp"]))
	}
	
	if (verbose) {
		print(paste("GetCaptureAsSegs .... processing", nrow(glad.mat), "segments"))
	}
	##  check number of total probes in each segment
	n.seg <- nrow(glad.mat)
	n.capseg.probes <- foreach(i = 1:n.seg, .combine=c) %dopar% {
		sum(GetProbeIx(i, capseg.d))
	}
		
	# too.small.ix = n.snp.probes + n.cn.probes + n.capseg.probes > 2
	# glad.mat <- glad.mat[, ]
	
	n.seg <- nrow(glad.mat)
	
	if (verbose) {
		print(paste("GetCaptureAsSegs n.seg = ",n.seg))
	}
	
	capseg.res = foreach(i = 1:n.seg) %dopar% {
		capseg.ix <- GetProbeIx(i, capseg.d)
		raw.capseg.seg <- capseg.d[capseg.ix, "Intensity", drop=FALSE]
		if (verbose) {
			print(paste("Processing Segment: ", i, "Capseg size = ", dim(raw.capseg.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"]))
		}
		h.capseg.d <- t(raw.capseg.seg)
		h.capseg.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=capseg.d[capseg.ix, 2])
		list(h.capseg.d=h.capseg.d, h.capseg.annot=h.capseg.annot)
	}
	h.capseg.d = lapply(capseg.res, "[[", "h.capseg.d")
	h.capseg.annot = lapply(capseg.res, "[[", "h.capseg.annot")
	
	# Germline Het Allele Counts
	GetHetIx <- function(i, het.dat) {
		(het.dat$Chromosome == glad.mat[i, "Chromosome"]) & (het.dat$Start_position >= glad.mat[i, "Start.bp"]) & (het.dat$Start_position <= glad.mat[i, "End.bp"])
	}
	gh.res = foreach(i = 1:n.seg) %dopar% {
		gh.ix <- GetHetIx(i, germline.hets)
		raw.gh.seg <- germline.hets[gh.ix, c("i_t_ref_count", "i_t_alt_count"), drop=FALSE]
		rownames(raw.gh.seg) = apply(germline.hets[gh.ix, c("Chromosome", "Start_position", "Hugo_Symbol"), drop=FALSE], 1, function(x) paste(trim(x), collapse="_"))
		
		if (verbose) {
			print(paste("Processing Segment: ", i, "Germline Het size = ", dim(raw.gh.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"], ifelse(dim(raw.gh.seg)[1]==0, ": No HETS", "")))
		}
		gh.wes.allele.d <- t(raw.gh.seg)
		gh.wes.allele.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=germline.hets$Start_position[gh.ix], hugo.symbol=germline.hets$Hugo_Symbol[gh.ix], 
				ref.allele=germline.hets$Reference_Allele[gh.ix], alt.allele=germline.hets$Tumor_Seq_Allele1[gh.ix], dbSNP=germline.hets$dbSNP[gh.ix])
		list(gh.wes.allele.d=gh.wes.allele.d, gh.wes.allele.annot=gh.wes.allele.annot)
	}
	gh.wes.allele.d = lapply(gh.res, "[[", "gh.wes.allele.d")
	gh.wes.allele.annot = lapply(gh.res, "[[", "gh.wes.allele.annot")
	
	return(list(h.capseg.d=h.capseg.d, h.capseg.annot=h.capseg.annot,  gh.wes.allele.d=gh.wes.allele.d, gh.wes.allele.annot=gh.wes.allele.annot))
}

GetArrayAsSegs <- function(glad.mat, snp.d, cn.d, snp.gt.p, dbSNP.annot, impute.gt, verbose=FALSE) {
	GetProbeIx <- function(i, probe.dat) {
		(probe.dat[, 1] == glad.mat[i, "Chromosome"]) &
					((probe.dat[, 2] >= glad.mat[i, "Start.bp"]) &
					(probe.dat[, 2] <= glad.mat[i, "End.bp"]))
	}
	
	if (verbose) {
		print(paste("GetArrayAsSegs .... processing", nrow(glad.mat), "segments"))
	}
	##  check number of total probes in each segment
	n.seg <- nrow(glad.mat)
	n.snp.probes <- foreach(i = 1:n.seg, .combine=c) %dopar% {
		sum(GetProbeIx(i, snp.d)) 
	}
	n.cn.probes <- foreach(i = 1:n.seg, .combine=c) %dopar% {
		sum(GetProbeIx(i, cn.d))
	}
		
	if (verbose) {
		print(paste(nrow(glad.mat), " segs with > 2 SNPs and > 2 CN probes", sep=""))
	}
	n.seg <- nrow(glad.mat)
	
	if (verbose) {
		print(paste("get_AS_segs_ n.seg = ",n.seg))
	}
	
	
	snp.res = foreach(i = 1:n.seg) %dopar% {
		snp.ix <- GetProbeIx(i, snp.d)
		raw.snp.seg <- snp.d[snp.ix, c(3, 4), drop=FALSE]
		if (verbose) {
			print(paste("Processing Segment: ", i, "SNP size = ", dim(raw.snp.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"]))
		}
		
		h.snp.gt.p <- snp.gt.p[snp.ix, , drop=FALSE]
		rownames(h.snp.gt.p) <- rownames(snp.d)[snp.ix]
		
		h.snp.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=snp.d[snp.ix, 2]) 
		
		if (impute.gt) {
			dbSNP.ix <- match(rownames(snp.d)[snp.ix], rownames(dbSNP.annot))
			h.snp.annot[["dbSNP"]] <- dbSNP.annot[dbSNP.ix, , drop=FALSE] 
		}
		
		list(h.snp.d=t(raw.snp.seg), h.snp.gt.p=h.snp.gt.p, h.snp.annot=h.snp.annot)
	}
	h.snp.d = lapply(snp.res, "[[", "h.snp.d")
	h.snp.gt.p = lapply(snp.res, "[[", "h.snp.gt.p")
	h.snp.annot = lapply(snp.res, "[[", "h.snp.annot")

	
	# CN Probes
#	for (i in 1:n.seg) {
	cn.res = foreach(i = 1:n.seg) %dopar% {
		cn.ix <- GetProbeIx(i, cn.d)
		raw.cn.seg <- cn.d[cn.ix, "Intensity", drop=FALSE]
		if (verbose) {
			print(paste("Processing Segment: ", i, "CN size = ", dim(raw.cn.seg)[1],
							"start = ", glad.mat[i, "Start.bp"], "end = ",
							glad.mat[i, "End.bp"], "Chromosome = ",
							glad.mat[i, "Chromosome"]))
		}
		h.cn.d <- t(raw.cn.seg)
		h.cn.annot <- list(chr=as.integer(glad.mat[i, "Chromosome"]), pos=cn.d[cn.ix, 2])
		list(h.cn.d=h.cn.d, h.cn.annot=h.cn.annot)
	}
	h.cn.d = lapply(cn.res, "[[", "h.cn.d")
	h.cn.annot = lapply(cn.res, "[[", "h.cn.annot")

	
	return(list(h.snp.d=h.snp.d, h.cn.d=h.cn.d, h.snp.gt.p=h.snp.gt.p, h.snp.annot=h.snp.annot,
		h.cn.annot=h.cn.annot))
}


CalibrateAsDat <- function(dat, clusters.fn, calb.type="d1", 
		clusters.file.parser=DefaultClustersFileParser,
		verbose=FALSE) {
	if (verbose) {
		print(paste("Loading the birdseed clusters file: ", clusters.fn))
	}
	
	clusters <- clusters.file.parser(clusters.fn, verbose=verbose)
	
	dat <- dat[intersect(rownames(dat), rownames(clusters)), ]
	clusters <- clusters[rownames(dat), ]
	calb <- CalibrateAsSnps(clusters, dat, calb.type)
	
	return(calb)
}

CalibrateAsSnps <- function(clusters, dat, calb.type) {
	a.bg <- clusters[, "BB.a"]
	b.bg <- clusters[, "AA.b"]
	
	a.1 <- clusters[, "AB.a"]
	b.1 <- clusters[, "AB.b"]
	
	a.2 <- clusters[, "AA.a"]
	b.2 <- clusters[, "BB.b"]
	
	a.d1 <- a.1 - a.bg
	a.d2 <- a.2 - a.1
	
	b.d1 <- b.1 - b.bg
	b.d2 <- b.2 - b.1
	
	if (calb.type == "avg") {
		a.scale <- (a.d1 + a.d2) / 2
		b.scale <- (b.d1 + b.d2) / 2
	} else if (calb.type == "d1") {
		a.scale <- a.d1
		b.scale <- b.d1
	} else if (calb.type == "d2") {
		a.scale <- a.d2
		b.scale <- b.d2
	}
	
	a.calb <- (dat[,"A"] - a.bg) / a.scale
	b.calb <- (dat[,"B"] - b.bg) / b.scale
	
	cdat <- cbind(a.calb, b.calb)
	rownames(cdat) <- rownames(dat)
	cdat[!is.finite(cdat)] <- NA 
	
	return(cdat)
}



