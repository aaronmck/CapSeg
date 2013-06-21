## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


AbsolutePostProcess <- function(res, seg.dat) {
	# res <- iams.res

	# seg.dat <- res[["seg.dat"]]
	add <- data.frame(f=res[["capture.em.fit"]][["wes.f"]][, "f.hat"], tau=res[["capture.em.fit"]][["delta.tau"]][, "tau"], sigma.tau=res[["capture.em.fit"]][["cap.e.mu"]][, "sigma3"], mu.minor=res[["capture.em.fit"]][["cap.e.mu"]][, 'mu1'], sigma.minor=res[["capture.em.fit"]][["cap.e.mu"]][,"sigma1"], mu.major=res[["capture.em.fit"]][["cap.e.mu"]][, "mu2"], sigma.major=res[["capture.em.fit"]][["cap.e.mu"]][, "sigma2"], stringsAsFactors=F)
	out <- cbind(seg.dat, add)
	return(out)
}

ArrayWesConcordanceStat <- function(res) {
	# res <- iams.res

	require(mnormt)
	unatten.snp.mu <- t(apply(res[["capture.em.fit"]][["delta.tau"]], 1, function(x) AffyGetMeans(x[1], x[2])))

	n.segs = length(res[["capture.em.fit"]][["cap.e.mu"]][,"mu1"])
	snp.mu <- c(unatten.snp.mu[,1], unatten.snp.mu[,2], unatten.snp.mu[,3] )
	wes.mu <- c(AffyInvAtten(res[["capture.em.fit"]][["cap.e.mu"]][,"mu1"], res[["capture.em.fit"]][["theta"]][["at.capseg"]]),
			AffyInvAtten(res[["capture.em.fit"]][["cap.e.mu"]][,"mu2"], res[["capture.em.fit"]][["theta"]][["at.capseg"]]),
			AffyInvAtten(res[["capture.em.fit"]][["cap.e.mu"]][,"mu3"], res[["capture.em.fit"]][["theta"]][["at.capseg"]]) )

	snp.sigma = rep((res[["capture.em.fit"]][["delta.tau.sd"]][, "tau"]^2 + res[["capture.em.fit"]][["delta.tau.sd"]][, "delta"]^2)^(1/2) / 2, 3)

	wes.sigma = c( res[["capture.em.fit"]][["cap.e.mu"]][,"sigma1"], res[["capture.em.fit"]][["cap.e.mu"]][,"sigma2"], res[["capture.em.fit"]][["cap.e.mu"]][,"sigma3"] )

	idx = complete.cases(snp.mu) & complete.cases(wes.mu) & complete.cases(snp.sigma) & complete.cases(wes.sigma)
	snp.mu <- snp.mu[idx]
	wes.mu <- wes.mu[idx]
	snp.sigma <- snp.sigma[idx]
	wes.sigma <- wes.sigma[idx]

	stat <- sum(sapply(seq_along(wes.mu), function(i) {
		# i = 7
		center <- c(snp.mu[i], wes.mu[i])
		eval.pt <- rep( 0.5 * sum(center), 2)
		dmnorm(eval.pt, mean=center, varcov=matrix(c(snp.sigma[i], 0, 0, wes.sigma[i]), ncol=2), log=T )
		}))

	return(stat)
}

CalcConf = function(mu, sigma, conf=.95) qnorm(c((1-conf)/2, (1+conf)/2), mean=mu, sd=sigma)

GetGCContent = function(seg.dat, verbose=FALSE) {

	if (genome.build == "hg19") {
		require(BSgenome.Hsapiens.UCSC.hg19)
	} else if (genome.build == "hg18") {
		require(BSgenome.Hsapiens.UCSC.hg18)
	} else {
		print(paste(genome.build, "not supported for GC content.  Not returning any info"))
		return(NULL)
	}
	gc.perc = foreach(r=1:nrow(seg.dat), .combine=c) %dopar% {
				if (verbose) loopStatus(r, step=100)
				chr.name = gsub("24", "Y", gsub("23", "X", paste("chr", seg.dat$Chromosome[r], sep="")))
				d = DNAString(Hsapiens[[chr.name]])[ seg.dat$Start.bp[r]:seg.dat$End.bp[r] ]
				gc.count = letterFrequency(d, letters="GC")
				gc.count/length(d)
			}
	return(gc.perc)
}

LoadCached = function(statement, overwrite, res.fn, mod.name) {
	if (!overwrite) {
		if (file.exists(res.fn)) {
			print(paste(mod.name, " result already computed.  Returning cached version."))
			res = readRDS(res.fn)
			return(res)
		} else {
			print(paste("No cached result found for ", mod.name, ".  Recomputing and saving."))
			statement
		}
	} else {
		print(paste("Not looking for cached result for ", mod.name, ".  Computing from scratch and saving."))
		statement
	}

}

TruncateData = function(res, chr) {
	if (is.null(chr)) {
		return(res)
	}
	trunc.list = c("h.snp.d", "h.cn.d", "h.capseg.d", "h.snp.gt.p", "h.snp.annot", "h.cn.annot", "h.capseg.annot", "gh.wes.allele.d", "gh.wes.allele.annot")
	keep.idx = c()
	for (i in 1:length(res[["as.res"]][["h.seg.dat"]][["h.snp.d"]])) {
		if (res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["chr"]] == chr) keep.idx = c(keep.idx, i)
	}
	for (elem in trunc.list) {
		res[["as.res"]][["h.seg.dat"]][[elem]] <- res[["as.res"]][["h.seg.dat"]][[elem]][keep.idx]
	}

	res[["seg.dat"]][["seg.info"]] = res[["seg.dat"]][["seg.info"]][res[["seg.dat"]][["seg.info"]]$Chromosome == chr, ]

#	res[["capture.em.fit"]][["e.mu"]] <- res[["capture.em.fit"]][["e.mu"]][keep.idx, ]
#	res[["capture.em.fit"]][["wes.f"]]  <- res[["capture.em.fit"]][["wes.f"]][keep.idx]
#	res[["capture.em.fit"]][["het.phase.log.p"]] <- res[["capture.em.fit"]][["het.phase.log.p"]][keep.idx]
#	res[["capture.em.fit"]][["snp.clust.p"]] <- res[["capture.em.fit"]][["snp.clust.p"]][keep.idx]

	return(res)
}

DownSample = function(res, frac=1/10, min = 50 ) {
	## FIXME:  This doesn't work still
	out.res = res
	for (i in 1:length(res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]]) ) {loopStatus(i, step=1)
#		i = 1
		if(length(res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["pos"]]) * frac < min |
				length(res[["as.res"]][["h.seg.dat"]][['h.cn.annot']][[i]][["pos"]]) * frac < min) {
			# Segment already small
			next
		} else { # Seg down sampled
			snp.idx = 1:length(res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["pos"]])
			snp.keep.idx = sample(snp.idx, floor(frac * length(snp.idx)), replace = F)
			cn.idx = 1:length(res[["as.res"]][["h.seg.dat"]][['h.cn.annot']][[i]][["pos"]])
			cn.keep.idx = sample(cn.idx, floor(frac * length(cn.idx)), replace = F)
			out.res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["pos"]] = res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["pos"]][snp.keep.idx]
			out.res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["dbSNP"]] = res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["dbSNP"]][snp.keep.idx, , drop=F]
			out.res[["as.res"]][["h.seg.dat"]][['h.cn.annot']][[i]][["pos"]] = res[["as.res"]][["h.seg.dat"]][['h.cn.annot']][[i]][["pos"]][cn.keep.idx]
			out.res[["as.res"]][["h.seg.dat"]][['h.snp.d']][[i]] = res[["as.res"]][["h.seg.dat"]][['h.snp.d']][[i]][ , snp.keep.idx, drop=F]
			out.res[["as.res"]][["h.seg.dat"]][['h.cn.d']][[i]] = res[["as.res"]][["h.seg.dat"]][['h.cn.d']][[i]][ , cn.keep.idx, drop=F]
			out.res[["as.res"]][["h.seg.dat"]][['h.snp.gt.p']][[i]] = res[["as.res"]][["h.seg.dat"]][['h.snp.gt.p']][[i]][snp.keep.idx, , drop=F]
		}

	}
	return(out.res)
}



PrintHapSegStartMessage <- function(chars.per.line=73) {
  cat(rep(" ", chars.per.line), "\n", sep="")
  cat(rep("-", chars.per.line), "\n", sep="")
  cat("Data is loaded...running HapSeg:\n")
  cat(rep("-", chars.per.line), "\n", sep="")
}

LogAdd <- function(x) {

  ##  Calculates log(sum(exp(x)))  without "leaving" log space
  if (is.vector(x)) {
    mix <- which.max(x)
    max <- x[mix]
    res <- max + log(sum(exp(x - max )))
  }

  if (is.matrix(x)) {
    mv <- apply(x, 1, max)
    res <- mv + log(rowSums(exp(x - mv)))
  }

  return(res)
}

CheckGenomeBuild <- function(genome.build) {
  return(genome.build %in% c("hg18", "hg19"))
}

CreateTmpDir <- function(RESULTS.DIR) {
  tmp.dir <- file.path(RESULTS.DIR, "tmp")
  dir.create(tmp.dir, recursive=TRUE, showWarnings=FALSE)
  return(tmp.dir)
}

ReadCol <- function(fn, col.name, save.rownames=FALSE) {
  if (!file.exists(fn)) {
    stop("File does not exist: ", fn)
  }

  header <- read.delim(fn, nrow=1, as.is=TRUE, header=FALSE)
  if (nrow(header) < 1) {
    stop("Unable to read header information for ", fn)
  }

  col.idx <- match(col.name, header[1, ])
  if (is.na(col.idx)) {
    stop("Column ", col.name, " does not exist in file ", fn)
  }
  if (save.rownames) {
    col.idx <- c(1, col.idx)
    row.names <- 1
  } else {
    row.names <- NULL
  }
  col.idx.str <- paste(col.idx, collapse=",")
  cut.str <- paste("cut -f", col.idx.str, " ", fn, sep="")

  df <- read.delim(pipe(cut.str), row.names=row.names, as.is=TRUE)

  return(df)
}

GetChromPlotDir <- function(chrom, RESULTS.DIR) {
  return(file.path(RESULTS.DIR, paste("chr", chrom, sep="")))
}

