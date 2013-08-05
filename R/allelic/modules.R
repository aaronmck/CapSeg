
InitAndMergeSmall = function(chr.truncate=NULL, cached=F, fn="iams.res.rds") {
	print("Init and Merge Small")	
	
	LoadCached(iams.res <- ExtractProbeDat(array.name, capseg.sample.name, genome.build, use.pop, use.normal,
					normal, impute.gt, adj.atten, platform, 
					seg.fn, snp.fn, cn.fn, capseg.probe.fn, capseg.seg.fn,
					drop.x, drop.y, calls.fn, mn.sample,
					calibrate.data=calibrate.data,
					clusters.fn=clusters.fn, verbose=verbose,
					snp.file.parser=snp.file.parser,
					cn.file.parser=cn.file.parser,
					clusters.file.parser=clusters.file.parser), 
			cached=cached, res.fn=file.path(results.dir, "extracted.probe.dat.rds"), mod.name="Init and Merge Small")
	
#	saveRDS(iams.res, file=file.path(results.dir, "extracted.probe.dat.rds"))
#	iams.res = readRDS(file.path(results.dir, "extracted.probe.dat.rds"))
	
	TruncateData = function(res, chr) {
		trunc.elem = c("h.snp.d", "h.cn.d", "h.capseg.d", "h.snp.gt.p", "h.snp.annot", "h.cn.annot", "h.capseg.annot")
		kill.idx = c()
		for (i in 1:length(iams.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]])) {
			if (res[["as.res"]][["h.seg.dat"]][['h.snp.annot']][[i]][["chr"]] != chr) kill.idx = c(kill.idx, i)
		}
		for (elem in trunc.elem) {
			res[["as.res"]][["h.seg.dat"]][[elem]][kill.idx] <- NULL
		}
		return(res)
	}
	
	DownSample = function(res, frac=1/10, min = 50 ) {
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
	
	if (!is.null(chr.truncate) ) {
		iams.res = TruncateData(iams.res, chr.truncate)
	}
	
	if (merge.small) {
		merge.res <- JoinSmallSegsExtreme(iams.res[["as.res"]], min.seg.size, verbose=verbose)   	  
		iams.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]] <- merge.res[["h.snp.d"]]
		iams.res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]] <- merge.res[["h.snp.gt.p"]]
		iams.res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]] <- merge.res[["h.snp.annot"]]
		iams.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]] <- merge.res[["h.cn.d"]]
		iams.res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]] <- merge.res[["h.cn.annot"]]
		iams.res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]] <- merge.res[["h.capseg.d"]]
		iams.res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]] <- merge.res[["h.capseg.annot"]]
	}
	
	
	iams.res[["em.fit"]][["theta"]] <- InitThetaExtreme(array.name, prev.theta.fn, iams.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]], iams.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]], verbose)
	iams.res[["use.eps"]] <- ifelse((merge.close || impute.gt), 5e-3, 1e-4)
	

	iams.res[["em.fit"]] <- HscrSegFitExtreme(iams.res[["as.res"]][["h.seg.dat"]], iams.res[["em.fit"]][["theta"]], eps=iams.res[["use.eps"]], out.p=out.p, 
			force.diploid=force.diploid, verbose=verbose)
	
	saveRDS(iams.res, file=file.path(results.dir, fn))
	return(iams.res)
}

MergeCloseAndFit = function(iams.res, cached=F, fn="mcaf.res.rds") {
	
	print("Merge Close and Fit")
	CheckCaching(cached, file.path(results.dir, "image_at_mcaf_res.rda"), "Merge close and Fit")
	
	if (merge.close == TRUE) {
		mcaf.res = iams.res
		## use fit of error model to merge segments
		
		mrg.res <- JoinCloseSegsExtreme(h.d = list(snp=mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]], cn=mcaf.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]]), 
				h.snp.gt.p=mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]], 
				h.probe.annot=list(snp=mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]], cn=mcaf.res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]]), 
				theta=mcaf.res[["em.fit"]][["theta"]], force.diploid, out.p,merge.thresh=seg.merge.thresh, verbose=verbose)
		
		mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]] <- mrg.res[["h.d"]][["snp"]]
		mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]] <- mrg.res[["h.snp.gt.p"]]
		mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]] <- mrg.res[["h.probe.annot"]][["snp"]]
		mcaf.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]] <- mrg.res[["h.d"]][["cn"]]
		mcaf.res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]] <- mrg.res[["h.probe.annot"]][["cn"]]
		mcaf.res[["seg.dat"]][["merged.loci"]] <- mrg.res[["merged.loci"]]
		mcaf.res[["seg.dat"]][["final.merge.prob"]] <- mrg.res[["final.merge.prob"]]
		
		## update fit with merged segs 
		mcaf.res[["use.eps"]] <- ifelse(impute.gt, 1e-2, 1e-4)
		
		mcaf.res[["em.fit"]] <- HscrSegFit(mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]], mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]], 
				mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]], mcaf.res[["em.fit"]][["theta"]], mcaf.res[["as.res"]][["h.seg.dat"]][["segmentation.information"]],
				eps=mcaf.res[["use.eps"]], out.p=out.p, force.diploid=force.diploid, verbose=verbose)
		
		mcaf.res[["seg.dat"]][["seg.expected.phase"]] <- mcaf.res[["em.fit"]][["seg.expected.phase"]]
		
		
	} else {
		mcaf.res = iams.res
		mcaf.res[["seg.dat"]][["merged.loci"]] <- NA
		mcaf.res[["seg.dat"]][["final.merge.prob"]] <- matrix(NA, ncol=2,nrow=length(iams.res[["h.snp.d"]]) - 1)
		mcaf.res[["seg.dat"]][["seg.expected.phase"]] <- iams.res[["em.fit"]][["seg.expected.phase"]]
		
	}
	
	saveRDS(mcaf.res, file=file.path(results.dir, "mcaf.res.rds"))
	return(mcaf.res)
}

ImputeGT = function(cached=F, fn="gt.res.rds") {
	
		
		ps.res <- PhaseSnps(mcaf.res[["em.fit"]][["snp.clust.p"]], mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]], 
				mcaf.res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]], platform, tmp.dir, plate.name, phased.bgl.dir, verbose=verbose)
		## G_hat_1
		gt.res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]] <- ps.res[["h.snp.gt.p"]]  
		gt.res[["seg.dat"]][["h.switch.ix"]] <- ps.res[["h.switch.ix"]]
		
		## finalize fit with phased SNPs  
		gt.res[["em.fit"]] <- HscrSegFit(gt.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]], gt.res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]], 
				gt.res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]], gt.res[["em.fit"]][["theta"]], eps=1e-4, out.p=out.p, force.diploid=force.diploid, verbose=verbose)
		gt.res[["em.fit"]][["theta"]][["loglik"]] <- gt.res[["em.fit"]][["loglik"]]
		gt.res[["em.fit"]][["theta"]][["p.snp.cond.out"]] <- 1 / 25
		
	return(gt.res) 
}


