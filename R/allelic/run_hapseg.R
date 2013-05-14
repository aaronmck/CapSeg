## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

RunHapSeg <- function(plate.name, array.name, seg.fn, snp.fn, 
                      genome.build, RESULTS.DIR, platform, use.pop, impute.gt,
                      plot.segfit, merge.small, merge.close, min.seg.size,
                      normal, out.p, seg.merge.thresh, use.normal, adj.atten,
                      phased.bgl.dir, force.diploid=normal, drop.x=FALSE,
                      drop.y=TRUE, calls.fn=NULL, mn.sample=NULL,
                      out.file=NULL, calibrate.data=FALSE,
                      clusters.fn=NULL, prev.theta.fn=NULL, 
                      snp.file.parser=DefaultSnpFileParser,
                      clusters.file.parser=DefaultClustersFileParser,
                      verbose=FALSE) {

  if (!CheckGenomeBuild(genome.build)) {
    stop("Unsupported genome build: ", genome.build)
  }

  ## Note that we're not using the built in R tmpdir() due to
  ## potential space issues
  tmp.dir <- CreateTmpDir(RESULTS.DIR)
  on.exit(try(unlink(tmp.dir, recursive=TRUE), silent=TRUE), add=TRUE)

  iams.res <- InitAndMergeSmall(array.name, genome.build, use.pop, use.normal,
                                normal, impute.gt, adj.atten, platform,
                                seg.fn, snp.fn, drop.x, drop.y,
                                calls.fn, mn.sample, min.seg.size, merge.close,
                                out.p, merge.small, force.diploid,
                                calibrate.data, clusters.fn,
                                prev.theta.fn, snp.file.parser,
                                clusters.file.parser, verbose)
  
						
  if (merge.close == TRUE) {
    ## use fit of error model to merge segments
    mcaf.res <- MergeCloseAndFit(iams.res, out.p, seg.merge.thresh,
                                 impute.gt, platform, force.diploid=force.diploid,
                                 verbose=verbose)
  } else {
    if (verbose) {
      print("Merged Loci: nulled ")
    }
    seg.dat[["merged.loci"]] <- NA
    seg.dat[["final.merge.prob"]] <- matrix(NA, ncol=2,
                                            nrow=length(iams.res[["h.d"]]) - 1)
    seg.dat[["seg.expected.phase"]] <- iams.res[["em.fit"]][["seg.expected.phase"]]
    mcaf.res <- iams.res
    mcaf.res[["seg.dat"]] <- seg.dat
  }

  if (impute.gt == TRUE) {    
    gt.res <- ImputeGt(mcaf.res, platform, tmp.dir, plate.name, out.p,
                       phased.bgl.dir, verbose=verbose)
    gt.res <- PostImputeSegFit(gt.res, out.p, platform,
                               force.diploid=force.diploid, verbose=verbose)
  } else {
    gt.res <- mcaf.res
  }

  seg.dat <- FinishHapsegAndSave(gt.res, plate.name, normal, array.name,
                                 RESULTS.DIR, platform, out.file)

  ## final step - plot all of the data
  if (plot.segfit) {
    DoPlots(seg.dat, RESULTS.DIR, platform, verbose=verbose)
  } 

  return(TRUE)
}

