## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

PostProcessPlate <- function(plate.name, hapseg.files, 
                             plate.res.dir, calls.fn=NULL, verbose=FALSE) {
  if (!file.exists(plate.res.dir)) {
    dir.create(plate.res.dir, recursive=TRUE)
  }
  
  pdf.fn <- file.path(plate.res.dir,
                      paste(plate.name, ".segplots.pdf", sep=""))
  PlotEmSegFits(plate.name, hapseg.files, pdf.fn)

  fns <- file.path(plate.res.dir,
                   paste(plate.name,
                         c(".stats_summary.pdf",
                           ".SEG-MERGE_summary.pdf",
                           ".stats.RData"), sep=""))
  
  SummarizePlateStats(plate.name, hapseg.files, fns[1],
                      fns[2], fns[3], calls.fn=calls.fn,
                      verbose=verbose)

  return(TRUE)
}

  
