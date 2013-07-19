#' This script performs the basics of the allelic capseg
#'
#' @author Aaron McKenna <aaron@broadinstitute.org>
#' @license MIT License, see http://www.opensource.org/licenses/mit-license.html
#' @note current version was created on March 19th, 2012

library(optparse)
options(keep.source = TRUE)
options(error = quote({
  sink(file="error.txt");
  dump.frames();
  print(attr(last.dump,"error.message"));
  traceback();
  sink();
  q()}))
# require("genomic_plot.R")
option.list <- list(
               # we include the R dat files from
               make_option(c("--output.dir"),help="the output directory",default="./"),
               make_option(c("--probe.file"),help="the tumor exome coverage: lanes as rows, exome targets (or baits) as rows",default="blank"),
               make_option(c("--segment.file"),help="the list of targets we captured in sequencing",default="blank"),
               make_option(c("--coverage.file"),help="Any sex chromosomes; each with by tangent normalized by itself. Seperate a list with a comma, no spaces",default="blank"),
               make_option(c("--bam.name"),help="bam file",default="blank"),
               make_option(c("--source.directory"),help="where to include source code from",default="blank"),
               make_option(c("--bam.to.sample"),help="the bam file to sample name mapping",default="blank"),
               make_option(c("--seg.merge.thresh"),help="the threshold to merge segments: segments smaller than this size with be merged with the appropriate neighbor",default="blank"),
               make_option(c("--min.seg.size"),help="the probability threshold for Bayesian segment merging",default="blank")
               )
opt <- parse_args(OptionParser(option_list=option.list))

# not sure if these should be parameters
drop.x= TRUE
drop.y= TRUE
verbose= TRUE

# where to put our output
base.output.dir= opt$output.dir
capseg.probe.fn= opt$probe.file # "/xchip/cga2/bryanh/HAPSEG/hapseg_extreme/brain.mets/capseg.results/PB-PB0036-TM-NT-SM-2O41H-SM-2O41I.tsv"
capseg.seg.fn= opt$segment.file # "/xchip/cga4/home/peleg/METS_CapSeg77TM_NT/cbs/PB-PB0036-TM-NT-SM-2O41H-SM-2O41I.seg.txt"
germline.het.fn= opt$coverage.file # "/xchip/cga4/home/peleg/METS_germline/hetReadCounts/PB-PB0036-Tumor-SM-2O41H.cov"
bam.file = opt$bam.name
code.dir = opt$source.directory
sample.to.bam = opt$bam.to.sample
min.seg.size= opt$min.seg.size
seg.merge.thresh = opt$seg.merge.thresh

sample.table <- read.delim(sample.to.bam,stringsAsFactors=F)
print(summary(sample.table))

sample.table <- sample.table[!is.na(sample.table$BAM),]

print(summary(sample.table))
print(bam.file)
SID = sample.table[sample.table$BAM == basename(bam.file),"Sample"]
print(SID)
if (is.na(SID) | length(SID) < 1) {
    print("Unable to find the sample, please check the inputs")
    quit(status=1)
}

dddd = lapply(list.files(code.dir, pattern="\\.R$", full.name=T), source)

SetPlatformSpecificFuncs("WES")

library(doMC)
registerDoMC(1)

library(numDeriv)  ## needed for 'hessian()'

# lookup the sample name from the sample mapping file
# save.image("allelic.data")
AllelicCapseg( capseg.probe.fn, capseg.seg.fn, germline.het.fn, SID, base.output.dir, min.seg.size, drop.x, drop.y, seg.merge.thresh, verbose )

