#' This script performs the basics of the allelic capseg
#'
#' @author Aaron McKenna <aaron@broadinstitute.org>
#' @license MIT License, see http://www.opensource.org/licenses/mit-license.html
#' @note current version was created on March 19th, 2012

library(optparse)
# require("genomic_plot.R")
option.list <- list(
               # we include the R dat files from
               make_option(c("--output.dir"),help="the output directory",default="./"),
               make_option(c("--probe.file"),help="the tumor exome coverage: lanes as rows, exome targets (or baits) as rows",default="blank"),
               make_option(c("--segment.file"),help="the list of targets we captured in sequencing",default="blank"),
               make_option(c("--coverage.file"),help="Any sex chromosomes; each with by tangent normalized by itself. Seperate a list with a comma, no spaces",default="blank")
               make_option(c("--sample.name"),help="Sample name",default="blank")
)
opt <- parse_args(OptionParser(option_list=option.list))


drop.x= TRUE
drop.y= TRUE
min.seg.size= 10
verbose= TRUE
base.output.dir= "AllelicCapseg_results/Brain_mets_AllelicCapseg_03.05"
capseg.probe.fn= "/xchip/cga2/bryanh/HAPSEG/hapseg_extreme/brain.mets/capseg.results/PB-PB0036-TM-NT-SM-2O41H-SM-2O41I.tsv"
capseg.seg.fn= "/xchip/cga4/home/peleg/METS_CapSeg77TM_NT/cbs/PB-PB0036-TM-NT-SM-2O41H-SM-2O41I.seg.txt"
germline.het.fn= "/xchip/cga4/home/peleg/METS_germline/hetReadCounts/PB-PB0036-Tumor-SM-2O41H.cov"
SID= "PB-PB0036-TM-NT-SM-2O41H-SM-2O41I"
#library(HAPSEG)

CODE_DIR = "~scarter/HAPSEG/"
lapply(list.files(CODE_DIR, pattern="\\.R$", full.name=T), source)

SetPlatformSpecificFuncs("WES" )

library(doMC)
registerDoMC(1)

library(numDeriv)  ## needed for 'hessian()'


AllelicCapseg( capseg.probe.fn, capseg.seg.fn, germline.het.fn, SID, base.output.dir, min.seg.size, drop.x, drop.y, verbose )
