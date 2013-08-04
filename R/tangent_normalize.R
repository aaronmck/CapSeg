#' This script processes raw data from the exome coverage tool and produces
#' segmentation calls and other data needed for the HapSeg Algorithm for
#' each of the samples seen in the data
#'
#' @author Aaron McKenna <aaron@broadinstitute.org>
#' @license MIT License, see http://www.opensource.org/licenses/mit-license.html
#' @note current version was created on March 19th, 2012

library(optparse)
# require("genomic_plot.R")
option.list <- list(
		# we include the R dat files from
               make_option(c("--sample.table"),help="the file containing the mapping individuals to their tumor and normal bams",default="blank"),
               make_option(c("--normal.lane.data"),help="the normal exome coverage: lanes as rows, exome targets (or baits) as rows",default="./"),
               make_option(c("--tumor.lane.data"),help="the tumor exome coverage: lanes as rows, exome targets (or baits) as rows",default="blank"),
               make_option(c("--target.list"),help="the list of targets we captured in sequencing",default="blank"),
               make_option(c("--script.dir"),help="where we can find the wesseg scripts - where you placed the checked out tool into",default="blank"),
               make_option(c("--normal.sample.bams"),help="the file containing the mapping of the sample names (for normals) to bam files",default="blank"),
               make_option(c("--tumor.sample.bams"),help="the file containing the mapping of the sample names (for tumors) to bam files",default="blank"),
               make_option(c("--output.location"),help="where to write output files to - the segmentation results plus any graphs",default="blank"),
               make_option(c("--tangent.database.location"),help="the directory of tangent planes to normalize against; this directory should contain only tangent planes",default="blank"),
               make_option(c("--output.tangent.database"),help="the directory where we put the output tangent data",default="blank"),
               make_option(c("--build"),help="are we running with hg18 and hg19",default="blank"),
               make_option(c("--analysis.set.name"),help="what was the name of the analysis set",default="blank"),
               make_option(c("--bait.factor"),help="the bait factor data file",default="blank"),
               make_option(c("--bam.file.listing"),help="the listing of bam files, by tumor and by normal",default="blank"),
               make_option(c("--signal.files"),help="the sample name to signal file",default="blank"),
               make_option(c("--use.histo.data"),help="should we use historical data",default="blank"),
               make_option(c("--debug"),help="dump lots of debugging data to the <output_dir>/debug directory",default="blank")
)
opt <- parse_args(OptionParser(option_list=option.list))

# get the libraries and R code we rely on
library("plyr")
library("corpcor")

# make sure that we don't check names when we load up matrices - this replaces special characters with dots, which causes
# a lot of problems
options(check.names=F)

# now source the path to some of our utilities
code.dir = opt$script.dir
source(paste(code.dir,"/R/data_loading.R",sep="")) #,verbose=T)
source(paste(code.dir,"/R/output_data.R",sep="")) #,verbose=T)
source(paste(code.dir,"/R/command_line_utils.R",sep="")) #,verbose=T)
source(paste(code.dir,"/R/genomic_plot.R",sep="")) #,verbose=T)
source(paste(code.dir,"/R/data_transformations.R",sep="")) #,verbose=T)

# post process the file
args <- post.process.commandline(opt)

# -------------------------------------------------------------------------------
# load the command line arguments into local variables (with some post-processing)
# -------------------------------------------------------------------------------
sample.table      			<- args$sample.table
normal.lane.data				<- args$normal.lane.data
tumor.lane.data					<- args$tumor.lane.data
target.list.file				<- args$target.list
try.use.cached.data				<- args$use.cache
cached.location 				<- args$cache.location
wes.working.dir 				<- args$script.dir
normal.lanes.to.samples.file	<- args$normal.sample.bams
tumor.lanes.to.samples.file		<- args$tumor.sample.bams
output.location					<- args$output.location
tangent.database.location		<- args$tangent.database.location
tangent.database.output		    <- args$output.tangent.database
build.version					<- args$build
analysis.set.name				<- args$analysis.set.name
by.lane 					    <- args$bylane
lane.parallel 					<- args$parallel != "blank"
bait.factor.file 				<- args$bait.factor
bam.file.listing                <- args$bam.file.listing
signal.files                    <- args$signal.files
use.histo.data                  <- toupper(args$use.histo.data) == "TRUE"
debug                           <- toupper(args$debug) == "TRUE"
sex.calls                       <- args$sex.calls
vcf.calls.file                  <- args$call.database


# some constants to use
removeBadBaitsAndLanes = T # TRUE
optimize.bf = F # optimize those baits!
calibrate.against.others = T


if (build.version == "hg19") {
  sex.chromosomes = c("X","Y") # this should get loaded based on the platform
} else if (build.version == "hg18") {
  sex.chromosomes = c("chrX","chrY") # this should get loaded based on the platform
} else if (build.version == "mm9") {
  sex.chromosomes = c("chrX","chrY") # this should get loaded based on the platform
} else {
  print(paste("Unable to work with build version provided:",build.version,"expected hg18, hg19, mm9"))
}
# our epsilon value - used to make sure we're not producing log(0) calls
epsilon <- .Machine$double.eps * 10^6
# options(error=dump.frames)

# create the output directory and the cache directory if needed, and setup some debug logging locations (used only if debug == T)
create.dir.if.missing(output.location)
if (debug) {
    debug.location = paste(output.location,"debug",sep="/")
    create.dir.if.missing(debug.location)
    #sink(paste(output.location,"debug","debugging_log.txt",sep="/"))
    save.image(paste(debug.location,".parameters.Rdata",sep="/")) # save off a copy of the parameters we used for the run
}

# load the tumor and the normal sample information
if (debug) print("Loading the sample tables from the file on disk...")

tumorSampleTable = read.table(tumor.lanes.to.samples.file,header=T,stringsAsFactors=F,sep="\t")
normalSampleTable = read.table(normal.lanes.to.samples.file,header=T,stringsAsFactors=F,sep="\t")

# get the baits (targets)
if (debug) print(paste("loading the baits from the file",target.list.file))
baits <- read.csv(target.list.file,header=TRUE,colClasses=c("character","character","integer","integer"))
baits <- data.frame(baits[!duplicated(baits$name),])

# load the big data - the csv files of coverage; let the user know how long this is taking
if (debug) print(paste("Starting to load the data at",format(Sys.time(), "%a %b %d %H:%M:%S %Y")))

tumor.data <- load.exome.data(tumor.lane.data)
normal.data <- load.exome.data(normal.lane.data)

if (debug) print(paste("Intersecting the normal and tumor bait lists, normal data has",nrow(normal.data),"rows, tumor data has",nrow(tumor.data),"rows"))

# cut the rownames to the intersect of what's available in the normal and the tumor
target.intersect <- intersect(rownames(normal.data),rownames(tumor.data))
if (debug) print(paste("intersection of the tumor and normal has ",length(target.intersect),"rows"))

normal.data <- intersect.tumor.normal.targets(normal.data,target.intersect)
tumor.data <- intersect.tumor.normal.targets(tumor.data,target.intersect)
baits.filtered <- baits[is.element(baits$name,target.intersect),]

# quality control the normal coverage
qc.report.directory = paste(output.location,"/normal_tissue_qc_files/",sep="")
create.dir.if.missing(qc.report.directory)
normal.data = qc.normal.samples(normal.data, tumor.data, baits.filtered, qc.report.directory, allowed.normal.dev=0.3, allowed.normal.arm.dev=0.5)

# create sex calls (guesses) for each individual
sex.calls = create.sex.assignments(tumor.data,normal.data,sample.table,baits.filtered,"Y","X")
# and write this table out
write.table(sex.calls,file=paste(output.location,"sex.calls.txt",sep="/"),sep="\t",quote=F,row.names=F)

# load up our bait factor
bait.factor <- read.delim(bait.factor.file)
bait.factor[bait.factor[,1]<=0,2] = epsilon
bait.factor = bait.factor[is.element(rownames(bait.factor),rownames(tumor.data)),]

# do the non-sex based normalization
# mean center the tumor and normal samples
tumor.data = sweep(tumor.data,2,apply(tumor.data,2,mean),"/")

# split into females and males
male.tumors <- tumor.data[,is.element(colnames(tumor.data),sex.calls$sample[!as.logical(sex.calls$is.female)])]
female.tumors <- tumor.data[,is.element(colnames(tumor.data),sex.calls$sample[as.logical(sex.calls$is.female)])]
male.normals <- normal.data[,is.element(colnames(normal.data),sex.calls$sample[!as.logical(sex.calls$is.female)])]
female.normals <- normal.data[,is.element(colnames(normal.data),sex.calls$sample[as.logical(sex.calls$is.female)])]

male.calibrated <- calibrate.and.pi.tumors(male.tumors,male.normals)
female.calibrated <- calibrate.and.pi.tumors(female.tumors,female.normals)

# put together the calibrated data
calibrated.tumors <- cbind(male.calibrated,female.calibrated)

# save off the data to a Rdata object for later
save.off.processed.data(log2(normal.data),log2(tumor.data),calibrated.tumors,baits,paste(tangent.database.output,build.version,sep="/"),analysis.set.name,build.version)

# output the raw data and plots for each sample
output.and.plot.data(calibrated.tumors,tumor.data,baits,output.location,signal.files)

# output metrics for each sample
create.sample.metric.files(calibrated.tumors,log.tumors,output.location)
