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
		make_option(c("--normal.lane.data"),help="the normal exome coverage: lanes as rows, exome targets (or baits) as rows",default="./"),
		make_option(c("--tumor.lane.data"),help="the tumor exome coverage: lanes as rows, exome targets (or baits) as rows",default="blank"),
		make_option(c("--target.list"),help="the list of targets we captured in sequencing",default="blank"),
		make_option(c("--use.cache"),help="should we use the cached data from a previous run (if available)?",default="blank"),
		make_option(c("--cache.location"),help="where we should look for cached data (if use_cached_data is set), and where we should save data to post-processing",default="blank"),
		make_option(c("--script.dir"),help="where we can find the wesseg scripts - where you placed the checked out tool into",default="blank"),
		make_option(c("--normal.sample.to.lanes.file"),help="the file containing the mapping of the read groups to the sample names (for normals)",default="blank"),
		make_option(c("--tumor.sample.to.lanes.file"),help="the file containing the mapping of the read groups to the sample names (for tumors)",default="blank"),
                make_option(c("--normal.to.bam"),help="the file matching the normal bam to the sample",default="blank"),
		make_option(c("--tumor.to.bam"),help="the file matching for tumor bam to the sample",default="blank"),
		make_option(c("--output.location"),help="where to write output files to - the segmentation results plus any graphs",default="blank"),
		make_option(c("--tangent.database.location"),help="the directory of tangent planes to normalize against; this directory should contain only tangent planes",default="blank"),
		make_option(c("--build"),help="are we running with hg18 and hg19",default="blank"),
		make_option(c("--analysis.set.name"),help="what was the name of the analysis set",default="blank"),
		make_option(c("--bylane"),help="is the data coming in by lane? (if not it should be by sample)",default="blank"),
		make_option(c("--parallel"),help="should we merge lanes to samples in parallel",default="blank"),
                make_option(c("--bait.factor"),help="the bait factor data file",default="blank"),
                make_option(c("--bam.file.listing"),help="the listing of bam files, by tumor and by normal",default="blank"),
                make_option(c("--signal.files"),help="the sample name to signal file",default="blank"),
                make_option(c("--histo.data"),help="should we use historical data?",default="blank")
               make_option(c("--debug"),help="dumb lots of debugging data to the <output_dir>/debug directory",default="blank")
)
opt <- parse_args(OptionParser(option_list=option.list))

# get the libraries and R code we rely on
library("plyr")
library("entropy")
library("corpcor")

# make sure that we don't check names when we load up matrices - this replaces special characters with dots, which causes
# a lot of problems
options(check.names=F)

# now source the path to some of our utilities
wes.working.dir = opt$script.dir
source(paste(wes.working.dir,"/R/extract_bait_data.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/data_loading.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/output_data.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/command_line_utils.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/genomic_plot.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/debug_functions.R",sep="")) #,verbose=T)

# post process the file
args <- post.process.commandline(opt)

# -------------------------------------------------------------------------------
# load the command line arguments into local variables (with some post-processing)
# -------------------------------------------------------------------------------
normal.lane.data				<- args$normal.lane.data
tumor.lane.data					<- args$tumor.lane.data
target.list.file				<- args$target.list
try.use.cached.data				<- args$use.cache
cached.location 				<- args$cache.location
wes.working.dir 				<- args$script.dir
normal.lanes.to.samples.file	                <- args$normal.sample.to.lanes.file
tumor.lanes.to.samples.file		        <- args$tumor.sample.to.lanes.file
output.location					<- args$output.location
tangent.database.location		        <- args$tangent.database.location
build.version					<- args$build
analysis.set.name				<- args$analysis.set.name
by.lane 					<- args$bylane
lane.parallel 					<- args$parallel != "blank"
bait.factor.file 				<- args$bait.factor
bam.file.listing                                <- args$bam.file.listing
tumor.to.bam                                    <- args$tumor.to.bam
normal.to.bam                                   <- args$normal.to.bam
signal.files                                    <- args$signal.files
histo.data                                      <- toupper(args$histo.data) == "TRUE"
debug                                           <- toupper(args$debug) == "TRUE"


# some constants to use
removeBadBaitsAndLanes = T # TRUE
optimize.bf = F # optimize those baits!
calibrate.against.others = T

# our epsilon value - used to make sure we're not producing log(0) calls
epsilon <- .Machine$double.eps * 10^6
# options(error=dump.frames)



# create the output directory and the cache directory if needed, and setup some debug logging locations (used only if debug == T)
dir.create(output.location,recursive=T)
dir.create(cached.location,recursive=T)
debug.location = paste(output.location,"debug",sep="/")

if (debug) {
    dir.create(debug.location,recursive=T)
    sink(paste(output.location,"debug","debugging_log.txt",sep="/"))
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

# load the big data data - the csv files of coverage; let the user know how long this is taking
if (debug) print(paste("Starting to load the data at",format(Sys.time(), "%a %b %d %H:%M:%S %Y")))

# create the output directory and the cache directory if needed
dir.create(output.location,recursive=T)
dir.create(cached.location,recursive=T)

tumor.data <- load.exome.data(tumor.lane.data)
normal.data <- load.exome.data(normal.lane.data)

if (debug) print(paste("Intersecting the normal and tumor bait lists, normal data has",nrow(normal.data),"rows, tumor data has",nrow(tumor.data),"rows"))

# cut the data to the intersect of the data
target.intersect <- intersect(rownames(normal.data),rownames(tumor.data))
if (debug) print(paste("intersection of the tumor and normal has ",length(target.intersect),"rows"))

normal.data <- normalize.tumor.normal.targets(normal.data,target.intersect)
tumor.data <- normalize.tumor.normal.targets(tumor.data,target.intersect)

# load up our bait factor
bait.factor <- read.delim(bait.factor.file)
bait.factor[bait.factor[,1]<=0,2] = epsilon
bait.factor = bait.factor[is.element(rownames(bait.factor),rownames(tumor.data)),]

td.mean = apply(tumor.data,2,mean)
tumor.data = tumor.data / td.mean
tumor.data[abs(tumor.data) < epsilon] = epsilon
log.tumors = data.frame(log2(tumor.data))

nd.mean = apply(normal.data,2,mean)
normal.data = normal.data / nd.mean
normal.data[abs(normal.data) < epsilon] = epsilon
log.normals = data.frame(log2(normal.data))

# do the initial block normalization
bgs.center = rep(mean(colMeans(log.normals)),nrow(log.normals)) # rep(0.0,nrow(normal.data)) # log.normals[,ncol(log.normals)] # rep(0.0,nrow(normal.data))
names(bgs.center) <- rownames(log.normals)
log.normals = data.frame(log.normals[,1:ncol(log.normals)] - bgs.center)
log.tumors = data.frame(log.tumors - bgs.center)

print("Data loaded and means and log values calculated...")
#load up the whole data set into the tangent normalization process, and calibrate each tumor against the matrix
if (histo.data) {
    print("loading the historical data...")
    log.normals = additional.normals(tangent.database.location,log.normals,build.version,analysis.set.name)
}

target.intersect <- intersect(rownames(log.normals),rownames(log.tumors))
print(dim(log.normals))
print(dim(log.tumors))
log.normals <- normalize.tumor.normal.targets(log.normals,target.intersect)
log.tumors <- normalize.tumor.normal.targets(log.tumors,target.intersect)

print("About to perform the SVD (PI)")
pseudo.inverse.norm  <- pseudoinverse(data.matrix(log.normals))
calibrated.tumors <- calibrate.tumors(data.matrix(log.tumors),data.matrix(pseudo.inverse.norm),data.matrix(log.normals),bgs.center,first=TRUE)
colnames(calibrated.tumors) <- colnames(tumor.data)

# save off the data to a Rdata object for later
save.off.processed.data(log.normals,log.tumors,calibrated.tumors,baits,paste(tangent.database.location,build.version,sep="/"),analysis.set.name,build.version)

# output the raw data and plots for each sample
output.and.plot.data(calibrated.tumors,tumor.data,baits,output.location,signal.files)

# output metrics for each sample
create.sample.metric.files(calibrated.tumors,log.tumors,output.location)
