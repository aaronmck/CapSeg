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
               make_option(c("--output.tangent.database"),help="the directory where we put the output tangent data",default="blank"),
               make_option(c("--build"),help="are we running with hg18 and hg19",default="blank"),
               make_option(c("--analysis.set.name"),help="what was the name of the analysis set",default="blank"),
               make_option(c("--bylane"),help="is the data coming in by lane? (if not it should be by sample)",default="blank"),
               make_option(c("--parallel"),help="should we merge lanes to samples in parallel",default="blank"),
               make_option(c("--bait.factor"),help="the bait factor data file",default="blank"),
               make_option(c("--bam.file.listing"),help="the listing of bam files, by tumor and by normal",default="blank"),
               make_option(c("--signal.files"),help="the sample name to signal file",default="blank"),
               make_option(c("--use.histo.data"),help="should we use historical data",default="blank"),
               make_option(c("--debug"),help="dump lots of debugging data to the <output_dir>/debug directory",default="blank"),
               make_option(c("--sex.calls"),help="the sex chromosome assignment",default="blank")
)
opt <- parse_args(OptionParser(option_list=option.list))

# get the libraries and R code we rely on
library("plyr")
library("corpcor")

# make sure that we don't check names when we load up matrices - this replaces special characters with dots, which causes
# a lot of problems
options(check.names=F)

# now source the path to some of our utilities
wes.working.dir = opt$script.dir
source(paste(wes.working.dir,"/R/data_loading.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/output_data.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/command_line_utils.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/genomic_plot.R",sep="")) #,verbose=T)
source(paste(wes.working.dir,"/R/data_transformations.R",sep="")) #,verbose=T)

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
normal.lanes.to.samples.file	<- args$normal.sample.to.lanes.file
tumor.lanes.to.samples.file		<- args$tumor.sample.to.lanes.file
output.location					<- args$output.location
tangent.database.location		<- args$tangent.database.location
tangent.database.output		    <- args$output.tangent.database
build.version					<- args$build
analysis.set.name				<- args$analysis.set.name
by.lane 					    <- args$bylane
lane.parallel 					<- args$parallel != "blank"
bait.factor.file 				<- args$bait.factor
bam.file.listing                <- args$bam.file.listing
tumor.to.bam                    <- args$tumor.to.bam
normal.to.bam                   <- args$normal.to.bam
signal.files                    <- args$signal.files
use.histo.data                  <- toupper(args$use.histo.data) == "TRUE"
debug                           <- toupper(args$debug) == "TRUE"
sex.calls                       <- args$sex.calls
vcf.calls.file                  <- args$call.database

# some constants to use
removeBadBaitsAndLanes = T # TRUE
optimize.bf = F # optimize those baits!
calibrate.against.others = T
debug = TRUE
# our epsilon value - used to make sure we're not producing log(0) calls
epsilon <- .Machine$double.eps * 10^6
sex.chromosomes = c("X","Y")

# create the output directory and the cache directory if needed, and setup some debug logging locations (used only if debug == T)
dir.create(output.location,recursive=T)
dir.create(cached.location,recursive=T)
debug.location = paste(output.location,"debug",sep="/")

dir.create(debug.location,recursive=T)
save.image(paste(debug.location,".parameters.Rdata",sep="/")) # save off a copy of the parameters we used for the run

# load the tumor and the normal sample information
print("Loading the sample tables from the file on disk...")

tumorSampleTable = read.table(tumor.lanes.to.samples.file,header=T,stringsAsFactors=F,sep="\t")
normalSampleTable = read.table(normal.lanes.to.samples.file,header=T,stringsAsFactors=F,sep="\t")

# get the baits (targets)
if (debug) print(paste("loading the baits from the file",target.list.file))
baits <- read.csv(target.list.file,header=TRUE,colClasses=c("character","character","integer","integer"))
baits <- data.frame(baits[!duplicated(baits$name),]) # make sure there are no duplicates

# load the big data data - the csv files of coverage; let the user know how long this is taking
if (debug) print(paste("Starting to load the data at",format(Sys.time(), "%a %b %d %H:%M:%S %Y")))

tumor.data <- load.exome.data(tumor.lane.data)
normal.data <- load.exome.data(normal.lane.data)

if (debug) print(paste("Intersecting the normal and tumor bait lists, normal data has",nrow(normal.data),"rows, tumor data has",nrow(tumor.data),"rows"))

# cut the data to the intersect of the data
target.intersect <- intersect(rownames(normal.data),rownames(tumor.data))
print(paste("intersection of the tumor and normal has",length(target.intersect),"rows"))

normal.data <- intersect.tumor.normal.targets(normal.data,target.intersect)
tumor.data <- intersect.tumor.normal.targets(tumor.data,target.intersect)

# load up our bait factor
print("Subsetting the data")
baits.subset <- baits[is.element(baits$name,rownames(normal.data)),]
bait.factor <- read.delim(bait.factor.file)
bait.factor[bait.factor[,1]<=0,2] = epsilon
bait.factor = bait.factor[is.element(rownames(bait.factor),rownames(tumor.data)),]

# mean center the tumor and normal samples
tumor.data = sweep(tumor.data,2,apply(tumor.data,2,mean),"/")

# get the sex assignments; a list of each sex with their associated columns
sex.columns = process.sex.assignments(sex.calls,colnames(tumor.data))

# order both our tumors and normals; females then males
tumor.data <- cbind(tumor.data[!sex.columns],tumor.data[sex.columns])
normal.data <- cbind(normal.data[!sex.columns],normal.data[sex.columns])

# order the sex data
sex.columns <- sex.columns[order(sex.columns)]

# get a list of the split sex - autosomal data
tumor.split = split.out.sex.chromosomes(tumor.data,sex.chromosomes,baits)
normal.split = split.out.sex.chromosomes(normal.data,sex.chromosomes,baits)

# now that we're processing in pieces, we have to save this in a list format
processed.data = list()

# now normalize coverage across each of the split data (split by sex/autosomal data)
for (n in c("autosome",sex.chromosomes)) {
    print(paste("Starting analysis for chromosome",n))
    log.tumor = mean.center.log.transform(tumor.split[[n]])
    log.normal = mean.center.log.transform(normal.split[[n]])

    # do the initial block normalization
    bgs.center = rep(mean(colMeans(log.normal)),nrow(log.normal)) # rep(0.0,nrow(normal.data)) # log.normals[,ncol(log.normals)] # rep(0.0,nrow(normal.data))
    names(bgs.center) <- rownames(log.normal)
    log.normal = data.frame(log.normal[,1:ncol(log.normal)] - bgs.center,check.names=F)
    log.tumor = data.frame(log.tumor - bgs.center,check.names=F)

    print(paste("Data loaded and means and log values calculated...",n))
    #load up the whole data set into the tangent normalization process, and calibrate each tumor against the matrix
    if (use.histo.data) {
        print("loading the historical data...")
        log.normal = additional.normals(tangent.database.location,log.normal,build.version,analysis.set.name)
    }

    target.intersect <- intersect(rownames(log.normal),rownames(log.tumor))
    log.normal <- intersect.tumor.normal.targets(log.normal,target.intersect)
    log.tumor <- intersect.tumor.normal.targets(log.tumor,target.intersect)

    print("About to perform the SVD (PI)")

    # if we're a sex chromosome, split the two sample piles by their sex assignment
    if (n != "autosome") {
        calibrated.tumor = data.frame()
        # split out the two piles, male and female

        if (sum(!sex.columns) > 0) {
            print("FEMALE")
            females.ln = log.normal[,!sex.columns]
            females.lt = log.tumor[,!sex.columns]
            pseudo.inverse.norm  <- pseudoinverse(data.matrix(females.ln))
            calibrated.tumor <- calibrate.tumors(data.matrix(females.lt),data.matrix(pseudo.inverse.norm),data.matrix(females.ln),bgs.center,first=TRUE)
            colnames(calibrated.tumor) <- colnames(females.lt)
        }

        if (sum(sex.columns) > 0) {
            print("HERE")
            males.ln = log.normal[,sex.columns]
            males.lt = log.tumor[,sex.columns]
            pseudo.inverse.norm  <- pseudoinverse(data.matrix(males.lt))
            calibrated.tumor <- cbind(calibrated.tumor,calibrate.tumors(data.matrix(males.lt),data.matrix(pseudo.inverse.norm),data.matrix(males.ln),bgs.center,first=TRUE))
            print(dim(calibrated.tumor))
            print(length(c(colnames(calibrated.tumor)[!sex.columns],colnames(males.lt))))
            colnames(calibrated.tumor) <- c(colnames(calibrated.tumor)[!sex.columns],colnames(males.lt))
        }
    }
    else {
        pseudo.inverse.norm  <- pseudoinverse(data.matrix(log.normal))
        calibrated.tumor <- calibrate.tumors(data.matrix(log.tumor),data.matrix(pseudo.inverse.norm),data.matrix(log.normal),bgs.center,first=TRUE)
        colnames(calibrated.tumor) <- colnames(tumor.data)
    }
    # save off each piece of our work as we go
    tn = calibrated.tumor# log.tumor) # ,log.normal)
    processed.data[[n]] = tn
}

# save off the data to a Rdata object for later -- not anymore
# save.off.processed.data(log.normals,log.tumors,calibrated.tumors,baits,paste(tangent.database.output,build.version,sep="/"),analysis.set.name,build.version)
calibrated.tumors <- rbind(data.frame(processed.data[["autosome"]],check.names=F),data.frame(processed.data[["X"]],check.names=F),data.frame(processed.data[["Y"]],check.names=F))
rownames(calibrated.tumors) <- rownames(tumor.data)

# output the raw data and plots for each sample
output.and.plot.data(calibrated.tumors,tumor.data,baits,output.location,signal.files)

# output metrics for each sample
create.sample.metric.files(calibrated.tumors,mean.center.log.transform(tumor.data),output.location)
