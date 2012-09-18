# run the individual segmentation
# ---------------------------------------------------------------------------
# Aaron
#
# Jan 3rd, 2011
#
# This script processes raw data from the exome coverage tool and produces
# segmentation calls and other data needed for the HapSeg Algorithm for
# each of the samples seen in the data
#
# ---------------------------------------------------------------------------


# get the command line options, and make sure we got the right length
library(optparse)
# require("genomic_plot.R")
option_list <- list(
                    # we include the R dat files from
                    make_option(c("--bam.file.name"),help="our bam file name",default="./"),
                    make_option(c("--bam.to.sample.file"),help="bam to sample file",default="error"),
                    make_option(c("--signal.file"),help="the signal directory",default="error"),
                    make_option(c("--r.dir"),help="the r source directory",default="error"),
                    make_option(c("--output.filename"),help="the output directory",default="error")
                    )
opt <- parse_args(OptionParser(option_list=option_list))

source(paste(opt$r.dir,"/R/genomic_plot.R",sep="")) #,verbose=T)


# load the arguments
bam.name = opt$bam.file.name
bam.to.sample.file = opt$bam.to.sample.file
output.location = opt$output.location
output.filename = opt$output.filename
signal.file = opt$signal.file

# use vega for the segmentation - it has worked better than CBS for our applications
library("DNAcopy")

                                        # load the arguments
print(paste("bam.name=",bam.name))
print(paste("output.location=",output.filename))
print(paste("bam.to.sample.file=",bam.to.sample.file))
print(paste("signal.file",signal.file))
#save.image(".parameters.seg.rData")

# check that we got the signal file for our bam file; otherwise we don't run
bam.to.sample.table <- read.table(bam.to.sample.file,sep="\t",header=T,check.names=T)
sample.name = sub(" ",".",sub("-",".",bam.to.sample.table[bam.to.sample.table[,1]==bam.name,2]))
print(paste("Sample name",sample.name))

# if we dont have a sample if the bam to sample table, get out!
if (length(sample.name) < 1) {
	print(paste("unable to find a sample for the bam file",bam.name))
	q(status=1)
}

# load up the signal data
print(paste("going to load signal file ",signal.file))

if (!file.exists(signal.file)) {
    print(paste("the signal file doesn't exist",signal.file))
}

# load up the tumor data
signal.data <- read.table(signal.file,sep="\t",header=T,check.names=F)

# ------------------------------------------------ writeSegmentationTable -------------------------------------------------
# try to get segmentations for all of the chromosomes in a samples - return the result table we expect for HAPSEG
createSegmentationForAllContig <- function(tumor.data,sample.name) {
                sg <- tumor.data[,2:5]
                sg[,1] <- sub("chr","",sg[,1])
                sg <- sg[complete.cases(sg),]
                seg <- data.frame(vega(sg, unique(sg[,1]),min_region_size=100,beta=0.5))
		return(cbind(sample.name,seg[,1:5]))
}

createSegmentationForAllContig2 <- function(tumor.data,sample.name) {
               temp.data <- tumor.data[is.finite(tumor.data[,5]),]

               # run the segmentation algorithm - CBS
               cna_data <- CNA(genomdat=as.matrix(temp.data[,5]), chrom=temp.data$contig, maploc=as.matrix(temp.data$start+(temp.data$stop-temp.data$start)), data.type=c("logratio"), sampleid=colnames(tumor.data)[5])
               smoothed.CNA.object <- smooth.CNA(cna_data)
               results <- segment(smoothed.CNA.object) # ,alpha=0.1,nperm=5000)
               return (results$output)
          }
print("Data has been loaded...")

# get the sample name
sample.name <- colnames(signal.data)[5]

# announce the segmentation
print(paste("Running the segmentation for sample...",sample.name))

# the base name - create the directory if it doesn't exist
segmentation.dir = unlist(strsplit(output.filename,"/"))
segmentation.dir = paste(segmentation.dir[1:length(segmentation.dir)-1],collapse="/")
print(paste("Making segmentation directory",segmentation.dir))
dir.create(segmentation.dir,recursive=T)

# segment the data
file.name <- output.filename

# make the segmentation data using cbs
results <- createSegmentationForAllContig2(signal.data,sample.name)

# median center the segments
results <- cbind(results[,1:5],results[,6]-median(results[,6]))

# write the segments table
file.name = sub(" ","_",file.name)
print(paste("Writing to segmentation file name ",file.name))
write.table(results,file=file.name,quote=F,row.names=F,col.names=c("Sample","Chromosome","Start","End","Num_Probes","Segment_Mean"),sep="\t")

