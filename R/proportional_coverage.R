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
options(stringsAsFactors=F)

option_list <- list(
                    # we include the R dat files from
                    make_option(c("--intervals"),help="the input file of targets we pulled coverage for (bed file)",default="./"),
                    make_option(c("--output.file"),help="the output file of proportional coverage",default="./"),
                    make_option(c("--output.cr.stat.file"),help="the output file of the cr.stats for each lane",default="./"),
                    make_option(c("--output.column.sums"),help="the column sums output",default="./"),
                    make_option(c("--input.file"),help="the input file of coverage",default="./")
                    )

# add the options
opt <- parse_args(OptionParser(option_list=option_list))
#save.image(".parameters.pc.Rdata")

# load the coverage
print(paste("Loading file",opt$input.file,"from disk"))
input <- read.delim(opt$input.file,stringsAsFactors=F,check.names=F)
input.colnames <- colnames(input)
intervals <- read.delim(opt$intervals,header=F,stringsAsFactors=F)

# intervals <- read.delim(interv,header=F)
colnames(intervals) <- c("contig","start","stop","names")

# this is very contrived, but it saves a ton of time over merge
input <- merge(intervals,input,by.x="names",by.y="row.names",all.x = T,sort=F)

input <- data.frame(input[,5:ncol(input)])
colnames(input) <- input.colnames
rownames(input) <- intervals$names
input[is.na(input)] = 0

# column sums
clSums = data.frame(colSums(input))
clSums[clSums==0] = 1 # if any of the columns have a sum of zero, set it to one

# output the column sum data
colnames(clSums) <- c("sum")
write.table(clSums,file=opt$output.column.sums,sep="\t",quote=F)

write.table(input,file=paste(opt$output.file,"midway",sep="."),sep="\t",quote=F)
# calculate the copy ratio stats for each lane
input.stat = sweep(input,2,colSums(input),"/")

cr.stat <- data.frame(apply(input,2,function(x) {median(abs(diff(x)))}))
colnames(cr.stat) <- c("cr.stat")

# output the cr stats
write.table(cr.stat,file=opt$output.cr.stat.file,sep="\t",quote=F,row.names=T)

# output the proportional table
write.table(input,file=opt$output.file,sep="\t",quote=F)
