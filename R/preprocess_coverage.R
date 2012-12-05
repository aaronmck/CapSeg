library(bigmemory)

# load each of the coverage matricies in order, getting the median of each of the rows
coverage.matricies = List()

# open up the list of coverage files
normal.coverage.files = read.delim("/Users/aaron/Desktop/CapSeg/PR_Mouse_Jacks_SCLC_Capture/results/normalCoverageFiles.txt")
tumor.coverage.files = read.delim("/Users/aaron/Desktop/CapSeg/PR_Mouse_Jacks_SCLC_Capture/results/tumorCoverageFiles.txt")

load.coverage.data <- function(data.row) {
    print(paste("loading",data.row[2]))
    coverage <- read.big.matrix(data.row[2],
                                sep="\t",
                                header=TRUE,
                                type="double",
                                row.names=1,
                                backing.file=paste(data.row[2])
}

# now cycle through each, building up the master matrix
coverage.matricies = apply(normal.coverage.files,1,load.coverage.data)

# now get the total coverage, and divide by it

# now get the median of each lane, and divide by that

# now get the CR stat: the med(abs(dif(x))), and drop the worst lanes

# now collapse to the sample level matrix

# now output the data

