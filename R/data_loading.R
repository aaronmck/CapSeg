
#' load the exome data and process it down to a data matrix
#' @param data the data csv file name
#' @param baits the baits name
#' @param removeBadBaitsAndLanes remove the bad baits and lanes from the matrix
#' @param minimumBaitOverlap we require at least this ratio of overlap between the bait names and the pulldown capture names
#' @return a data matrix of the coverage
#' @keywords exome coverage
#'
load.exome.data = function(
		data, # the data matrix produced by the rest of capseg
		baits
) {
	return(read.delim(data,check.names=F,stringsAsFactors=F))
}

#' load the exome data and process it down to a data matrix
#' @param data the data csv file name
#' @param baits the baits name
#' @param removeBadBaitsAndLanes remove the bad baits and lanes from the matrix
#' @param minimumBaitOverlap we require at least this ratio of overlap between the bait names and the pulldown capture names
#' @return a data matrix of the coverage
#' @keywords exome coverage
#'
load.float.exome.data = function(
		data, # the data matrix produced by the rest of capseg
		baits) {

	# this next part is a bit hacky, but its a speedup and a huge memory win over the default read.csv
	line.count = as.numeric(system(paste("wc -l",data,"| cut -f 1 -d ' '"),intern = TRUE))
	col.count = as.numeric(system(paste("head -1",data,"| tr ',' '\n' | wc -l"),intern = TRUE))
	col.types = c("character",rep("numeric", times=col.count-1))

	# load up the raw data - there shouldn't be duplicate row names anymore
	raw <- read.table(data,sep=",",row.names=1,header=TRUE,stringsAsFactors=F,check.names=F,nrows=line.count,comment.char = "",colClasses=col.types)

	# raw2 <- join(baits,raw,by = (baits$name,rownames(raw)),type="left")
	return(raw)
}

# for the tumor samples, return a
#' @param data the data csv file name
#' @param baits the baits name
#' @param removeBadBaitsAndLanes remove the bad baits and lanes from the matrix
#' @param minimumBaitOverlap we require at least this ratio of overlap between the bait names and the pulldown capture names
#' @return a data matrix of the coverage
#' @keywords exome coverage
#'
process.sex.assignments <- function(sex.calls,tumor.samples) {
    # open the sex calls file
    s.calls <- read.delim(sex.calls)
    is.male = s.calls$call == "MALE"
    return(is.male[order(match(s.calls$sample,tumor.samples))])
}

