## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

ReadGladMat <- function(glad.fn, sample.name=NULL, drop.x=TRUE, drop.y=TRUE, glad.log=FALSE, type = "snp", verbose=FALSE) {

					
  sp.res <- strsplit(glad.fn, ".", fixed=TRUE)[[1]]
  ln <- length(sp.res)

  if (verbose) {
    print(paste("loading glad mat", glad.fn))
  }
  glad.mat <- read.delim(glad.fn, row.names=NULL, as.is=TRUE)
  glad.mat[, 2] <- gsub("X", "23", glad.mat[, 2])
  glad.mat[, 2] <- gsub("Y", "24", glad.mat[, 2])
  glad.mat[, 2] <- as.numeric(glad.mat[, 2])
  
  if (verbose) {
    print("summary pre-processing")
  }

  glad.mat[, 2] <- sub("chr", "", glad.mat[, 2])
  
  ## enforce standard (GCM) names   
  colnames(glad.mat) <- c( "Sample", "Chromosome", "Start.bp", "End.bp", "Num.Probes",
                          "Seg.CN" )

  ## version with one merged GLAD file
  if (!is.null(sample.name)) {
    if (verbose) {
      print(paste("Extracting out sample", sample.name,
                  "from the seg file"))
    }
      
    idx <- glad.mat[, 1] == sample.name 
    if (sum(idx) == 0) {
      ## this means that none of the rows contained our sample
      stop("Unable to find any data in the sample column matching the sample name, ",
           "even though it was specified.  Check your segmentation file")
    }
    glad.mat <- glad.mat[idx, ] 
    if (nrow(glad.mat) == 0) {
      stop()
    }

    if (verbose) {
      print("read GLAD mat: ")
      print(dim(glad.mat))
    }
  }
  
  dropChroms <- character()
  if (drop.x) {
    dropChroms <- c(dropChroms, c("23", "X"))
  }
  if (drop.y) {
    dropChroms <- c(dropChroms, c("24", "Y"))
  }
  if (length(dropChroms) > 0) {
    glad.mat <- DropChromosomes(glad.mat, dropChroms)
  }
    
  seg.lens <- (glad.mat[, "End.bp"] - glad.mat[, "Start.bp"]) / 1e7
  glad.mat <- cbind(glad.mat, seg.lens)
  colnames(glad.mat)[ncol(glad.mat)] <- "length"

  data <- list(glad=glad.mat)
  seg.info <- list()
  NSEG <- nrow(glad.mat)
  cnames <- c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "copy_num")

  seg.info <- matrix(NA, nrow=NSEG, ncol=length(cnames))
  colnames(seg.info) <- cnames
  cols <- c("Chromosome", "Start.bp", "End.bp") 
  seg.info[, cols[1]] <- glad.mat[, cols[1]]
  seg.info[, cols[2]] <- glad.mat[, cols[2]]
  seg.info[, cols[3]] <- glad.mat[, cols[3]]
  
  if(glad.log == TRUE & type == "snp") {
    seg.info[, "copy_num"] <- 2^data[["glad"]][, "Seg.CN"]
  } else if (glad.log == TRUE & type == "capseg") {
	  seg.info[, "copy_num"] <- 2^(data[["glad"]][, "Seg.CN"] + 1)
  } else {
    seg.info[, "copy_num"] <- data[["glad"]][, "Seg.CN"]
  }
   
  seg.info[, "length"] <- data[["glad"]][, "length"]
  seg.info[, "n_probes"] <- data[["glad"]][, "Num.Probes"]
  
  seg.info <- matrix(as.numeric(seg.info), nrow=nrow(seg.info), ncol=ncol(seg.info))
  colnames(seg.info) <- cnames
  if (verbose) {
    print(seg.info[1:10, 1])
  }
  return(list(seg.info=seg.info))
}

DropChromosomes <- function(glad.mat, chromNames) {
  idx <- which(glad.mat[, "Chromosome"] %in% chromNames)
  if (length(idx) > 0) {
    glad.mat <- glad.mat[-idx, ]
  }
  return(glad.mat)
}
