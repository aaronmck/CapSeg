DefaultGermlineHetFileParser = function(germline.het.fn) 
{
   dat = read.delim(germline.het.fn, stringsAsFactors=FALSE, check.names=FALSE, blank.lines.skip=TRUE, comment.char="#")
		
   if (!all((dat$Start_position == dat$End_position))) stop (paste("There's something wrong with the Germline Het Table."))

   dat$Chromosome = as.character(gsub("Y", "24", gsub("X", "23", dat$Chromosome)))

   return(dat)
}

DefaultCapsegFileParser = function(capseg.probe.fn, drop.y=TRUE) {
	
	capseg.d = read.delim(capseg.probe.fn, strings=F)
	rownames(capseg.d) <- capseg.d[,1]
	capseg.d[,1] <- NULL
	
	capseg.d[,5] <- capseg.d[,2] + floor((capseg.d[,3] - capseg.d[,2])/2)
	capseg.d = capseg.d[,c(1, 5, 2, 3, 4)] 
	names(capseg.d) <- c('Chromosome', "Center.bp", "Start.bp", "End.bp", "Intensity")
	capseg.d[["Intensity"]] <- 2^(capseg.d[["Intensity"]] + 1)
	capseg.d[["Chromosome"]] = gsub("Y", "24", gsub("X", "23", capseg.d[["Chromosome"]]))
		
	if (drop.y) {
		capseg.d = capseg.d[capseg.d$Chromosome != "24", ]
	}
	if (drop.x) {
		capseg.d = capseg.d[capseg.d$Chromosome != "23", ]
	}
	return(capseg.d)
}

DefaultSnpFileParser <- function(snp.fn, verbose=FALSE) {
  snp.try <- try(load(snp.fn), silent=TRUE)
  ## FIXME: Some of Scott's old data uses RData inputs instead of
  ## text snp files. These are already diced per-sample and don't need
  ## to be processed separately
  if ((!inherits(snp.try, "try-error")) && (exists("dat")) && (identical(colnames(dat), c("A", "B")))) {
    if (verbose) {
      print("Using RData format SNP file")
    }
    allele.data <- dat
  } else {
    if (verbose) {
      print("Using text format SNP file")
    }
    allele.data <- ReadDicedSnpPipelineFile(snp.fn, sample.name)
  }
  
  return(allele.data)
}


DefaultCnFileParser <- function(cn.fn, verbose=FALSE) {
	if (verbose) {
		print("Using text format CN file")
	}
	cn.data <- ReadDicedCnPipelineFile(cn.fn)
	
	return(cn.data)
}


FullSnpFileParser <- function(snp.fn, sample.name, verbose=FALSE) {
  allele.data <- SnpToBin(snp.fn, sample.name)
  
  return(allele.data)
}

DefaultClustersFileParser <- function(clusters.fn, verbose=FALSE) {
  ## If the load() works, it pulls in 'clusters'
  try.val <- suppressWarnings(try(load(clusters.fn), silent=TRUE))
  if (inherits(try.val, "try-error")) {
    clusters <- read.delim(clusters.fn, row.names=1)
    if (verbose) {
      print("Loaded the clusters file as text")
    }
  } else {
    if (verbose) {
      print("Loaded the clusters file as RData")
    }
  }
  
  return(clusters)
}

BirdseedClustersFileParser <- function(clusters.fn, verbose=FALSE) {
  tmp.fn <- tempfile()
  tr.cmd <- paste("cat ", clusters.fn, " | tr ';' \\\t | tr ' ' \\\t > ",
                  tmp.fn, sep="")
  system(tr.cmd)
  clusters <- read.delim(tmp.fn, header=0, row.names=1, 
                    colClasses=c("character", rep("numeric", 18)))
  clusters <- clusters[, c(2, 8, 14, 15, 9, 3) - 1]
  colnames(clusters) <- c("AA.a", "AB.a", "BB.a", "BB.b", "AB.b", "AA.b")
  rownames(clusters) <- sub("-2$", "", rownames(clusters))
  
  return(clusters)
}

