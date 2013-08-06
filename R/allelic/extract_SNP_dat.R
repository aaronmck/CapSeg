## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.




ExtractCaptureDat <- function(capseg.probe.fn=NULL, seg.dat, germline.het.fn, drop.x, drop.y, verbose=FALSE, germline.het.file.parser=DefaultGermlineHetFileParser) {
  
  capseg.d = DefaultCapsegFileParser(capseg.probe.fn, drop.x=drop.x, drop.y=drop.y)
  germline.hets = germline.het.file.parser(germline.het.fn)
  
  h.seg.dat <- GetCaptureAsSegs(seg.dat, capseg.d, germline.hets, verbose=verbose)

  return(list(seg.dat=seg.dat, as.res=list(h.seg.dat=h.seg.dat)))
}


ExtractArrayDat <- function(array.name, genome.build, use.pop, use.normal,
  normal, impute.gt, adj.atten, platform, 
  seg.dat, snp.fn, cn.fn,
  calls.fn, mn.sample,
  platform.vals.fn,
  calibrate.data=FALSE, clusters.fn=NULL, verbose=FALSE,
  snp.file.parser=DefaultSnpFileParser,
  cn.file.parser=DefaultCnFileParser,
  clusters.file.parser=DefaultClusterFileParser) {

  # data(list=GetPlatformDataName(platform), package="HAPSEG")
  # if (! genome.build %in% names(platform.annots)) {
  #   stop("Unsupported genome build, ", genome.build,
  #     ", for this platform: ", platform)
  # }
  # platform.vals <- mget(c("snp.freqs", "dbSNP.annot", "snp.annot", "post.birdseed.calibration"), env=platform.annots[[genome.build]])
  if (verbose) print(paste("Using plaform values from:", platform.vals.fn, "with platform set to:", platform))
  platform.vals <- readRDS(platform.vals.fn)

  snp.freqs <- platform.vals[["snp.freqs"]]
  dbSNP.annot <- platform.vals[["dbSNP.annot"]]
  snp.annot <- platform.vals[["snp.annot"]]
  post.birdseed.calibration <- platform.vals[["post.birdseed.calibration"]]


  snp.data <- snp.file.parser(snp.fn, sample.name=array.name, verbose=verbose)

  ## the allele frequency data 
  if (calibrate.data) {
    if (is.null(clusters.fn)) {
      stop("Calibrating data, but there is no clusters file!")
    }
    if (verbose) {
      print("Calibrating data")
    }
    snp.d <- CalibrateAsDat(snp.data, clusters.fn, 
      clusters.file.parser=clusters.file.parser,
      verbose=verbose)
   } else {
    if (!verbose) {
     print("Not calibrating data")
   }
   snp.d <- snp.data
 }

  #### Comment out to not do cross-hybridization correction
  snp.d <- PostCalibrateAsDat(snp.d, post.birdseed.calibration[["snp.tx"]],
   adj.atten=adj.atten, verbose=verbose)
  ####
  
  nas <- apply(is.na(snp.d), 1, sum)
  snp.d <- snp.d[nas == 0, ]
  
  # No calibration is needed for cn probes because they have already been calibrated in SnpPipeline
  # Just remove NAs
  
  cn.d <- cn.file.parser(cn.fn)
  nas <- apply(is.na(cn.d), 1, sum)
  cn.d <- cn.d[nas == 0, ,drop=F] 
  
  as.res <- GetArrayAlleleSegData(snp.d, cn.d, snp.annot, glad.mat=seg.dat,
   snp.freqs, use.pop, use.normal,
   normal, NA, dbSNP.annot, impute.gt,
   calls.fn, mn.sample, verbose=verbose)
  
  return(list(seg.dat=seg.dat, as.res=as.res))
}


ExtractProbeDat <- function(array.name, capseg.sample.name, genome.build, use.pop, use.normal,
  normal, impute.gt, adj.atten, 
  seg.dat, snp.fn=NULL, cn.fn,  capseg.probe.fn=NULL, germline.het.fn, tumor.sample.barcode,
  calls.fn, mn.sample, drop.x, drop.y,
  platform.vals.fn,
  calibrate.data=FALSE, clusters.fn=NULL, verbose=FALSE,
  snp.file.parser=DefaultSnpFileParser,
  cn.file.parser=DefaultCnFileParser,
  capseg.file.parser=DefaultCapsegFileParser,
  germline.het.file.parser=DefaultGermlineHetFileParser,
  clusters.file.parser=DefaultClusterFileParser) {
	## Returns:
	## seg.dat - segmentation data file
	## as.res - allele specific segmentation data from all probe types

  platform <- "SNP_6.0"
  # data(list=GetPlatformDataName(platform), package="HAPSEG")
  # if (! genome.build %in% names(platform.annots)) stop("Unsupported genome build, ", genome.build, ", for this platform: ", platform)
  # platform.vals <- mget(c("snp.freqs", "dbSNP.annot", "snp.annot", "post.birdseed.calibration"), env=platform.annots[[genome.build]])
  if (verbose) print(paste("Using plaform values from:", platform.vals.fn, "with platform set to:", platform))
  platform.vals <- readRDS(platform.vals.fn)

  snp.freqs <- platform.vals[["snp.freqs"]]
  dbSNP.annot <- platform.vals[["dbSNP.annot"]]
  snp.annot <- platform.vals[["snp.annot"]]
  post.birdseed.calibration <- platform.vals[["post.birdseed.calibration"]]
  

  # seg.dat = CreateMergedSegFile( array.seg.fn , capseg.seg.fn)

  if (!is.null(snp.fn)){
    snp.data <- snp.file.parser(snp.fn, sample.name=array.name, verbose=verbose)

  	## the allele frequency data 
  	if (calibrate.data) {
  		if (is.null(clusters.fn)) {
  			stop("Calibrating data, but there is no clusters file!")
  		}
  		if (verbose) {
  			print("Calibrating data")
  		}
  		snp.d <- CalibrateAsDat(snp.data, clusters.fn, 
        clusters.file.parser=clusters.file.parser,
        verbose=verbose)
     } else {
      if (!verbose) {
       print("Not calibrating data")
     }
     snp.d <- snp.data
   }

  	#### Comment out to not do cross-hybridization correction
  	snp.d <- PostCalibrateAsDat(snp.d, post.birdseed.calibration[["snp.tx"]],
     adj.atten=adj.atten, verbose=verbose)
  	####
  	
  	nas <- apply(is.na(snp.d), 1, sum)
  	snp.d <- snp.d[nas == 0, ]
  	
  	# No calibration is needed for cn probes because they have already been calibrated in SnpPipeline
  	# Just remove NAs
  	
  	cn.d <- cn.file.parser(cn.fn)
  	nas <- apply(is.na(cn.d), 1, sum)
  	cn.d <- cn.d[nas == 0, ,drop=F] 
     
  }
	
  if (!is.null(capseg.probe.fn)) {
  	capseg.d <- capseg.file.parser(capseg.probe.fn, drop.x, drop.y)
  	germline.hets <- germline.het.file.parser(germline.het.fn)
      
  }

	as.res <- GetAlleleSegData(snp.d, cn.d, capseg.d, germline.hets, snp.annot, glad.mat=seg.dat,
   snp.freqs, use.pop, use.normal,
   normal, NA, dbSNP.annot, impute.gt,
   calls.fn, mn.sample, verbose=verbose)
	
	return(list(seg.dat=seg.dat, as.res=as.res))
}


CreateMergedSegFile <- function(array.seg.fn, array.name, capseg.seg.fn, capseg.sample.name, drop.x, drop.y, verbose=FALSE) {
  # Output seg file is a seg file from the union of breakpoints of the input seg files.

  if (!is.null(array.seg.fn)) {
    array.seg.dat <- ReadGladMat(array.seg.fn, array.name, glad.log=TRUE, drop.x=drop.x, drop.y=drop.y, type="snp", verbose=verbose)[[1]]
  }
  if (!is.null(capseg.seg.fn)) {
    capseg.seg.dat <- ReadGladMat(capseg.seg.fn, capseg.sample.name, glad.log=TRUE, drop.x=drop.x, drop.y=drop.y, type="capseg", verbose=verbose)[[1]]
  } 


  
  if (!is.null(array.seg.fn) & !is.null(capseg.seg.fn)) {
    seg.dat = MergeSegs(array.seg.dat, capseg.seg.dat)
  } else if (is.null(array.seg.fn)) {
    seg.dat = capseg.seg.dat
  } else if (is.null(capseg.seg.fn)) {
    seg.dat = array.seg.dat
  } else {
    stop("Both capseg and snp seg files are NULL.")
  }

  return(seg.dat[,c("Chromosome", "Start.bp", "End.bp")])
}



ExtractSnpDat <- function(sample.name, genome.build, use.pop, use.normal,
  normal, impute.gt, adj.atten, platform, seg.fn,
  snp.fn, drop.x, drop.y, calls.fn, mn.sample,
  calibrate.data=FALSE, clusters.fn=NULL, verbose=FALSE,
  snp.file.parser=DefaultSnpFileParser, 
  clusters.file.parser=DefaultClusterFileParser) {
  ## Returns:
  ## seg.dat - segmentation data file
  ## as.res - allele specific segmentation data

  data(list=GetPlatformDataName(platform), package="HAPSEG")
  if (! genome.build %in% names(platform.annots)) {
    stop("Unsupported genome build, ", genome.build,
     ", for this platform: ", platform)
  }
  platform.vals <- mget(c("snp.freqs", "dbSNP.annot", "snp.annot",
    "post.birdseed.calibration"),
  env=platform.annots[[genome.build]])
  snp.freqs <- platform.vals[["snp.freqs"]]
  dbSNP.annot <- platform.vals[["dbSNP.annot"]]
  snp.annot <- platform.vals[["snp.annot"]]
  post.birdseed.calibration <- platform.vals[["post.birdseed.calibration"]]
  
  ## read the segmentation data in
  seg.dat <- ReadGladMat(seg.fn, sample.name, glad.log=TRUE,
   drop.x=drop.x, drop.y=drop.y, verbose=verbose)


  allele.data <- snp.file.parser(snp.fn, verbose=verbose)
  
  ## the allele frequency data 
  if (calibrate.data) {
    if (is.null(clusters.fn)) {
      stop("Calibrating data, but there is no clusters file!")
    }
    if (verbose) {
      print("Calibrating data")
    }
    as.d <- CalibrateAsDat(allele.data, clusters.fn, 
     clusters.file.parser=clusters.file.parser,
     verbose=verbose)
    } else {
      if (!verbose) {
        print("Not calibrating data")
      }
      as.d <- allele.data
    }


  #### Commented out to not do cross-hybridization correction
#  as.d <- PostCalibrateAsDat(as.d, post.birdseed.calibration[["snp.tx"]],
#                             adj.atten=adj.atten, verbose=verbose)
	####
	
  nas <- apply(is.na(as.d), 1, sum)
  as.d <- as.d[nas == 0, ]

  as.res <- GetAlleleSegData(as.d, snp.annot, seg.dat[["seg.info"]],
   snp.freqs, use.pop, use.normal,
   normal, NA, dbSNP.annot, impute.gt,
   calls.fn, mn.sample, verbose=verbose)

  return(list(seg.dat=seg.dat, as.res=as.res))
}

GetPlatformDataName <- function(platform) {
  platform <- switch(platform,
   "SNP_6.0"="6.0",
   "SNP_250K_STY"="sty",
   stop("Unsupported platform: ", platform))
  return(paste("platform", platform, sep="_"))
}

ReadDicedSnpPipelineFile <- function(fn, sample.name) {
  snp.file <- read.delim(fn, as.is=TRUE)
  ## There's an extra header line in here
  snp.file <- snp.file[-1, ]
  rownames(snp.file) <- snp.file[, 1]
  snp.file <- snp.file[, -1]
  colnames(snp.file) <- c("A", "B")
  ## Filter out the non-SNP markers
  snp.file <- snp.file[grep("^SNP", rownames(snp.file)), , drop=FALSE]
  for (i in colnames(snp.file)) {
    snp.file[, i] <- as.numeric(snp.file[, i])
  }
  return(snp.file)
}

ReadDicedCnPipelineFile <- function(cn.fn) {
	# Bryan Hernandez
	
	cn <- read.delim(cn.fn, as.is=TRUE)
	## There's an extra header line in here
	cn <- cn[-1, ]
	rownames(cn) <- cn[, 1]
	cn <- cn[, -1, drop=F]
	colnames(cn) <- c("Intensity")
	## Filter out the non-CN markers
	cn <- cn[grep("^CN", rownames(cn)), , drop=FALSE]
	cn[,1] <- as.numeric(cn[,1])
	return(cn)
}


## FIXME: Supercedes the FIXME below. This needs to be replaced with
## diced input. This is temporary
ReadSnpPipelineFile <- function(fn, sample.name) {
  snp.file <- read.delim(fn, as.is=TRUE)
  cols <- paste(sample.name, c("alleleA", "alleleB"), sep="_")
  return(snp.file[, cols])
}

## FIXME: replace split.snp.file with this (Scott's older version).
## Both will be replaced by diced input
SnpToBin <- function(fn, col.name) {
  snp <- ReadCol(fn, col.name, save.rownames=TRUE)
  ## Filter out the non-SNP markers
  snp <- snp[grep("^SNP", rownames(snp)), , drop=FALSE]
  a.idx <- grep("-A", rownames(snp))
  b.idx <- grep("-B", rownames(snp))

  a.snps <- as.numeric(snp[a.idx, ])
  b.snps <- as.numeric(snp[b.idx, ])

  names(a.snps) <- gsub("-A", "", rownames(snp)[a.idx])
  names(b.snps) <- gsub("-B", "", rownames(snp)[b.idx])

  if (!setequal(rownames(a.snps), rownames(b.snps))) {
    stop("-A and -B snps do not match")
  }

  df <- cbind(a.snps, b.snps)
  colnames(df) <- c("A", "B")

  return(df)
}

## FIXME: can this be made better?
split.snp.file <- function(SNP.FN) {
  tbl <- read.table(SNP.FN,skip=2,header=F)
  
  colnames(tbl) <- c("SNP","signal")
  tbl <- tbl[substr(tbl$SNP,1,3)=="SNP",]

  ## check that we have even piles
  if (nrow(tbl) %% 2 != 0) {
    print(paste("We have an uneven number of rows in the SNP file (of rows starting with SNP): ",nrow(tbl)))
    q(status=1)
  }

  a.col = tbl[seq(1,nrow(tbl),2),]
  b.col = tbl[seq(2,nrow(tbl),2),]

  ## make sure the sizes are the same at least
  if (nrow(a.col) != nrow(b.col)) {
    print("uneven number of rows in the data: ",nrow(a.col),nrow(b.col))
    q(status=1)
  }
  
  ## should we really dig in and check our results?
  check.results = FALSE # turning this on is really expensive
  if (check.results) {
    count <- 0
    for (i in seq(1,nrow(a.col))) {
      val.a = substr(a.col[i,1],1,length(a.col[i,1])-2)
      val.b = substr(b.col[i,1],1,length(b.col[i,1])-2)
      if (val.a != val.b) {
        print(paste("Unable to match column",i,"val a =",val.a,", val b =",val.b))
        q(status=1)
      }
      count <- count + 1
      if(count %% 1000 == 0) {print(paste("Checked",count,"records"))}
    }
  }

  ## make the final return data frame - the lapply to get the names is really slow....
  ##data.row.names <- unlist(lapply(a.col[,1],function(x){return(substr(x,1,nchar(x)-2))}))
  qq <- strsplit(as.character(a.col[,1]), "-")
  data.row.names <- sapply(qq, function(x){paste(x[1],x[2],sep="-")})
  ret <- cbind(a.col[,2],b.col[,2])
  rownames(ret) <- c(data.row.names)
  colnames(ret) <- c("A","B")
  return(ret)
}

