## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CallBeagle <- function(chr, allele.tab, dbSNP.annot, platform,
                       tmp.dir, plate.name, phased.bgl.dir, verbose=FALSE) {  
  kNImp <- 0
  kNIter <- 10
  
  obs.gt.fn <- file.path(tmp.dir, paste(plate.name, "GT", sep="."))
  markers.fn <- file.path(tmp.dir, paste(plate.name, "markers", sep="."))
  
  panel.fn <- file.path(path.expand(phased.bgl.dir),
                        paste("chr", chr, ".phased.bgl", sep=""))
  WriteGtProbs(chr, obs.gt.fn, allele.tab, dbSNP.annot, verbose=verbose)  
  WriteMarkers(dbSNP.annot, markers.fn)
  
  phased.fn <- file.path(tmp.dir, paste(plate.name, plate.name, "GT", "phased",
                                        "gz", sep="."))
  gprob.fn <- file.path(tmp.dir, paste(plate.name, plate.name, "GT", "gprobs",
                                       "gz", sep="."))
  stdout.fn <- file.path(tmp.dir, paste(plate.name, "BEAGLE", "stdout", sep="."))

  gzFiles <- c(phased.fn, gprob.fn)
  
  i <- 1
  fail <- FALSE

  beagle.fn <- file.path("java", "beagle.jar")
  beagle.loc <- system.file(beagle.fn, package="HAPSEG")

  if (file.exists(beagle.loc) == FALSE) {
    stop("BEAGLE jar does not exist")
  }
  
  while (TRUE) {
    ## call BEAGLE
    beagle.tmp.dir <- file.path(tmp.dir, "BEAGLE")
    dir.create(beagle.tmp.dir, recursive=TRUE)
    ## FIXME: Look into calling this via rJava
    beagle.call <- paste("java -Xmx3500m ",
                         "-Djava.io.tmpdir=", beagle.tmp.dir,
                         " -jar ", beagle.loc, " like=",
                         obs.gt.fn, " phased=", panel.fn,
                         " lowmem=true markers=", markers.fn,
                         " missing=? out=", file.path(tmp.dir, plate.name),
                         " nimputations=", kNImp,
                         " niterations=", kNIter, " > ", stdout.fn, sep="")
    system(beagle.call)

    if (all(file.exists(gzFiles))) {
      break
    }

    if (verbose) {
      print(paste("RETRY BEAGLE call #", i, sep=""))
    }

    i <- i + 1
    if (i > 5) {
      fail <- TRUE 
      break
    }
  }
  
  if (!fail) {
    ## read BEAGLE output
    new.allele.tab <- ReadBeagleOutput(phased.fn, gprob.fn, dbSNP.annot)
    
    ## remove misc BEAGLE output files
    RemoveBeagleFiles(phased.fn, gprob.fn, stdout.fn, plate.name,
                      chr, tmp.dir)
  } else {
    new.allele.tab <- NA 
  }

  for (rm.fn in c(obs.gt.fn, markers.fn,
                  file.path(tmp.dir, paste(plate.name, "log", sep=".")))) {
    file.remove(rm.fn)
  }

  return(list(phase.snp.clust.p=new.allele.tab, fail=fail))
}

RemoveBeagleFiles <- function(phased.fn, gprob.fn, stdout.fn, plate.name,
                              chr, tmp.dir) {
  chrStr <- paste("chr", chr, sep="")
  rmFiles <- c(phased.fn, gprob.fn, stdout.fn,
               file.path(tmp.dir, paste(plate.name, plate.name, "GT", "r2", sep=".")),
               file.path(tmp.dir, paste(plate.name, chrStr, "phased", "bgl",
                                        "phased", "gz", sep=".")))
  for (rm.fn in rmFiles) {
    file.remove(rm.fn)
  }

  return(TRUE)
}

ReadBeagleOutput <- function(phased.fn, gprob.fn, dbSNP.annot) {
  beagle.output <- read.delim(phased.fn, row.names=NULL, sep=" ",
                              colClasses=c("character"))
  snp.ids <- beagle.output[, 2]
  beagle.output <- beagle.output[, c(3:ncol(beagle.output))]
  n.segs <- ncol(beagle.output) / 2
  
  gprob <- read.delim(gprob.fn, row.names=NULL, sep=" ",
                      colClasses=c(rep("character", 3),
                        rep("numeric", 3 * n.segs)))
  ## rs IDs are the same as the .phased file, so remove them from gprobs
  gprob <- gprob[, c(2:ncol(gprob))]
  
  snp.ix <- match(snp.ids, dbSNP.annot[, "dbSNP RS ID"])
  
  beagle.output <- beagle.output[snp.ix[!is.na(snp.ix)], ]
  gprob <- gprob[snp.ix[!is.na(snp.ix)], ]
  n.snp <- nrow(beagle.output)
  
  allele.tab <- array(NA, dim=c(n.segs, n.snp, 4))
  
  ac <- which(colnames(dbSNP.annot) == "Allele A")
  bc <- which(colnames(dbSNP.annot) == "Allele B")
  
  for (i in 1:n.segs) {
    snp.clust.p <- matrix(0, nrow=n.snp, ncol=4)

    ## BB, BA, AB, AA   
    ix <- (i - 1) * 2 + 1
    gt <- beagle.output[, c(ix:(ix + 1))]
    gp <- gprob[, 2 + c(ix:(ix + 2))]
    
    ## read phase for het SNPs
    ab.ix <- (dbSNP.annot[, ac] == gt[, 1]) & (dbSNP.annot[, bc] == gt[, 2])
    ba.ix <- (dbSNP.annot[, bc] == gt[, 1]) & (dbSNP.annot[, ac] == gt[, 2])

    ## BB
    snp.clust.p[, 1] <- gp[, 3]  
    ## AA
    snp.clust.p[, 3] <- gp[, 1] 
    ## het
    snp.clust.p[, c(2, 4)] <- gp[, 2]
    ## BA
    snp.clust.p[ba.ix, 2] <- gp[ba.ix, 2]
    snp.clust.p[ba.ix, 4] <- 0
    ## AB
    snp.clust.p[ab.ix, 4] <- gp[ab.ix, 2]
    snp.clust.p[ab.ix, 2] <- 0
    
    allele.tab[i, , ] <- snp.clust.p 
  }

  return(allele.tab)
}

WriteGtProbs <- function(chr, out.fn, allele.tab, dbSNP.annot, verbose=FALSE) {
  samples <- dimnames(allele.tab)[[1]]

  if (verbose) {
    print(nrow(dbSNP.annot))
  }
  
  ## dbSNP and A/B allele base
  snp.dat <- dbSNP.annot[ ,c("dbSNP RS ID", "Allele A", "Allele B")]
  n.row <- dim(allele.tab)[1]
  
  cols <- c("ID", "allele_A", "allele_B")
  
  for (i in 1:n.row) {
    dat <- allele.tab[i, , ]    
    hgt <- matrix(0, nrow=nrow(dat), ncol=3)
    hgt[, 2] <- rowSums(dat[, c(2, 4)], na.rm=TRUE)
    hgt[, 1] <- dat[, 3]
    hgt[, 3] <- dat[, 1]
    
    if (sum(hgt) > 0) {
      snp.dat <- cbind(snp.dat, hgt)
    }
    
    cols <- c(cols, rep(i, 3))
  }
  
  snp.dat <- snp.dat[rowSums(snp.dat[, c(4:ncol(snp.dat))]) > 0, ] 
  snp.dat[, c(4:ncol(snp.dat))] <- round(snp.dat[, c(4:ncol(snp.dat))], 4)
  
  colnames(snp.dat) <- cols
  
  return(write.table(snp.dat, file=out.fn, sep=" ", quote=FALSE,
                     row.names=FALSE))
}   

WriteMarkers <- function(dbSNP.annot, out.fn) {
  dat <- dbSNP.annot[, c("dbSNP RS ID", "pos", "Allele A", "Allele B")]
  return(write.table(dat, file=out.fn, sep=" ", quote=FALSE, row.names=FALSE,
                     col.names=FALSE))
}
