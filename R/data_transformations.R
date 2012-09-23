library(plyr)

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
  baits
  ) {

  # this next part is a bit hacky, but its a speedup and a huge memory win over the default read.csv
  line.count = as.numeric(system(paste("wc -l",data,"| cut -f 1 -d ' '"),intern = TRUE))
  col.count = as.numeric(system(paste("head -1",data,"| tr ',' '\n' | wc -l"),intern = TRUE))
  col.types = c("character",rep("numeric", times=col.count-1))

  # load up the raw data - there shouldn't be duplicate row names anymore
  raw <- read.table(data,sep=",",row.names=1,header=TRUE,stringsAsFactors=F,check.names=F,nrows=line.count,comment.char = "",colClasses=col.types)

  # raw2 <- join(baits,raw,by = (baits$name,rownames(raw)),type="left")
  return(raw)
}

#' given a coverage matrix, remove all of the bad lanes, merge them to samples, and remove bad baits
#' @param lane.data the lane coverage data, normalized and calibrated
#' @param samples.to.lanes information about the sample names and their associated lanes
#' @param bait.factor the bait factor, the median coverage per bait (or target) in the pull-down
#' @param report.dir where to put our PDF with the eCDFs of CR.stats
#' @param log.bait.factor.cut any bait below this threshold is removed from the bait matrix
#' @param cr.stat.lane.cut if a lanes copy ratio statistic is above this threshold, cut it from the matrix
#' @param max.lane.pct.cut the maximum percentage of the lanes we're allowed to cut from a sample (given the cut threshold); if more exceed this threshold we keep them anyway
#' @return a data matrix of the coverage collapsed down to samples (from lanes), and bad baits removed
#' @keywords exome coverage
#'
collapse.exome.coverage = function(lane.data,samples.to.lanes,bait.factor,report.dir="UNK",log.bait.factor.cut=-12.0,cr.stat.lane.cut=0.45,max.lane.pct.cut=0.5,parallel=FALSE) {
  sample.data <- NA
  print("Calculating the CR stat")
  cr.stats <- apply(lane.data,2,function(x) {return(mean(abs(diff(x))))})
  data.colnames <- NULL

  # if we're running in parallel, make sure to include the library, and use the foreach to execute
  if (parallel) {
    print("Running parallel...")

    # setup the multithreaded libraries
    library(foreach)
    library(doMC)
    registerDoMC()
    print(getDoParWorkers())

    # now ask foreach to parallelize out the calculation
    sample.data <- foreach(tumor=unique(samples.to.lanes[,1]), .combine='data.frame') %dopar% aggregate.sample(lane.data,tumor,samples.to.lanes,cr.stats,quantile(cr.stats)["75%"],max.lane.pct.cut)

  } else {
    print("Running single-threaded...")
    # now for each sample, find the lanes that comprise the coverage, remove bad lanes, and merge back the results
    for (target.sample in unique(samples.to.lanes[,1])) {
      if (is.na(sample.data)) {
        aggregate.sample(lane.data,target.sample,samples.to.lanes,cr.stats,cr.stat.lane.cut,max.lane.pct.cut)
      } else {
        sample.data <- data.frame(sample.data,aggregate.sample(lane.data,target.sample,samples.to.lanes,cr.stats,quantile(cr.stats)["75%"],max.lane.pct.cut))
      }
    }
  }
  return(sample.data)
}



#' given the coverage for a sample (a collection of lanes) and the CR.stat cutoff, throw out bad
#' lanes, returning the merged coverage as a vector of baits (or targets)
#'
#' @param lane.data the coverage data
#' @param cr.stats the copy ratio stats
#' @param target.sample the sample to merge
#' @param samples.to.lanes the mapping of samples to their lanes
#' @param cr.stat.lane.cut the copy ratio stat cutoff to remove lanes
#' @param max.lane.pct.cut the maximum fraction (as a value between 0 and 1) of lanes that we're willing to drop from a sample
#'
#' @return the coverage of lanes that passed QC, merged to a sample level value
#'
#' @keywords aggregate sample lane cr.stat
#'
aggregate.sample <- function(lane.data,target.sample,samples.to.lanes,cr.stats,cr.stat.lane.cut,max.lane.pct.cut) {
    print(target.sample)
  specific.lanes <- samples.to.lanes[which(samples.to.lanes[,1]==target.sample),2]
  data.lanes <- data.frame(lane.data[,match(specific.lanes,colnames(lane.data))])
  sample.cr.stat <- cr.stats[match(specific.lanes,colnames(lane.data))]
  data.lanes <- data.frame(data.lanes[,!colSums(is.na(data.lanes)) > 0])
  sample.cr.stat <- sample.cr.stat[!is.na(sample.cr.stat)]
  # if we have more than a single lane
  if (length(specific.lanes) > 1 & nrow(data.lanes) > 1) {
      #print(sample.cr.stat)
      #print(cr.stat.lane.cut)
    cr.stat.not.exceeded <- sum(sample.cr.stat <= cr.stat.lane.cut)

    data.lanes = data.lanes[,order(sample.cr.stat,decreasing=F)]
    #print(dim(data.lanes))
    sample.cr.stat = sample.cr.stat[order(sample.cr.stat)]
    #print(ncol(data.lanes))
    #print(cr.stat.not.exceeded)
    #print(max.lane.pct.cut)
    if ((ncol(data.lanes)-cr.stat.not.exceeded)/ncol(data.lanes)>max.lane.pct.cut) {
      print(paste("(MAX-THROW) Threw away lanes with a median CR.stat of",median(cr.stats[seq(max.lane.pct.cut*ncol(data.lanes)+1,ncol(data.lanes))])))
      data.lanes = data.frame(data.lanes[,seq(1,ceiling(max.lane.pct.cut*ncol(data.lanes)))])
    } else {

      print(paste("Threw away lanes with a median CR.stat of",median(cr.stats[seq(cr.stat.not.exceeded+1,ncol(data.lanes))])))
      data.lanes = data.lanes[,seq(1,cr.stat.not.exceeded)]
    }
    print(paste("sample",target.sample,"exceeded =",ncol(data.lanes)-cr.stat.not.exceeded,", size",nrow(data.lanes),"rows by",ncol(data.lanes),"columns"))

  }
  if (length(specific.lanes) >= 1 & nrow(data.lanes) > 1) {
    ret = data.frame(apply(data.lanes,1,median))
    colnames(ret) = c(target.sample)
    # print(dim(ret))
    return(ret)
  } else {
    print("returning NULL!")
    return(NULL);
  }
}


#' cut down to intersecting baits
#' @param normal.data the normal data matrix
#' @param intersecting.baits
#' @return the pseudo inverse (in log 2 space) of the normal matrix
#' @keywords exome coverage tangent pseudo inverse
#'
normalize.tumor.normal.targets = function(normal.data,intersecting.baits) {
  rn <- rownames(normal.data)
  cl <- colnames(normal.data)
  subset = is.element(rownames(normal.data),intersecting.baits)
  rn <- rn[subset]
  normal.data = data.frame(normal.data[subset,])
  rownames(normal.data) <- rn
  colnames(normal.data) <- cl
  return(normal.data)
}


#' given a pseudo inverse matrix and a case matrix, return the multiplication
#' @param log.tumor the tumor data in log space
#' @param pseudo.inverse the inverse matrix
#' @param normals.cr the copy ratio calibrated normal data in log2 space
#' @return the pseudo inverse (in log 2 space) of the normal matrix
#' @keywords exome coverage tangent pseudo inverse
#'
calibrate.tumors <- function(log.tumor,pseudo.inverse,log.normals.cr, center, first=TRUE) {
	# let them know when we're starting the run
	# print(paste("starting the normalization at",format(Sys.time(), "%a %b %d %H:%M:%S %Y")))

        matchings <- matrix(NA,ncol=ncol(log.normals.cr),nrow=ncol(log.tumor))
        tumor.matrix <- matrix(NA,ncol=ncol(log.tumor),nrow=nrow(log.tumor))
        dimnames(tumor.matrix) <- dimnames(log.tumor)

	log.median.tumor <- apply(log.tumor,2,median)
        colnames(matchings) <- colnames(log.normals.cr)
        rownames(matchings) <- colnames(log.tumor)

	# now run the target tumor lanes, and see what happens
	for (i in seq(1,ncol(log.tumor))) {
		# output where we are
		cat(paste(i,",",sep=""))

                #med.log.tumor <- median(log.tumor[,i])

		# now get a cancer sample

                if (first)
                    tumor.medianed <- log.tumor[,i]#  - center
                else
                    tumor.medianed <- log.tumor[,i]

		# find the f value
                per.sample = t(tumor.medianed) %*% t(pseudo.inverse)
		f.value <- log.normals.cr %*% t(per.sample)

		# do the subtraction
		result <- t(t(tumor.medianed - (f.value)))

		tumor.matrix[,i] <- result
	}

	# ending the normalization run
	print(paste("Ending the normalization at",format(Sys.time(), "%a %b %d %H:%M:%S %Y")))
	return(tumor.matrix)
}



#' output and plot the data
#' @param tumor.matrix the post calibrated tumor matrix data
#' @param tumors.cr the pre-calibrated, log2 normalized CR data for the tumor
#' @param baits the baits list
#' @param output.location where to write the data to
#' @keywords output plot cr
#'
output.and.plot.data <- function(tumor.matrix,tumors.cr,baits,output.location,signal.files) {
  # now output and plot the data for each sample
  print("Writing out the data and plotting each sample...")
  raw.dir.name <- paste(output.location,"/signal/",sep="")
  dir.create(raw.dir.name,recursive=T)
  plot.dir.name <- paste(output.location,"/plots/",sep="")
  plot.nt.dir.name <- paste(output.location,"/plots_nt/",sep="")
  dir.create(plot.dir.name,recursive=T)
  dir.create(plot.nt.dir.name,recursive=T)

  # output the raw data
  iter.bait.names <- intersect(baits$name,rownames(tumor.matrix))
  if (length(iter.bait.names) < 1) {
    print("Unable to find an intersection between the target names.  Make sure that you have matching bait files (csv, interval_file, bed) in your setup.")
    q()
  }

  # open up the singal files
  sgn.files = read.delim(signal.files,check.names=F,stringsAsFactors=F)
  print(sgn.files)
  print(summary(tumor.matrix))
  for (name in colnames(tumor.matrix)) {
    if (sum(is.na(tumor.matrix[,name])) > 0) {
      print(paste("Skipping sample",name))
      next()
    }
    base.name <- paste(raw.dir.name,name,sep="/")

    # lookup the filename in sgn.files
    file.name = sgn.files[sgn.files$sample==name,"signal.file"]
    print(file.name)

    print(paste("Writing file name",file.name,"for tumor",name))
    #print(rownames(tumor.matrix)[1:100])
    output.data <- cbind(baits[is.element(baits$name,iter.bait.names),],tumor.matrix[is.element(rownames(tumor.matrix),iter.bait.names),name])
    #print(output.data[1:100,])
    write.table(output.data,file=file.name,quote=F,row.names=F,col.names=c(colnames(baits),name),sep="\t")
    base.name <- paste(plot.dir.name,name,sep="/")
    base.name.nt <- paste(plot.nt.dir.name,name,sep="/")
    file.name <- paste(base.name,"plot.png",sep=".")

    # file.name <- paste(,"plot.png",sep=".")
    file.name.nt <- paste(base.name.nt,"plot.png",sep=".")

    png(file.name,width=2000,height=700)
    # plot the data from the post tangent normalization
    plot_genome_data(output.data[,2],output.data[,3],output.data[,4],tumor.matrix[,name],name)
    dev.off()

    png(file.name.nt,width=2000,height=700)
    # plot the data from the pre-tangent normalization
    lg2 <- data.frame(log2(tumors.cr[,name]))
    colnames(lg2) <- name

    plot_genome_data(output.data[,2],output.data[,3],output.data[,4],lg2,name)
    dev.off()
  }
}

#' calculate the per-sample metrics, and output them to a single file
#' @param tumor.matrix the tumor matrix of calibrated data
#' @param log.tumor the logged (2) tumor data
#' @param output.location where to put the data
#' @keywords exome coverage metric
#'
create.sample.metric.files <- function(tumor.matrix,log.tumor,output.location) {
  # output metrics for each of the samples
  print("creating metrics for every sample...")
  dir.name <- paste(output.location,"/stats/",sep="")
  dir.create(dir.name,recursive=T)

  # compute the sample level stats
  sample.stats <- cbind(colnames(tumor.matrix),apply(tumor.matrix,2,function(x) { median(abs(diff(x))) }))
  colnames(sample.stats) <- c("sample","MAD")
  write.table(sample.stats,sep="\t",file=paste(dir.name,"samples.stats.txt",sep=""),quote=F,row.names=F)

  # compute the lane level stats
  lane.stats <- cbind(colnames(tumor.matrix),apply(tumor.matrix,2,function(x) { median(abs(diff(x))) }))
  colnames(lane.stats) <- c("lane","MAD")
  write.table(lane.stats,sep="\t",file=paste(dir.name,"lane.stats.txt",sep=""),quote=F,row.names=F)

  # compute the sample level stats
  sample.stats <- cbind(colnames(tumor.matrix),apply(log.tumor,2,function(x) { median(abs(diff(x))) }))
  colnames(sample.stats) <- c("sample","MAD")
  write.table(sample.stats,sep="\t",file=paste(dir.name,"samples.stats.nt.txt",sep=""),quote=F,row.names=F)

  # compute the lane level stats
  lane.stats <- cbind(colnames(tumor.matrix),apply(log.tumor,2,function(x) { median(abs(diff(x))) }))
  colnames(lane.stats) <- c("lane","MAD")
  write.table(lane.stats,sep="\t",file=paste(dir.name,"lane.stats.nt.txt",sep=""),quote=F,row.names=F)
}


# save off the data to a Rdata object for later
#' @param log.normals the normal data, uncalibrated, log2
#' @param log.tumors the tumor data, uncalibrated, log2
#' @param pseudo.inverse.norm, the pseudo inverse object
#' @param baits the bait information
#' @param cache.location where to store the data
#' @param analysis.set.name the analysis set name
save.off.processed.data <- function(log.normals,log.tumors,calibrated.tumors,baits,cache.location,analysis.set.name,build.version) {
  save.data <- list(log.normals=log.normals,log.tumors=log.tumors,calibrated.tumors=calibrated.tumors,baits=baits,analysis.set.name=analysis.set.name,version="5",build.version=build.version)
  output.file = paste(cache.location,"/",analysis.set.name,".v4.rData",sep="")
  save(save.data,file=output.file)
}

#' tangent normalize the tumors by the normals
#' @param normal.data the normal data matrix
#' @return the pseudo inverse (in log 2 space) of the normal matrix
#' @keywords exome coverage tangent pseudo inverse
#'

pseudo.invert.normals = function(log.normals.cr) {
  # any zeros become epsilons
  print("performing pseudo inverse")

}




#' make a big panel of normals out of the database that we have on disk
#' and use that to supliment our normals for the pseudo-inverse
#' @param normal.database where to find the normal files
#' @param our.normals this sets normals
#' @param build.ver the build version (hg19, hg18, mm9)
#' @param analysis.set the analysis set name
#' @param max.normals the max number of normals to use (the max column count for the matrix)
#' @param max.lost.markers.per.iteration the max percentange loss of targets per addition (the max we can lose each time we add a new matrix of normals)
#' @return a matrix containing the new normals
#' @keywords exome coverage tangent pseudo inverse
additional.normals = function(normal.database,our.normals,build.ver,analysis.set,max.normals=600,max.lost.markers.per.iteration=0.15) {
    normal.database = paste(normal.database,build.ver,sep="/")
    print(paste("looking for normal subspaces in the directory ",normal.database))

    our.normals <- data.frame(our.normals)
    for (file in list.files(normal.database)) {
        if(!grep("rData",file)) {
            print(paste("Skipping file",file,"since it's not an R data file"))
            next
        }
        # load up the normal -- load up the log'ed matrix
        fl = paste(normal.database,file,sep="/")
        load(fl)

        # find the overlaping names for the stored normal and the target normal matrix; if we're take more than a hit of max.lost.markers.per.iteration
        # percent of the markers, don't use this matrix'
        percent.overlap = sum(is.element(rownames(our.normals),rownames(save.data$log.normals)))/nrow(our.normals)

        # drop max.lost.markers.per.iteration
        if (1.0-percent.overlap > max.lost.markers.per.iteration) {
            print(paste(file,": The overlap between the normals and first normal vector was only",percent.overlap,",we require at least",1.0-max.lost.markers.per.iteration))
            next
        } else if (save.data$analysis.set.name==analysis.set) {
            print(paste("The set has the same analysis set name",analysis.set))
            next
        } else {
            print(paste("tangent calculation with",file,"since it has an overlap of ",percent.overlap))
        }

        if (ncol(save.data$log.normals)+ncol(our.normals)>max.normals) {
            save.data$log.normals = save.data$log.normals[,1:(max.normals-ncol(our.normals))]
        }

        # now combine the two
        if (ncol(save.data$log.normals > 0)) {
            inter <- intersect(rownames(our.normals),rownames(save.data$log.normals))
            our.normals <- cbind(our.normals[is.element(rownames(our.normals),inter),],save.data$log.normals[is.element(rownames(save.data$log.normals),inter),])
        }
    }
    return(our.normals)
}


#' a quick function to get the normal and the tumor names extracted
basename <- function(file) {
    sp = sapply(file,function(x) {nt = unlist(strsplit(x,"/")); return(nt[length(nt)])})
    return(sp)
}


#' given a location for other normal files, load additional normals
#'
#' @param normal.db.location where we can find our normal databases
#' @param set.normals the normals from our set
#' @param build.ver the build version (hg18, hg19, etc)
#' @param analysis.set the analysis set name
#' @param max.lost.markers.per.iteration how many markers we're willing to lose per iteration (max is plane.limit * max.list.markers)
#' @param sample.limit the max number of normal samples to include
load.additional.normals <- function(normal.db.location,
                                                    set.normals,
                                                    build.ver,
                                                    analysis.set,
                                                    max.lost.markers.per.iteration=0.10,
                                                    sample.limit=400) {

    # find each of the normals available - list the files in the directory
    normal.database = paste(normal.db.location,build.ver,sep="/")
    print(paste("looking for normal subspaces in the directory ",normal.database))

    historical.normals.distance <- list()
    median.of.samples <- apply(set.normals,1,median)

    # loop over all the normal files, looking for the closest sample.limit samples
    for (file in list.files(normal.database)) {
        # we have two reason for skippin a file - first that it's not an R data file,
        # and the second is if we have already processed two many planes
        if(!grep("rData",file)) {
            print(paste("Skipping file",file,"since it's not an R data file"))
            next
        }

        # load up the normal -- load up the log'ed matrix
        fl = paste(normal.database,file,sep="/")
        load(fl)

        # find the overlaping names for the stored normal and the target normal matrix; if we're take more than a hit of max.lost.markers.per.iteration
        # percent of the markers, don't use this matrix'
        percent.overlap <- sum(is.element(names(set.normals),rownames(save.data$log.normals)))/length(set.normals)

        # drop max.lost.markers.per.iteration
        if (1.0 - percent.overlap >= max.lost.markers.per.iteration) {
            print(paste("The overlap between the normals and original normal vector was only ",percent.overlap,", we require at least",max.lost.markers.per.iteration))
            next
        } else if (save.data$analysis.set.name==analysis.set) {
            print(paste("The set has the same analysis set name",analysis.set))
            next
        } else {
            print(paste("tangent calculation with",file,"since it has an overlap of ",percent.overlap))
        }
        print(paste("Using tangent plane in ",fl,"for calibration"))

        # add the first normal vector to the front
        baits.to.use <- intersect(names(first.nv),rownames(save.data$log.normals))
        loaded.normals <- save.data$log.normals[is.element(rownames(save.data$log.normals),baits.to.use),]
        set.normals <- set.normals[is.element(rownames(set.normals),baits.to.use),]

        #  now record the distance
        historical.normals.distance <- apply(loaded.normals - median.of.samples,2,function(x) {sum(abs(x))})
    }
    return(calibrated.tumors)
}

#' convert the names in the tumor matrix to the names we've passed in
#'
#' @param bamSampleListing the matching of bams, tumors and normals
#' @return the matching sample listings
#' @keywords exome coverage tangent pseudo inverse
#'
tumor.to.individual = function(tumor.bams.samples,normal.bams.samples,bam.sample.listing) {
    tbs = read.delim(tumor.bams.samples,header=T,check.names=F,stringsAsFactors=F)
    nbs = read.delim(normal.bams.samples,header=T,check.names=F,stringsAsFactors=F)
    bsl = read.delim(bam.sample.listing,header=T,check.names=F,stringsAsFactors=F)

    # now line up the tumor and normal bams with their samples
    ret <- cbind(tbs$Sample[match(basename(bsl$tumor_bam),tbs$BAM)],nbs$Sample[match(basename(bsl$normal_bam),nbs$BAM)])
    colnames(ret) <- c("tumor","normal")
    return(ret)
}
