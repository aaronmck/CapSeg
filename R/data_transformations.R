#' cut down to intersecting baits
#' @param normal.data the normal data matrix
#' @param intersecting.baits
#' @return the pseudo inverse (in log 2 space) of the normal matrix
#' @keywords exome coverage tangent pseudo inverse
#'
intersect.tumor.normal.targets = function(normal.data,intersecting.baits) {
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
        print(paste("looking at file ",file))
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
        if (!is.vector(save.data$log.normals) && ncol(save.data$log.normals) > 1) {
            inter <- intersect(rownames(our.normals),rownames(save.data$log.normals))
            our.normals <- cbind(our.normals[is.element(rownames(our.normals),inter),],save.data$log.normals[is.element(rownames(save.data$log.normals),inter),])
        }
    }
    print(dim(our.normals))
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
        print(paste("examining file ",file))

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
