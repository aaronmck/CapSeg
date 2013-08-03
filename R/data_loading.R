#' load the exome data and process it down to a data matrix
#' @param data the data csv file name
#' @param baits the baits name
#' @return a data matrix of the coverage
#' @keywords exome coverage
#'
load.exome.data = function(data, baits) {
	return(read.delim(data,check.names=F,stringsAsFactors=F))
}

#' QC normal samples
#' @param normal.data a normal data matrix, with coverage for an indentical set of lanes as the tumor data
#' @param tumor.data the tumor data matrix
#' @return a qc'ed normal matrix, with failed samples removed 
#' @keywords exome qc coverage failed
#'
qc.normal.samples = function(normal.data, tumor.data, baits.filtered, qc.report.directory, allowed.normal.dev=0.3, allowed.normal.arm.dev=0.5, autosomes=as.character(seq(1,20))) {
  
  normal.mean = as.data.frame(cbind(apply(normal.data,2,mean),"NORMAL"),stringsAsFactors=F)
  rownames(normal.mean) = paste(colnames(normal.data),"Normal",sep=".")
  
  tumor.mean = as.data.frame(cbind(apply(tumor.data,2,mean),"TUMOR"),stringsAsFactors=F)
  rownames(tumor.mean) = paste(colnames(tumor.data),"Tumor",sep=".")
  
  total.mean = as.data.frame(rbind(normal.mean,tumor.mean),stringsAsFactors=F)
  colnames(total.mean) <- c("mean","type")
  rownames(total.mean) <- c(paste(colnames(normal.data),"Normal",sep="."),paste(colnames(tumor.data),"Tumor",sep="."))
  total.mean = transform(total.mean, mean = as.numeric(mean))
  
  normal.cr.plt = ggplot(total.mean,aes(x=mean,y=rownames(total.mean),col=type)) + geom_point() + theme_bw()
  ggsave(normal.cr.plt,file=paste(qc.report.directory,"normal_copy_ratio.jpg",sep="/"))
  
  # now exclude any normal where it's CR value is way out of line, one plus or minus allowed.normal.dev
  normal.excluded = total.mean[total.mean$type=="NORMAL" & abs(1.0 - total.mean$mean) > allowed.normal.dev, ]
  normal.data <- normal.data[,!is.element(paste(colnames(normal.data),"Normal",sep="."),rownames(normal.excluded))]
  
  # now mean center the data for the arm level check -- if we see greater deviation that expected at the arm level, drop that sample
  normal.data = sweep(normal.data,2,apply(normal.data,2,mean),"/")
  
  # now exclude any normals where the arm level copy ratio of any arm exceeds our threshold
  total.mean <- cbind(total.mean,0.0)
  colnames(total.mean) <- c("mean","type","arm.dev")
  for (i in autosomes) {
    chrome.baits <- baits.filtered$name[baits.filtered$contig==i]
    sample.dev <- apply(normal.data[is.element(rownames(normal.data),chrome.baits),],2,mean)
    #print(i)
    #print(sample.dev)
    for (i in names(sample.dev)) {
      #print(sample.dev[i])
      if (abs(1.0 - sample.dev[i]) > total.mean[which(rownames(total.mean)==paste(i,"Normal",sep=".")),"arm.dev"]) {
        #print(paste("replace ",total.mean[which(rownames(total.mean)==paste(i,"Normal",sep=".")),"arm.dev"],sample.dev[i]))
        total.mean[which(rownames(total.mean)==paste(i,"Normal",sep=".")),"arm.dev"] = abs(1.0 - sample.dev[i])
      }
    }
  }
  
  for (i in autosomes) {
    chrome.baits <- baits.filtered$name[baits.filtered$contig==i]
    sample.dev <- apply(tumor.data[is.element(rownames(tumor.data),chrome.baits),],2,mean)
    #print(i)
    #print(sample.dev)
    for (i in names(sample.dev)) {
      #print(sample.dev[i])
      if (abs(1.0 - sample.dev[i]) > total.mean[which(rownames(total.mean)==paste(i,"Tumor",sep=".")),"arm.dev"]) {
        #print(paste("replace ",total.mean[which(rownames(total.mean)==paste(i,"Tumor",sep=".")),"arm.dev"],sample.dev[i]))
        total.mean[which(rownames(total.mean)==paste(i,"Tumor",sep=".")),"arm.dev"] = abs(1.0 - sample.dev[i])
      }
    }
  }
  arm.movement = ggplot(total.mean,aes(arm.dev,fill=type)) + geom_density(alpha=0.8)
  ggsave(normal.cr.plt,file=paste(qc.report.directory,"normal_arm_ratio.jpg",sep="/"))
  
  # throw out any normals with arm level events that deviate from CR 1 by more than 0.5
  
  normal.excluded = total.mean[total.mean$type=="NORMAL" & total.mean$arm.dev >= allowed.normal.arm.dev, ]
  normal.data <- normal.data[,!is.element(paste(colnames(normal.data),"Normal",sep="."),rownames(normal.excluded))]
  return(normal.data)
}

# for the samples with normals, guess the sex of the patient.  If there's no normal, guess from the tumor
#' @param tumor.data the tumor matrix
#' @param normal.data the normal matrix
#' @param sample.table the normal matrix
#' @param filtered.baits the list of baits that are in the tumor and the normal
#' @param male.chr the male chromosome
#' @param female.chr the female chromosome
#' @return a sex assignment for each individual, with some information about which bam it was called from and the male to female ratio
#' @keywords sex assignment 
#'
create.sex.assignments <- function(tumor.data,
                                   normal.data,
                                   sample.table,
                                   filtered.baits,
                                   male.chr="Y",
                                   female.chr="X") {
  male.coverage.targets = filtered.baits$name[filtered.baits$contig == male.chr]
  female.coverage.targets = filtered.baits$name[filtered.baits$contig == female.chr]
  
  # for each sample find the male/female coverage.  In humans this should be either 1/1 or 0/2 (ish).  We
  # can then split or cluster to determine sex
  sex.assignments <- c()
  normal.male.coverage <- apply(normal.data[male.coverage.targets,],2,mean)
  normal.female.coverage <- apply(normal.data[female.coverage.targets,],2,mean)
  tumor.male.coverage <- apply(tumor.data[male.coverage.targets,],2,mean)
  tumor.female.coverage <- apply(tumor.data[female.coverage.targets,],2,mean)
  
  # remove any duplicated entries
  normal.male.coverage <- normal.male.coverage[!duplicated(names(normal.male.coverage))]
  normal.female.coverage <- normal.female.coverage[!duplicated(names(normal.female.coverage))]
  tumor.male.coverage <- tumor.male.coverage[!duplicated(names(tumor.male.coverage))]
  tumor.female.coverage <- tumor.female.coverage[!duplicated(names(tumor.female.coverage))]
  
  for (ind in sample.table$sample) {
    sample.ratio = 1.0 # our default will be in the female range
    call.type = "UNK"
    if (is.element(ind,colnames(normal.data))) {
      sample.ratio <- normal.male.coverage[ind] / normal.female.coverage[ind]
      call.type = "NORMAL"
    }
    else if (is.element(ind,colnames(tumor.data))) {
      sample.ratio <- tumor.male.coverage[ind] / tumor.female.coverage[ind]
      call.type = "TUMOR"
    }
    sex.assignments <- rbind(sex.assignments,c(ind,sample.ratio,sample.ratio<0.5,call.type))
  }
  ret = as.data.frame(sex.assignments)
  colnames(ret) <- c("sample","male.to.female","is.female","call.bam.source")   
  return(ret)
}

