
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
