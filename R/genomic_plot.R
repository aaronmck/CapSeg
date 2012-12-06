# plot data that is arranged genomicly
plot_genome_data <- function(contig,start.pos,stop.pos,value,sample,segs=NA,is.hg18=FALSE,max.plot=3.0,inflate.axis=1.0,label.inflation=1.0) {
	unique_contigs <- unique(contig)
	max_pos = 0
    max.plot = log2(max.plot)
    value[value > max.plot] = max.plot
	# find the maximum position
	for (unique_contig in unique_contigs) { max_pos <- max_pos + max(stop.pos[which(contig==unique_contig)]) }
	# plot the initial plot
	par(mar=c(5,5,4,2) + 0.1)   
	plot(0,ylim=c(-0.5,3.5),xlim=c(0,max_pos),main=sample,ylab="Tangent Normalized CR",xlab="Chromosome",col=rgb(0,0,0,0.0),axes=F,bty="n",cex.lab=label.inflation)
    # check if theres a difference between the contig naming
    if (is.hg18 & !is.na(segs)) {
    	print("Fixing hg18 names")
		  segs[,2] = paste("chr",as.character(segs[,2]),sep="")
    }

	cur_pos = 0
	contig_midpoints = c()
	cur_col = "green"
	for (unique_contig in unique_contigs) {
            points(x=start.pos[contig==unique_contig & value < max.plot]+cur_pos,y=2^(value[contig==unique_contig & value < max.plot]),pch=18,col=rgb(0,0,1,0.03))
            points(x=start.pos[contig==unique_contig & value >= max.plot]+cur_pos,y=2^(value[contig==unique_contig & value >= max.plot]),pch=18,col=rgb(1,0,0,0.09))
            segments(x0=cur_pos+(max(stop.pos[which(contig==unique_contig)]))*0.1,x1=cur_pos + (max(stop.pos[which(contig==unique_contig)]))*0.9,y0=0,lwd=4,col=cur_col)

            # if we have plot information, get it
            if (!is.na(segs)) {
                seg.regions = segs[segs[,2]==unique_contig,]
                segments(x0=seg.regions$Start+cur_pos,x1=seg.regions$End+cur_pos,y0=2^(seg.regions$Segment_Mean),y1=2^(seg.regions$Segment_Mean),col=rgb(1,0,0,0.9),lwd=2)
            }
            if (cur_col == "green") { cur_col="black"} else {cur_col="green"}
            contig_midpoints = c(contig_midpoints,cur_pos + 0.5*max(stop.pos[which(contig==unique_contig)]) )
            cur_pos <- cur_pos + max(stop.pos[which(contig==unique_contig)])
	}
	axis(at=c(0.0,0.5,1.0,1.5,2.0,2.5,3.0),labels=T,side=2,cex.axis=inflate.axis)
	axis(at=contig_midpoints,labels=unique_contigs,side=1,tick=F,line=F,outer=F,cex.axis=inflate.axis)
        abline(h=0.5,col=rgb(0,0,0,0.1))
        abline(h=1.5,col=rgb(0,0,0,0.1))
        abline(h=1.0,col=rgb(0,0,0,0.1))
        abline(h=2.0,col=rgb(0,0,0,0.1))
        abline(h=2.5,col=rgb(0,0,0,0.1))
        abline(h=3.0,col=rgb(0,0,0,0.1))
}
