## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

FullGenPlots <- function(seg.dat, plots.dir) {
  plot.fn = file.path( plots.dir, "FullGenome.jpeg") 
  jpeg(plot.fn, 12, 7, units="in", type="cairo", res=500, quality=100)
  par( mfrow=c(2,1), cex=0.75, las=1 )

  chr.dat = SetUpPlot("Total copy ratio", 0, 5, "", F)
  TotalCopyRatio(seg.dat, chr.dat)

  SetUpPlot("Fraction of alternate reads", 0, 1, "Chromosome", T)
  AltReadFrac(seg.dat, chr.dat)

  dev.off()
}

#Chromosome Background for data
SetUpPlot <- function(y.lab, y.min, y.max, x.lab, lab.chr){
  ## Via http://genome.ucsc.edu/cgi-bin/hgTables 
  ## HG19 Chr Lens
  chr.lens <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
        146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
        102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566)

  ## HG19 Cent Pos
  cent.pos <- c(121535434, 92326171, 90504854, 49660117, 46405641, 58830166, 58054331,
        43838887, 47367679, 58632012, 10104553, 39254935, 51644205, 34856694, 16000000,
        16000000, 17000000, 35335801, 22263006, 15460898, 24681782, 26369569, 11288129,
        13000000)
  chr.w <- chr.lens / sum(chr.lens)
  par( mar=c(3.1,3.6,0.1,0),mgp=c(2,-0.2,-1.1) )
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(y.min, y.max),
       xlab = "", ylab = "", main = "", xaxt = "n")
  mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
  mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)
  ww <- as.vector(rbind(chr.w, chr.w)) / 2
  chr.mids <- cumsum(ww)[(c(1:length(ww)) - 1) %% 2 == 0]
  if(lab.chr){  
  	lab.vals <- (c(1:length(chr.w)))
  	odd.ix <- lab.vals %% 2 == 1
    
  	mtext(text = lab.vals[odd.ix], side = 1, line = -0.45, at = chr.mids[odd.ix], 
    		las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  	mtext(text = lab.vals[!odd.ix], side = 1, line = 0, at = chr.mids[!odd.ix],
      	las = 1, cex = par("cex.axis") * par("cex") * 0.7)
  }  
  chr.offsets <- c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1)
  cent.pos <- (cent.pos)/sum(chr.lens) +
  	chr.offsets[c(1:(length(chr.offsets) - 1))]
    
  for (i in 1:(length(chr.offsets) - 1)) {
    use.col <- ifelse(i%%2 == 1, "grey90", "white")
    rect(xleft = chr.offsets[i], ybottom = y.min, xright = chr.offsets[i + 1], 
    		ytop = y.max, col = use.col, border = NA)
    lines(y = c(y.min, y.max), x = rep(cent.pos[i], 2), lty = 3, lwd = 0.5)
  }
  chr.dat = list(chr.w, chr.lens, chr.offsets)
  return( chr.dat )
}

## Plotting Copy Ratio 
TotalCopyRatio <- function(seg.dat, chr.dat){
  chr.w = chr.dat[[1]]
  chr.lens = chr.dat[[2]]
  chr.offsets = chr.dat[[3]]
  #stop()
  for( s in 1:length(seg.dat[["as.res"]][["h.seg.dat"]][[1]])) {
	 Theta = seg.dat[["capture.em.fit"]][["Theta"]]
   	 delta.tau = seg.dat[["capture.em.fit"]][["delta.tau"]]
	 atten.tau1 = AffyAtten(delta.tau[s, 2], Theta[["at.capseg"]])
	 chr <- seg.dat$as.res$h.seg.dat$h.capseg.annot[[s]][["chr"]]
	 short <- seg.dat[["as.res"]][['h.seg.dat']]
	 seg.crds <- short[["h.capseg.annot"]][[s]][["pos"]]
	 genome.crds <- chr.offsets[chr] + seg.crds / chr.lens[chr] * chr.w[chr]
	 
	 colors=c("coral", "dodgerblue")

	 points(genome.crds,short[["h.capseg.d"]][[s]][1,], col=colors[s%%2+1], pch=16, cex=0.2)
	 lines(range(genome.crds), rep(atten.tau1, 2), col="black", lwd=1, lty=2)
  }
}

## Plotting Alternate Read Fraction
AltReadFrac <- function(res, chr.dat){
   chr.w = chr.dat[[1]]
   chr.lens = chr.dat[[2]]
   chr.offsets = chr.dat[[3]]
   for( i in 1:length(res[["as.res"]][["h.seg.dat"]][[1]])) {
      h.seg.dat = res[["as.res"]][["h.seg.dat"]]
      Theta = res[['capture.em.fit']][["Theta"]]
      d = h.seg.dat[["gh.wes.allele.d"]][[i]]
      if(length(d)!=0){
	 chr = res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][["chr"]]
	 f.hat = res[["capture.em.fit"]][["wes.f"]][i,"f.hat"]
	 het.mat = res[["capture.em.fit"]][["het.phase.prob"]][[i]]
	 het.out = het.mat[,3]

	 het.phase.prob = het.mat[,1] / (1 - het.out) # renormalize to marg out outlier state

	 het.crds = h.seg.dat[["gh.wes.allele.annot"]][[i]][['pos']]
	 het.crds <- chr.offsets[chr] + het.crds / chr.lens[chr] * chr.w[chr]

	 seg.crds = res[["as.res"]][['h.seg.dat']][["h.capseg.annot"]][[i]][["pos"]]
	 seg.crds <- chr.offsets[chr] + seg.crds / chr.lens[chr] * chr.w[chr]

	 alt= d["alt",] 
	 ref= d["ref",] 

	 het.1.color <- "red"
	 het.2.color <- "blue"
	 mid.color <- "darkviolet"
	 pal <- colorRampPalette(c(het.1.color, mid.color, het.2.color))
	 cols <- pal(1000)
	 het.col =  cols[999 * het.phase.prob + 1]
	 out.ix = het.out > 0.5
	 het.col[out.ix] = "black"
	 use.pch = rep(16, length(het.out))
	 use.pch[out.ix] = 1

	 cov = alt+ref

	 het_CI_low = qbeta( 0.025, alt+1, ref+1 )
	 het_CI_high = qbeta( 0.975, alt+1, ref+1 )
	 
	 ## Plot Error Bars
	 #for( i in 1:length(het.crds)){
	 #    if( !out.ix[i] ) {next}
	 #	lines( x=rep(het.crds[i],2), y=c(het_CI_low[i],het_CI_high[i]), 
	 #		col=het.col[i], lwd=0.3)
	 #}

	 points( het.crds[!out.ix], (alt/cov)[!out.ix], 
			col=het.col[!out.ix], pch=use.pch[!out.ix], cex=0.2)
	 ## Plot Error Bars
	 #for( i in 1:length(het.crds)){
	 #    if( out.ix[i] ) {next}
	 #	lines( x=rep(het.crds[i],2), y=c(het_CI_low[i],het_CI_high[i]), 
	 #		col=het.col[i], lwd=0.3)
	 #}

	 points( het.crds[out.ix], (alt/cov)[out.ix], 
			col=het.col[out.ix], pch=use.pch[out.ix], cex=0.075)
		 
	 A1 = Theta[["f_skew"]]*f.hat
  	 A2 = 1 - (Theta[["f_skew"]]^-1 * f.hat)
  	 lines(range(seg.crds), rep(A1, 2), lwd=1, lty=2)
  	 lines(range(seg.crds), rep(A2, 2), lwd=1, lty=2)
	 } 
   }  
}
