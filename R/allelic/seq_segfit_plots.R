
PlotAllCapture <- function(res, save = F, plots.dir)
{
   CalcIntensityLims <- function(probe.type)
   {
       all.d = list()
       all.d[["capseg"]] = c(res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i]], res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i+1]])
       medians = lapply(all.d, median )

       if (!is.na(medians[[probe.type]]))
       {
          other.probe.types = setdiff(names(all.d), probe.type)
          d1 = abs(medians[[probe.type]] - medians[[other.probe.types[1] ]])
          d2 = abs(medians[[probe.type]] - medians[[other.probe.types[2] ]])
          if ( max(c(d1, d2)) > 30 | is.na(max(c(d1, d2))) ) {
             quant = quantile(all.d[[probe.type]], probs=seq(from=.01, 1, length.out=100))
             intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
          } else {
             quant = quantile(c(all.d[["capseg"]]), probs=seq(from=.01, 1, length.out=100))
             intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
          }
          if (intensity.lim[2] - intensity.lim[1] < 5 ) intensity.lim = c(mean(intensity.lim) - 2.5, mean(intensity.lim) + 2.5)
          if( intensity.lim[1] < 0 ) { intensity.lim[1] = 0 }  #scarter
          return(intensity.lim)
       } else {
          return(NULL)
       }
   }

   dir.create(plots.dir, recursive=TRUE)
   cat("Plotting segs:")
   for( i in 1:(length(res[["as.res"]][["h.seg.dat"]][[1]]) -1))
   {
      if (res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][["chr"]] != res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i+1]][["chr"]])
      {
         next
      }

      genomic.limits = range(c(res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i+1]][['pos']]))
      genomic.limits = genomic.limits / 1e6

      cat( paste(i, ".", sep=""))

      ilim = CalcIntensityLims("capseg")

      plot.fn = file.path( plots.dir, paste("Segs_", i,"_",i+1, ".jpeg", sep=""))
      jpeg(plot.fn, 7, 5, units="in", type="cairo", res=200, quality=100)

      par(mfrow=c(2,3))
      ScarterPar()
      PlotCapIntensities(i, res, genomic.limits, ilim, use.capseg.mean=F, draw.legend=F)
      PlotCapsegSegfit(i, res, d.col="coral", min=ilim[1], max=ilim[2], plot=T)
      PlotCapsegSegfit(i+1, res, d.col="dodgerblue", min=ilim[1], max=ilim[2], plot=T)

      Plot_allelic_fraction_vs_genome( i, res, genomic.limits )

      PlotFFit(i, res, conf=.95, plot=T )
      PlotFFit(i+1, res, conf=.95, plot=T )

#      Plot_het_AF_vs_cov(i, res)

      dev.off()
   }
}





PlotCapIntensities <- function(i, res, genomic.limits, intensity.lim, use.capseg.mean=T, draw.legend=T) {

   #########
   seg1.col = c("coral", "black")
   seg2.col = c("dodgerblue", "black")
   seg.chr = res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][["chr"]]
   Theta = res[["em.fit"]][["Theta"]]
   delta.tau = res[["em.fit"]][["delta.tau"]]

   pos1 = res[["as.res"]][['h.seg.dat']][["h.capseg.annot"]][[i]][["pos"]] / 1e6
   pos2 = res[["as.res"]][['h.seg.dat']][["h.capseg.annot"]][[i+1]][["pos"]] / 1e6
   int = c(res[["as.res"]][['h.seg.dat']][["h.capseg.d"]][[i]], res[["as.res"]][['h.seg.dat']][["h.capseg.d"]][[i+1]])
   seg12 = data.frame(intensity=int, position=c(pos1, pos2), col=c(rep(seg1.col[1], length(pos1)), rep(seg2.col[1], length(pos2))), stringsAsFactors=F)

   XLAB = paste("Chromosome ", seg.chr, " position (MB)", sep="")
   if (nrow(seg12) == 0) {
      plot(1, type="n", xlab=XLAB, ylab="Total copy ratio", xlim=genomic.limits, ylim=intensity.lim, main=paste("Segs", i, "and", i+1))
      return()
   }

   seg12$intensity = pmin(intensity.lim[2], pmax(intensity.lim[1], seg12$intensity))

   main = ifelse(!is.null(res[["seg.dat"]][["seg.info"]][i,"GC.Content"]),
         paste("Capseg Chr:", seg.chr, "Segments: ", i, "and", i+1,
               "\n Seg1 GC: ", round(res[["seg.dat"]][["seg.info"]][i,"GC.Content"], 2),
               "Seg2 GC: ", round(res[["seg.dat"]][["seg.info"]][i+1,"GC.Content"], 2)),
         paste("Capseg Chr:", seg.chr, "Segments: ", i, "and", i+1))

   plot(seg12$position, seg12$intensity, pch=19, col=seg12$col, cex=.4, main=main, xlab=XLAB, ylab="Total copy ratio", ylim = intensity.lim, xlim=genomic.limits)

   atten.tau1 = AffyAtten(delta.tau[i, 2], Theta[["at.capseg"]])
   atten.tau2 = AffyAtten(delta.tau[i+1, 2], Theta[["at.capseg"]])
   lines(range(pos1), rep(atten.tau1, 2), col=seg1.col[2], lwd=1, lty=2)
   lines(range(pos2), rep(atten.tau2, 2), col=seg2.col[2], lwd=1, lty=2)

   if (draw.legend) {
      legend("topright", legend=c("SNP-derived Mu3", "Capseg Mean"), lty=c(1, 4), lwd=3)
   }
}



PlotCapsegSegfit <- function(i, res, d.col, min=NULL, max=NULL, plot=T)
{
   h.d = res[['as.res']][["h.seg.dat"]][["h.capseg.d"]][[i]]
   tau = res[['em.fit']][["delta.tau"]][i, 2]
   Theta = res[['em.fit']][["Theta"]]
#   array.name = basename(RESULTS.DIR)
   array.name = ""

   if (length(h.d) == 0) {
      if (plot) {
         plot(1, type="n")
      }
      return(h.d)
   }

   if (is.null(min) & is.null(max)) {
      quant = quantile(h.d, probs=seq(from=.01, 1, length.out=100))
      intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))

      if (intensity.lim[2] - intensity.lim[1] < 5 ) intensity.lim = c(mean(intensity.lim) - 2.5, mean(intensity.lim) + 2.5)

      min = intensity.lim[1]
      max = intensity.lim[2]
   }

   XLAB="Total copy ratio"
   if (plot) {
      if (length(h.d) > 0) {
         df <- pmax(min, pmin(max, h.d))
         hist(df, breaks=seq(from=min, to=max, length.out=100), freq=FALSE, main="", xlab=XLAB, xlim=c(min, max), border=d.col, col=d.col)
         atten.tau = AffyAtten(tau, Theta[["at.capseg"]])
         abline( v=atten.tau, lty=2 )
      } else {
         plot(1, type="n", main=array.name, xlab=XLAB, xlim=c(min, max))
            }

      if(!is.null(tau)) {
         xgl = 1001
         x <- seq(min, max, length.out=xgl)
         y = ExomeDFunc(x, tau, Theta)
         lines(x, y, lwd=2, col="green")
      }
   } else {
      return(df)
   }
}



PlotFFit <- function( i, res, conf=.95, plot=FALSE )
{
   # i = 3
   # res = tot.res

   h.seg.dat = res[["as.res"]][["h.seg.dat"]]
   wes.f = res[["em.fit"]][["wes.f"]][i,]
   Theta = res[['em.fit']][["Theta"]]
   d = h.seg.dat[["gh.wes.allele.d"]][[i]]
   het.phase.prob = res[["em.fit"]][["het.phase.prob"]][[i]]

   f.hat = wes.f["f.hat"]
   f.H0.p = wes.f["p.H0"]
   f.H1.p = wes.f["p.H1"]

   post.col = "blue"
   f.hat.col = "orange"
   seg.col = "purple"
   conf.col = "darkgrey"
   f.loc.col = "black"

   if (length(d) == 0 )
   {
#      array.name = basename(RESULTS.DIR)
      array.name = ""
      plot(0, type="n", main=paste("No Data for seg", i, "\n", array.name), cex.main=.75)
      return()
   }

   alt = d["alt",]
   ref = d["ref",]
#   alt.phase.prob = CapturePhaseProb( alt, ref, f.hat, Theta)[,1]
   e.alpha.seg = 1 + sum(alt * het.phase.prob) + sum(ref * (1 - het.phase.prob))
   e.beta.seg = 1 + sum(ref * het.phase.prob) + sum(alt * (1 - het.phase.prob))
   e.mode = (e.alpha.seg-1)/(e.alpha.seg + e.beta.seg - 2) ## mode loc
   d.mode = dbeta( e.mode, e.alpha.seg, e.beta.seg) ## density at mode

   seg.crds = paste(paste(round(range(h.seg.dat[["gh.wes.allele.annot"]][[i]][['pos']]) / 1e6, 2), collapse="-"), "MB")
   if (!is.finite(d.mode)) {
       print("D MODE not finite")
       return()
   }
   curve( dbeta(x, e.alpha.seg, e.beta.seg) , from=0, to=1, n=1001, add=FALSE, col= post.col, xlab=paste("Fraction of alternate reads \n",  seg.crds), ylab="Density", ylim=c(0,d.mode))

   lb = qbeta((1 - conf) / 2, e.alpha.seg, e.beta.seg)
   ub = qbeta((1 + conf) / 2, e.alpha.seg, e.beta.seg)

   abline(v=c(lb, ub), col=conf.col)
   abline(v=c(f.hat, 1-f.hat), col=f.hat.col, lwd=3, lty=3)

   A1 = Theta[["f_skew"]] * f.hat
   abline(v=A1, col=f.loc.col, lwd=3, lty=3)
#   abline(v= Theta[["f_skew"]] * (1-f.hat), col=f.loc.col, lwd=3, lty=3)
   A2 =  1 - (Theta[["f_skew"]]^-1 * f.hat)
   abline(v=A2, col=f.loc.col, lwd=3, lty=3)

   med.cov = median(colSums(d))

   x = seq(0,1,length=1001)
   running_sum = rep(0, length(x))
   for( k in 1:ncol(d) )
   {
      running_sum = running_sum + dbeta( x, d["alt",k]+1, d["ref",k]+1 ) / ncol(d)
   }
  ## scale two plots to fit on same axes..
   running_sum = running_sum * (d.mode/max(running_sum) * 1/3)
   lines( x, running_sum, col=seg.col )

   legend(cex = .7, "topright",
      legend=c("F Posterior of seg", paste("F.hat =", round(f.hat, 4)),
         "Summed het posteriors", paste(round(conf*100), '% interval'),
         "f_skew * f_hat"),
   fill=c(post.col, f.hat.col, seg.col, conf.col, f.loc.col), bty="n" )
   gc.string <- ifelse(!is.null(res[["seg.dat"]][["seg.info"]][i, "GC.Content"]), paste("GC Content: ", round(res[["seg.dat"]][["seg.info"]][i, "GC.Content"], 4), "\n"), "")
   title(paste("Seg:", i, "\n", ncol(d), "het snps \n", med.cov, "median coverage \n", gc.string,
         paste("H0 Prob:", round(f.H0.p, 5), "\nH1 Prob:", round(f.H1.p, 5))), cex.main = .75 )
}



Plot_het_AF_vs_cov = function(i, res)
{
   d = h.seg.dat[["gh.wes.allele.d"]][[i]]
   alt=d["alt",]
   ref= d["ref",]

   cov = alt+ref
   AF = alt / cov

   mn =  paste( ncol(d), "het snps", sep="" )
   plot( cov, AF, pch=16, cex=0.2, main=mn, xlab="Coverage", ylab="Fraction alt reads" )

}



Plot_allelic_fraction_vs_genome = function( i, res, genomic.limits )
{
   draw_fit_lines = function( pos, f.hat, Theta )
   {
      A1 = Theta[["f_skew"]]*f.hat
#      A2 = Theta[["f_skew"]]*(1-f.hat)
      A2 = 1 - (Theta[["f_skew"]]^-1 * f.hat)

      lines(range(pos), rep(A1, 2), lwd=1, lty=2)
      lines(range(pos), rep(A2, 2), lwd=1, lty=2)
   }

   h.seg.dat = res[["as.res"]][["h.seg.dat"]]
   Theta = res[['em.fit']][["Theta"]]
   d1 = h.seg.dat[["gh.wes.allele.d"]][[i]]
   d2 = h.seg.dat[["gh.wes.allele.d"]][[i+1]]
   d = cbind(d1,d2)
   seg.chr = res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][["chr"]]

   f.hat.1 = res[["em.fit"]][["wes.f"]][i,"f.hat"]
   f.hat.2 = res[["em.fit"]][["wes.f"]][i+1,"f.hat"]

   if (length(d1) == 0 )
   {
      array.name = ""
      plot(0, type="n", main=paste("No Data for seg", i, "\n", array.name), cex.main=.75)
      return()
   }

   if (length(d2) > 0 )
   {
      het.mat.1 = res[["em.fit"]][["het.phase.prob"]][[i]]
      het.mat.2 = res[["em.fit"]][["het.phase.prob"]][[i+1]]
      het.mat = rbind(het.mat.1, het.mat.2)
   }
   else
   {
      het.mat = res[["em.fit"]][["het.phase.prob"]][[i]]
   }

   het.out = het.mat[,3]
#   het.marg = het.mat[,c(1,2), drop=FALSE] * matrix( (1 / het.out), nrow=nrow(het.mat), ncol=2, byrow=FALSE )
#   het.phase.prob = het.marg[,1]

   het.phase.prob = het.mat[,1] / (1 - het.out) # renormalize to marg out outlier state


   crds1 = h.seg.dat[["gh.wes.allele.annot"]][[i]][['pos']] / 1e6
   crds2 = h.seg.dat[["gh.wes.allele.annot"]][[i+1]][['pos']] / 1e6
   het.crds = c(crds1, crds2)

   pos1 = res[["as.res"]][['h.seg.dat']][["h.capseg.annot"]][[i]][["pos"]] / 1e6
   pos2 = res[["as.res"]][['h.seg.dat']][["h.capseg.annot"]][[i+1]][["pos"]] / 1e6
   seg.crds = c(pos1, pos2)
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
   XLAB=paste("Chromosome ", seg.chr, " position (MB)", sep="")
   plot( 0, type="n", xlab=XLAB, ylab="Fraction of alternate reads", main="", ylim=c(0,1), xlim=range(genomic.limits) )

   het_CI_low = qbeta( 0.025, alt+1, ref+1 )
   het_CI_high = qbeta( 0.975, alt+1, ref+1 )

   for( i in 1:length(het.crds))
   {
      if( !out.ix[i] ) {next}
      lines( x=rep(het.crds[i],2), y=c(het_CI_low[i],het_CI_high[i]), col=het.col[i], lwd=0.5)
   }

   points( het.crds[!out.ix], (alt/cov)[!out.ix], col=het.col[!out.ix], pch=use.pch[!out.ix], cex=0.75 )

   for( i in 1:length(het.crds))
   {
      if( out.ix[i] ) {next}
      lines( x=rep(het.crds[i],2), y=c(het_CI_low[i],het_CI_high[i]), col=het.col[i], lwd=0.5)
   }

   points( het.crds[out.ix], (alt/cov)[out.ix], col=het.col[out.ix], pch=use.pch[out.ix], cex=0.75 )

   draw_fit_lines( pos1, f.hat.1, Theta )
   draw_fit_lines( pos2, f.hat.2, Theta )


if( FALSE )
{
  ## scale two plots to fit on same axes..
   scale = 1/max(cov)
   lines( seg.crds, cov*scale )
#   ax.ix = seq( 0, max(cov), length
   axis( at=cov*scale, labels=cov, label="cov", side=4 )
   mtext( text="Coverage", side=4, las=1 )
}

}


