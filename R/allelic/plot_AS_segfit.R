## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.





########################################
## Bryan's
########################################

## Plot SNP Model Segfit
PlotSnpSegfit <- function(i, res, min=NULL, max=NULL,  plot=T) {

# i = 1
# res = iams.res
# plot=T
   
   h.snp.d = res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i]]
   h.snp.gt.p = res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]][[i]]
   h.probe.annot = res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]]
   delta.tau = res[["em.fit"]][["delta.tau"]][i, ]
   mu = GetMeans(delta.tau[1], delta.tau[2])
   theta = res[["em.fit"]][["theta"]]
   array.name = basename(results.dir)
   
   if (ncol(h.snp.d) == 0) {
      if (plot) {
         plot(1, type="n", main = array.name)
      }
      return(h.snp.d)
   }
   
   if (is.null(min) & is.null(max)) {
      quant = quantile(colSums(h.snp.d), probs=seq(from=.01, 1, length.out=100))
      intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
      
      if (intensity.lim[2] - intensity.lim[1] < 5 ) intensity.lim = c(mean(intensity.lim) - 2.5, mean(intensity.lim) + 2.5)
      min = intensity.lim[1]
      max = intensity.lim[2]
   }
   
   het.col <- "green"
   hom.col <- "black"
   
   
   
   if (is.null(mu) ) {# Simple plot
      df.snp = h.snp.d
      df.snp[1,] <- pmin(max, pmax(min, df.snp[1,]))
      df.snp[2,] <- pmin(max, pmax(min, df.snp[2,]))
      if (plot) {
         hist(df.snp, breaks=seq(from=min, to=max, length.out=100), freq=FALSE, main="", xlab="Allelic copy-ratio", xlim=c(min, max))
         return(NULL)
      } else {
         return(df.snp)   
      }
   }
   
   d.snp = h.snp.d - theta[["bg"]]
   snp.gt.p = h.snp.gt.p
   atten.mu = Atten(mu, theta[["at"]])
   
   
   sigma.epsilon <- theta[["sigma.epsilon"]]
   sigma.eta <- theta[["sigma.eta"]]
   nu <- theta[["nu"]]
   het.cov <- theta[["het.cov"]]
   het.cov <- min(het.cov, (het.cov * atten.mu[1] * atten.mu[2]))
   mu0 <- 0
   sigma.h <- GetSigmaH(sigma.epsilon, sigma.eta)
   
   ct <- colSums(snp.gt.p)
   ct <- ct / sum(ct)
   ngt.t <- c((ct[1] + ct[3]),  (ct[2] + ct[4]))
   
   
   
# df.snp <- d.snp[((d.snp >= min) & (d.snp < max))]
   df.snp = d.snp
   df.snp[1,] <- pmin(max, pmax(min, d.snp[1,]))
   df.snp[2,] <- pmin(max, pmax(min, d.snp[2,]))
   
   if (plot) {
      hist(df.snp, breaks=seq(from=min, to=max, length.out=100), freq=FALSE, main="", xlab="Allelic copy-ratio", xlim=c(min, max))
#    hist(HTx(df.snp, sigma.epsilon, sigma.eta), breaks=100, freq=FALSE, main="", xlab="Allelic copy-ratio", xlim=c(min, max))
   }
   
   xgl <- 1001
   xg <- array(NA, dim=c(4, xgl))
   x <- seq(min(df.snp), max(df.snp), length.out=xgl)
   
   ## transform grid/ params
   x.tx <- HTx(x, sigma.epsilon, sigma.eta) 
   atten.mu.tx <- HTx(atten.mu, sigma.epsilon, sigma.eta)
   mu0.tx <- HTx(mu0, sigma.epsilon, sigma.eta)
   
   xg[1, ] <- ngt.t[1] * DFunc(x.tx, mu0.tx, sigma.h, nu) * HTxVar(x, sigma.epsilon, sigma.eta) 
   xg[2, ] <- ngt.t[2] * DFunc(x.tx, atten.mu.tx[1], sigma.h, nu) * HTxVar(x, sigma.epsilon, sigma.eta)
   xg[3, ] <- ngt.t[2] * DFunc(x.tx, atten.mu.tx[2], sigma.h, nu) * HTxVar(x, sigma.epsilon, sigma.eta)
   xg[4, ] <- ngt.t[1] * DFunc(x.tx, atten.mu.tx[3], sigma.h, nu) * HTxVar(x, sigma.epsilon, sigma.eta)
   
   csum <- apply(xg, 2, sum)
   
   scale <- 1 / 2
   
   if (plot) {
      for (i in 1:4) {
         lines(x, scale * xg[i, ], lwd=2, col=hom.col)
      }
      lines(x, scale * csum, lwd=2, col="coral")
      return(NULL)
   } else {
      return (df.snp)
   }
   
}

PlotCnSegfit <- function(i, res, min=NULL, max=NULL, plot=T) {
   
# i = 30
# res = iams.res
# plot=T

   
   h.cn.d = res[["as.res"]][["h.seg.dat"]][["h.cn.d"]][[i]]
   array.name = basename(results.dir)
   
   if (length(h.cn.d) == 0) {
      if (plot) {
         plot(1, type="n")
      }
      return(h.cn.d)
   }
   
   if (is.null(min) & is.null(max)) {
      quant = quantile(h.cn.d, probs=seq(from=.01, 1, length.out=100))
      intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
      if (intensity.lim[2] - intensity.lim[1] < 5 ) intensity.lim = c(mean(intensity.lim) - 2.5, mean(intensity.lim) + 2.5)
      
      min = intensity.lim[1]
      max = intensity.lim[2]
   }
   
   theta = res[["em.fit"]][["theta"]]
   tau = res[["em.fit"]][["delta.tau"]][i, 2]
   
   if (is.null(tau)) {
      df.cn <- pmax(min, pmin(max, h.cn.d))
      if (plot) {
         hist(df.cn, breaks=seq(from=min, to=max, length.out=100), freq=FALSE, main="", xlab="CN Probe Intensity", xlim=c(min, max))
         return(NULL)
      } else {
         return(df.cn)
      }
   }
   
   d.cn = h.cn.d[1,] - theta[["cn.bg"]]
   
   nu = theta[["nu"]]
   sigma.epsilon <- theta[["sigma.epsilon"]]
   sigma.eta <- theta[["sigma.eta"]] * theta[["sigma.eta.scale"]]
   
   sigma.h <- GetSigmaH(sigma.epsilon, sigma.eta)
   tau.tx = HTx(tau, sigma.epsilon, sigma.eta)
   
   
   df.cn <- pmax(min, pmin(max, d.cn))
   
   if (plot) {
#    hist(HTx(df.cn, sigma.epsilon, sigma.eta), breaks=100, freq=FALSE, main="", xlab="CN Probe Intensity", xlim=c(min, max))
      hist(df.cn, breaks=seq(from=min, to=max, length.out=100), freq=FALSE, main="", xlab="CN Probe Intensity", xlim=c(min, max))
   }
   
   xgl = 1001
   x <- seq(min(d.cn), max(d.cn), length.out=xgl)
   x.tx = HTx(x, sigma.epsilon, sigma.eta)
   y = DFunc(x.tx, tau.tx, sigma.h, nu) * HTxVar(x, sigma.epsilon, sigma.eta)
   
   if(plot) {
      lines(x, y, lwd=2, col="green")
      return(NULL)
   } else {
      return (df.cn)
   }
}


PlotSnpIntensities <- function(i, res, genomic.limits, intensity.lim, bgl.gt.p=NULL, subtitle="")   {

# res = iams.res
# i = 9
# bgl.gt.p = NA
# subtitle="test"
# genomic.limits = range(c(res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i+1]][['pos']],
#             res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i+1]][['pos']],
#             res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i+1]][['pos']]))
# genomic.limits = genomic.limits / 1e6 
# ScarterPar()
# probe.type = "snp"
# all.d = list()
# all.d[["snp"]] = c(colSums(iams.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i]]), colSums(iams.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i+1]]))
# all.d[["cn"]] = c(iams.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]][[i]], iams.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]][[i+1]])
# all.d[["capseg"]] = c(iams.res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i]], iams.res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i+1]])
# medians = lapply(all.d, median )
# if (!is.na(medians[[probe.type]]) ) {
#    other.probe.types = setdiff(names(all.d), probe.type)
#    d1 = abs(medians[[probe.type]] - medians[[other.probe.types[1] ]])
#    d2 = abs(medians[[probe.type]] - medians[[other.probe.types[2] ]])
#    if ( max(c(d1, d2)) > 30 | is.na(max(c(d1, d2))) ) {
#       quant = quantile(all.d[[probe.type]], probs=seq(from=.01, 1, length.out=100))
#       intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
#    } else {
#       quant = quantile(c(all.d[["snp"]], all.d[["cn"]], all.d[["capseg"]]), probs=seq(from=.01, 1, length.out=100))
#       intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
#    }
#    if (intensity.lim[2] - intensity.lim[1] < 5 ) intensity.lim = c(mean(intensity.lim) - 2.5, mean(intensity.lim) + 2.5)
# }  else {
#    plot(1, type="n", xlab="Genomic Position", ylab="Copy Ratio", xlim=genomic.limits, main="")
#    return()
#    
# }
   
   ######
   seg.chr = res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]][["chr"]]
   theta = res[["em.fit"]][["theta"]]
   h.tau = res[["em.fit"]][["e.mu"]]
   delta.tau = res[["em.fit"]][['delta.tau']]
   unatten1.mu = GetMeans(delta.tau[i, 1], delta.tau[i, 2])
   unatten2.mu = GetMeans(delta.tau[i+1, 1], delta.tau[i+1, 2])
   
   if (is.null(bgl.gt.p)) {
      clust.p = rbind( res[["em.fit"]][["snp.clust.p"]][[i]], res[["em.fit"]][["snp.clust.p"]][[i+1]]) 
   } else {
      clust.p = rbind( bgl.gt.p[[i]], bgl.gt.p[[i+1]])
   }
   
   
   seg1.d = PlotSnpSegfit(i, res, min = intensity.lim[1], max=intensity.lim[2], plot=F)
   seg2.d = PlotSnpSegfit(i+1, res, min = intensity.lim[1], max=intensity.lim[2], plot=F)
   ints = cbind(seg1.d, seg2.d)
   pos1 = res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]][["pos"]] / 1e6
   pos2 = res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i+1]][["pos"]] /1e6
   pos = c(pos1, pos2)
   seg12 = data.frame(t(ints), position=pos, stringsAsFactors=F)
   seg12 = seg12[order(seg12$position), ]
   seg12$position <- seg12$position 
   
   if (nrow(seg12) == 0) {
      plot(1, type="n", xlab="Genomic Position", ylab="Copy Ratio", xlim=genomic.limits, ylim=intensity.lim, main="No Data")
      return()
   }
   seg12[,1] = pmin(intensity.lim[2], pmax(intensity.lim[1], seg12[,1]))
   seg12[,2] = pmin(intensity.lim[2], pmax(intensity.lim[1], seg12[,2]))
   
   atten1.mu = Atten(unatten1.mu, theta[['at']])
   atten2.mu = Atten(unatten2.mu, theta[["at"]])
   
   main = ifelse(!is.null(res[["seg.dat"]][["seg.info"]][i,"GC.Content"]), 
            paste("SNP Chr:", seg.chr, "Segments: ", i, "and", i+1, 
                  "\n Seg1 GC: ", round(res[["seg.dat"]][["seg.info"]][i,"GC.Content"], 2), 
                  "Seg2 GC: ", round(res[["seg.dat"]][["seg.info"]][i+1,"GC.Content"], 2)),      
            paste("SNP Chr:", seg.chr, "Segments: ", i, "and", i+1, "\n", array.name))
   
   plot(0, type="n", main=main, xlab="Genomic Position (MB)", ylab="Copy Ratio", ylim = intensity.lim, xlim=genomic.limits)
   mtext(subtitle, cex=.7)

   
   pal <- GetHapsegMarkerPal()
   cols <- pal(1000)
   state.mat <- matrix(c(0,1,2,3), ncol=4, nrow=nrow(clust.p), byrow=TRUE)
   # I have no idea how these two lines work
   col.ix.1 <- round(rowSums(state.mat * clust.p[, c(3, 2, 4, 1)]) / 3 * 999, 0) + 1 
   col.ix.2 <- round(rowSums(state.mat * clust.p[, c(1, 4, 2, 3)]) / 3 * 999, 0) + 1 
   pcols.1 <- cols[col.ix.1]
   pcols.2 <- cols[col.ix.2]
   
   max.c <- apply(clust.p, 1, which.max)
   hom.ix <- which((max.c == 1) | (max.c == 3))
   pch = 19; cex = .4
   points(seg12$position[hom.ix], seg12[hom.ix, 1], col=pcols.1[hom.ix], pch=pch, cex=cex)
   points(seg12$position[hom.ix], seg12[hom.ix, 2], col=pcols.2[hom.ix], pch=pch, cex=cex)
   
   ## Het SNPs
   if (sum(max.c == 2) > 0) {
      het.ix <- which((max.c == 2) | (max.c == 4))
      
      points(seg12$position[het.ix], seg12[het.ix, 1], col=pcols.1[het.ix], pch=pch, cex=cex)
      points(seg12$position[het.ix], seg12[het.ix, 2], col=pcols.2[het.ix], pch=pch, cex=cex)
   }

   het.col = "green"
   hom.col = "black"
   seg.bp.col = "black"
   #seg 1
# abline(v=max(pos1), col=seg.bp.col, lty=3)
   
   lines(rep(max(c(max(pos1), min(pos2))), 2), range(intensity.lim), col=seg.bp.col, lty=3)
   lines(range(pos1), rep(atten1.mu[1], 2), col = het.col)
   lines(range(pos1), rep(atten1.mu[2], 2), col = het.col)
   lines(range(pos1), rep(atten1.mu[3], 2), col = hom.col)
   lines(range(pos1), rep(0, 2), col = hom.col)
   #seg 2
   lines(range(pos2), rep(atten2.mu[1], 2), col = het.col)
   lines(range(pos2), rep(atten2.mu[2], 2), col = het.col)
   lines(range(pos2), rep(atten2.mu[3], 2), col = hom.col)
   lines(range(pos2), rep(0, 2), col = hom.col)
   
   
}

PlotCnIntensities <- function(i, res, genomic.limits, intensity.lim) {

# res = iams.res
# i = 1
# genomic.limits = range(c(res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i+1]][['pos']],
#             res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i+1]][['pos']],
#             res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i+1]][['pos']]))
# genomic.limits = genomic.limits / 1e6 
# ScarterPar()
# probe.type = "cn"
# all.d = list()
# all.d[["snp"]] = c(colSums(iams.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i]]), colSums(iams.res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i+1]]))
# all.d[["cn"]] = c(iams.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]][[i]], iams.res[["as.res"]][["h.seg.dat"]][["h.cn.d"]][[i+1]])
# all.d[["capseg"]] = c(iams.res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i]], iams.res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i+1]])
# medians = lapply(all.d, median )
# if (!is.na(medians[[probe.type]]) ) {
#    other.probe.types = setdiff(names(all.d), probe.type)
#    d1 = abs(medians[[probe.type]] - medians[[other.probe.types[1] ]])
#    d2 = abs(medians[[probe.type]] - medians[[other.probe.types[2] ]])
#    if ( max(c(d1, d2)) > 30 | is.na(max(c(d1, d2))) ) {
#       quant = quantile(all.d[[probe.type]], probs=seq(from=.01, 1, length.out=100))
#       intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
#    } else {
#       quant = quantile(c(all.d[["snp"]], all.d[["cn"]], all.d[["capseg"]]), probs=seq(from=.01, 1, length.out=100))
#       intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
#    }
#    if (intensity.lim[2] - intensity.lim[1] < 5 ) intensity.lim = c(mean(intensity.lim) - 2.5, mean(intensity.lim) + 2.5)
# }  else {
#    plot(1, type="n", xlab="Genomic Position", ylab="Copy Ratio", xlim=genomic.limits, main="")
#    return()
#    
# }
   
   
   ########
   seg.chr = res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]][["chr"]]
   theta = res[["em.fit"]][["theta"]]
   delta.tau = res[["em.fit"]][["delta.tau"]]
   tau1 = delta.tau[i, 2]
   tau2 = delta.tau[i+1, 2]
   
   seg1.col = c("red", "darkred")
   seg2.col = c("blue", "lightblue")
   
   seg1.d = PlotCnSegfit(i, res, min = intensity.lim[1], max=intensity.lim[2], plot=F) 
   seg2.d = PlotCnSegfit(i+1, res, min = intensity.lim[1], max=intensity.lim[2], plot=F)
   ints = c(seg1.d, seg2.d)
   pos1 = res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i]][["pos"]] / 1e6
   pos2 = res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i+1]][["pos"]] / 1e6
   pos = c(pos1, pos2)
   seg12 = data.frame(intensity=ints, position=pos, col=c(rep(seg1.col[1], length(seg1.d)),  rep(seg2.col[1], length(seg2.d))), stringsAsFactors=F)
   seg12 = seg12[order(seg12$position), ]
   
   if (nrow(seg12) == 0) {
      plot(1, type="n", xlab="Genomic Position (MB)", ylab="Copy Ratio", xlim=genomic.limits, ylim=intensity.lim, main=paste("Seg", i, "and", i+1))
      return()
   }
   seg12$intensity = pmin(intensity.lim[2], pmax(intensity.lim[1], seg12$intensity))
   
   main = ifelse(!is.null(res[["seg.dat"]][["seg.info"]][i,"GC.Content"]), 
         paste("CN Chr:", seg.chr, "Segments: ", i, "and", i+1, 
               "\n Seg1 GC: ", round(res[["seg.dat"]][["seg.info"]][i,"GC.Content"], 2), 
               "Seg2 GC: ", round(res[["seg.dat"]][["seg.info"]][i+1,"GC.Content"], 2)),      
         paste("CN Chr:", seg.chr, "Segments: ", i, "and", i+1, "\n", array.name))

   plot(seg12$position, seg12$intensity, pch=19, col=seg12$col, cex=.4, main=main, xlab="Genomic Position (MB)", 
         ylab="Copy Ratio", ylim = intensity.lim, xlim=genomic.limits)
   
   lines(range(pos1), rep(tau1, 2), lwd=4, col=seg1.col[2] )
   lines(range(pos2), rep(tau2, 2), lwd=4, col=seg2.col[2] )
   legend("topright", legend=c("Seg1 Tau", "Seg2 Tau"), fill=c(seg1.col[2], seg2.col[2]))
}

PlotIntensityPlus <- function(i, res) {
# res = iams.res
# i = 9
   
   genomic.limits = range(c(res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i+1]][['pos']],
               res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.cn.annot"]][[i+1]][['pos']],
               res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i]][['pos']], res[["as.res"]][["h.seg.dat"]][["h.capseg.annot"]][[i+1]][['pos']]))
   genomic.limits = genomic.limits / 1e6
   
   
   CalcIntensityLims <- function(probe.type) {
      all.d = list()
      all.d[["snp"]] = c(res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i]][1,],
            res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i]][2,],
            res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i+1]][1,],
            res[["as.res"]][["h.seg.dat"]][["h.snp.d"]][[i+1]][2,])
      all.d[["cn"]] = c(res[["as.res"]][["h.seg.dat"]][["h.cn.d"]][[i]], res[["as.res"]][["h.seg.dat"]][["h.cn.d"]][[i+1]])
      all.d[["capseg"]] = c(res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i]], res[["as.res"]][["h.seg.dat"]][["h.capseg.d"]][[i+1]])
      medians = lapply(all.d, median )
      
      if (!is.na(medians[[probe.type]])) {
         other.probe.types = setdiff(names(all.d), probe.type)
         d1 = abs(medians[[probe.type]] - medians[[other.probe.types[1] ]])
         d2 = abs(medians[[probe.type]] - medians[[other.probe.types[2] ]])
         if ( max(c(d1, d2)) > 30 | is.na(max(c(d1, d2))) ) {
            quant = quantile(all.d[[probe.type]], probs=seq(from=.01, 1, length.out=100))
            intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
         } else {
            quant = quantile(c(all.d[["snp"]], all.d[["cn"]], all.d[["capseg"]]), probs=seq(from=.01, 1, length.out=100))
            intensity.lim = c(floor(quant[[10]]), ceiling(quant[[90]]))
         }
         if (intensity.lim[2] - intensity.lim[1] < 5 ) intensity.lim = c(mean(intensity.lim) - 2.5, mean(intensity.lim) + 2.5)
         return(intensity.lim)   
      } else {
         return(NULL)
      }
      
      
   }
   
   snp.ilim = CalcIntensityLims("snp")
   cn.ilim = CalcIntensityLims("cn")
   capseg.ilim = CalcIntensityLims("capseg")
   
   ScarterPar()
   par(mfrow=c(2, 4))
   PlotSnpIntensities(i, res, genomic.limits, intensity.lim=snp.ilim, BGL.GT.P, subtitle="Beagle Phased")
   PlotSnpIntensities(i, res, genomic.limits, intensity.lim=snp.ilim, subtitle="Hapseg Phased")
   PlotCnIntensities(i, res, genomic.limits, intensity.lim=cn.ilim)
   PlotCapIntensities(i, res, genomic.limits, intensity.lim=capseg.ilim)
   
   
   Plot2DAllele(i, res, main = "Seg 1")
   Plot2DAllele(i+1, res, main = "Seg 2")
   PlotFFit(i, res, plot=T)
   PlotFFit(i+1, res, plot=T)
   
}


PlotCdf <- function(seg1.sum, seg2.sum, mrg.seg.sum, probe.type, H0.ev, H1.ev) {
   plot(ecdf(seg1.sum), xlim=xlim, ylim=c(0,1), main="", xlab="", col=h0.cols, do.points=F)
   par(new=T)
   plot(ecdf(seg2.sum), xlim=xlim, ylim=c(0,1), main="", xlab="", col=h0.cols, do.points=F)
   par(new=T)
   plot(ecdf(mrg.seg.sum), xlim=xlim, ylim=c(0,1), main=paste(probe.type, "eCDF"), xlab=ifelse(probe.type == "cn", "Copy Ratio", "A + B"), col=h1.cols, do.points=F, lwd=2)
   legend("topleft", legend=c(paste("H0:", H0.ev), paste("H1: ", H1.ev)), fill=c(h0.cols, h1.cols) )
   
}

PlotEllipseScatter <- function(x, y, a, b, add=F, ...) {
   
   idx = complete.cases(a) & complete.cases(b) & complete.cases(x) & complete.cases(y)
   a = a[idx]
   b = b[idx]
   x=x[idx] 
   y=y[idx]
   
   CalcEllipseCol <- function(a, b){
      e.area = pi*a*b
      scale = quantile(e.area, 1:100/100)[90] / tan(.8)
      alpha = atan((e.area) / scale) / (pi / 2)
      grey.scale = rgb(alpha, alpha, alpha, (1-alpha)/2)
      return(grey.scale)
   }
   if (!add) {
      plot(0, type="n", ...)   
   }
   col = CalcEllipseCol(a, b)   
   min.rad = .025
   draw.ellipse(x=x, y=y, a=pmax(min.rad, a), b=pmax(min.rad, b), angle = 0, segment = c(0, 360), 
         arc.only = TRUE, deg = TRUE, nv = 100, border = NA, col = col, lty = 1, lwd = 1)
   
   pars = list(...)
   lines(c(0, max(pars[['xlim']])), c(0, max(pars[['ylim']])), lty=3, col="black", lwd=2)
   lines(c(0, max(pars[['ylim']])), c(0,0), lty=3, col="grey")
   lines(c(0,0), c(0, max(pars[['ylim']])), lty=3, col="grey")
   
}

PlotBivariateGaussian <- function(x.mu, y.mu, x.sigma, y.sigma, add=F, ...) {
   
   idx = complete.cases( x.sigma ) & complete.cases(y.sigma) & complete.cases(x.mu) & complete.cases(y.mu)
   x.sigma = x.sigma[idx]
   y.sigma = y.sigma[idx]
   x.mu=x.mu[idx] 
   y.mu=y.mu[idx]

   # N = 10
   # x.sigma = seq(1, 3, length.out=N)
   # y.sigma = seq(1, 5, length.out=N)
   # x.mu=seq(1, 10, length.out=N)
   # y.mu=seq(-5, 5, length.out=N)
   browser()
   pars = list(...)
   n.samples = 1e2
   require(mnormt)
   require(hexbin)
   pts = foreach(i=seq_along(x.mu), .combine=rbind) %dopar% {loopStatus(i,1 )
      rmnorm(n.samples, mean=c(x.mu[i], y.mu[i]), varcov=matrix(c(x.sigma[i], 0, 0, y.sigma[i]), ncol=2))
   }
   h = hexbin(pts[,1], pts[,2], xbins = 100, shape = 1,xlab = NULL, ylab = NULL, IDs = FALSE, xbnds=c(0,5), ybnds=c(0,5))
   # plot(h, xlab="", ylab="", xlim = pars[['xlim']], ylim = pars[['ylim']])
   # plot(h, xlab="", ylab="", xlim = c(0,5), ylim = c(0,5))
   plot(h, xlab="", ylab="")

   
   
   lines(c(0, max(pars[['xlim']])), c(0, max(pars[['ylim']])), lty=3, col="grey")
   lines(c(0, max(pars[['ylim']])), c(0,0), lty=3, col="grey")
   lines(c(0,0), c(0, max(pars[['ylim']])), lty=3, col="grey")
   
}

PlotWesSnpConcordance <- function(res, conf, type="unatten") {

   # conf = .7
   # res = iams.res
   # type = 'unatten'
   
   require(plotrix)
   lim = c(-0.5, 5)
   ScarterPar()
   par(mfcol=c(2,2))
   
   all.conf = list()
   unatten.snp.mu = t(apply(res[["em.fit"]][["delta.tau"]], 1, function(x) GetMeans(x[1], x[2])))
   
   
   #####
   #  WES minor v SNP minor
   #####
   
   n.segs = length(res[["em.fit"]][["cap.e.mu"]][,"mu1"])
   wes.minor = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu1"], res[["em.fit"]][["theta"]][["at.capseg"]]), 
         res[["em.fit"]][["cap.e.mu"]][,"mu1"])
   snp.minor = ifelse(rep(type=="unatten", n.segs),
         unatten.snp.mu[,1],
         Atten(unatten.snp.mu[,1], res[["em.fit"]][["theta"]][["at"]]))
   
   wes.conf = lapply(1:n.segs, function(i) CalcConf(mu = wes.minor[i], sigma=res[["em.fit"]][["cap.e.mu"]][i,"sigma1"], conf=conf) )
   snp.conf = lapply(1:n.segs, function(i) {
            sigma = (res[["em.fit"]][["delta.tau.sd"]][i, "tau"]^2 + res[["em.fit"]][["delta.tau.sd"]][i, "delta"]^2)^(1/2) / 2
            mu = snp.minor[i]
            CalcConf(mu=mu, sigma=sigma, conf=conf)
         } )
   all.conf[["minor"]] = cbind(wes=sapply(wes.conf, function(x) abs(diff(x))), snp=sapply(snp.conf, function(x) abs(diff(x))))
   
   PlotEllipseScatter(x=snp.minor, y=wes.minor, a=all.conf[["minor"]][,"snp"], b=all.conf[["minor"]][,"wes"], xlim=lim, ylim=lim, xlab="SNP", 
         ylab="WES", main=paste(type, " minor comparisons \n", round(100*conf), "% confidence interval"))
   
   
   #####
   # WES Major v SNP Major
   #####
   
   wes.major = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu2"], res[["em.fit"]][["theta"]][["at.capseg"]]), 
         res[["em.fit"]][["cap.e.mu"]][,"mu2"])
   snp.major = ifelse(rep(type=="unatten", n.segs),
         unatten.snp.mu[,2],
         Atten(unatten.snp.mu[,2], res[["em.fit"]][["theta"]][["at"]]))
   
   wes.conf = lapply(1:n.segs, function(i) CalcConf(mu = wes.major[i], sigma=res[["em.fit"]][["cap.e.mu"]][i,"sigma2"], conf=conf) )
   snp.conf = lapply(1:n.segs, function(i) {
            sigma = (res[["em.fit"]][["delta.tau.sd"]][i, "tau"]^2 + res[["em.fit"]][["delta.tau.sd"]][i, "delta"]^2)^(1/2) / 2
            mu = snp.major[i]
            CalcConf(mu=mu, sigma=sigma, conf=conf)
         } )
   all.conf[["major"]] = cbind(wes=sapply(wes.conf, function(x) abs(diff(x))), snp=sapply(snp.conf, function(x) abs(diff(x))))
   PlotEllipseScatter(x=snp.major, y=wes.major, a=all.conf[["major"]][,"snp"], b=all.conf[["major"]][,"wes"], xlim=lim, ylim=lim, xlab="SNP", 
         ylab="WES", main=paste(type, "major comparisons \n", round(100*conf), "% confidence interval"))   
   
   #####
   #  WES Hom v SNP Hom
   #####
   
   wes.hom = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu3"], res[["em.fit"]][["theta"]][["at.capseg"]]), 
         res[["em.fit"]][["cap.e.mu"]][,"mu3"])
   snp.hom = ifelse(rep(type=="unatten", n.segs),
         unatten.snp.mu[,3],
         Atten(unatten.snp.mu[,3], res[["em.fit"]][["theta"]][["at"]]))
   
   snp.conf = lapply(1:length(snp.hom), function(i) CalcConf(mu = snp.hom[i], sigma=res[["em.fit"]][["delta.tau.sd"]][i, "tau"], conf=conf) )
   wes.conf = lapply(1:length(wes.hom), function(i) CalcConf(mu = wes.hom[i] , sigma=res[["em.fit"]][["cap.e.mu"]][i,"sigma3"], conf=conf) )
   all.conf[["hom"]] = cbind(wes=sapply(wes.conf, function(x) abs(diff(x))), snp=sapply(snp.conf, function(x) abs(diff(x))))
   
   PlotEllipseScatter(x=snp.hom, y=wes.hom, a=all.conf[["hom"]][,"snp"], b=all.conf[["hom"]][,"wes"], xlim=lim, ylim=lim, xlab="SNP", 
         ylab="WES", main=paste(type, " hom comparisons \n", round(100*conf), "% confidence interval"))
   
   #####
   #  WES.mu3 - SNP.Mu3 v. delta
   #####
   
   y = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu3"], res[["em.fit"]][["theta"]][["at.capseg"]]) - res[["em.fit"]][["delta.tau"]][,2], 
         res[["em.fit"]][["cap.e.mu"]][,"mu3"] - Atten(res[["em.fit"]][["delta.tau"]][,2], res[["em.fit"]][["theta"]][["at"]]))
   delta = res[["em.fit"]][["delta.tau"]][,1]
   delta.conf = lapply(1:length(delta), function(i) CalcConf(mu = delta[i], sigma=res[["em.fit"]][["delta.tau.sd"]][i, "delta"], conf=conf) )
   y.conf = lapply(1:length(y), function(i) {
            snp.mu3.sigma = res[["em.fit"]][["delta.tau.sd"]][i, "tau"]
            wes.mu3.sigma = res[["em.fit"]][["cap.e.mu"]][i,"sigma3"] 
            sigma = (wes.mu3.sigma^2 + snp.mu3.sigma^2)^(1/2)
            mu = y[i]
            CalcConf(mu=mu, sigma=sigma, conf=conf)
         } )
   all.conf[["hom.diff.v.delta"]] = cbind(hom.diff=sapply(y.conf, function(x) abs(diff(x))), delta=sapply(delta.conf, function(x) abs(diff(x))))
   
   PlotEllipseScatter(x=delta, y=y, a=all.conf[["hom.diff.v.delta"]][,"delta"], b=all.conf[["hom.diff.v.delta"]][,"hom.diff"], xlim=c(-1,5), ylim=c(-3,3), 
         xlab="Delta", ylab="WES Hom - SNP Hom", main=paste("Delta Dependence on Concordance \n", round(100*conf), "% confidence interval"))
   
   
   # ######
   # # Plot 5: Distributions of error bars
   # ######
   
   # all.conf
   # for (i in 1:length(all.conf)) {
   #    i = 1
   #    x1 = all.conf[[i]][,1]
   #    x2 = all.conf[[i]][,2]
   #    plot.name = names(all.conf)[i]
   #    TwoHist (x1, x2, breaks="Sturges", main=plot.name)
   # }
   
}

PlotWesSnpConcordanceMinorMajor <- function(res, conf, type="unatten") {
   
   # conf = .7
   # res = iams.res
   # type = 'unatten'
   
   require(plotrix)
   lim = c(-0.5, 5)
   ScarterPar()
   par(mfcol=c(1,1))
   
   all.conf = list()
   unatten.snp.mu = t(apply(res[["em.fit"]][["delta.tau"]], 1, function(x) GetMeans(x[1], x[2])))
   
   
   #####
   #  WES minor v SNP minor
   #####
   
   n.segs = length(res[["em.fit"]][["cap.e.mu"]][,"mu1"])
   n.wes.na = sum(is.na(res[["em.fit"]][["cap.e.mu"]][,"mu1"]))
   not.modeled.idx = which(is.na(res[["em.fit"]][["cap.e.mu"]][,"mu1"]))
   frac.not.modeled = sum(sapply(not.modeled.idx, function(i) {
      res[["seg.dat"]][["seg.info"]][i, "End.bp"] - res[["seg.dat"]][["seg.info"]][i, "Start.bp"] + 1
      })) / 3.2e9
   wes.minor = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu1"], res[["em.fit"]][["theta"]][["at.capseg"]]), 
         res[["em.fit"]][["cap.e.mu"]][,"mu1"])
   snp.minor = ifelse(rep(type=="unatten", n.segs),
         unatten.snp.mu[,1],
         Atten(unatten.snp.mu[,1], res[["em.fit"]][["theta"]][["at"]]))
   
   wes.conf = lapply(1:n.segs, function(i) CalcConf(mu = wes.minor[i], sigma=res[["em.fit"]][["cap.e.mu"]][i,"sigma1"], conf=conf) )
   snp.conf = lapply(1:n.segs, function(i) {
            sigma = (res[["em.fit"]][["delta.tau.sd"]][i, "tau"]^2 + res[["em.fit"]][["delta.tau.sd"]][i, "delta"]^2)^(1/2) / 2
            mu = snp.minor[i]
            CalcConf(mu=mu, sigma=sigma, conf=conf)
         } )
   
   PlotEllipseScatter(x=snp.minor, y=wes.minor, a=all.conf[["minor"]][,"snp"], b=all.conf[["minor"]][,"wes"], xlim=lim, ylim=lim, xlab="SNP", 
         ylab="WES", main="", 
         add = F)
   
   
   #####
   # WES Major v SNP Major
   #####
   
   wes.major = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu2"], res[["em.fit"]][["theta"]][["at.capseg"]]), 
         res[["em.fit"]][["cap.e.mu"]][,"mu2"])
   snp.major = ifelse(rep(type=="unatten", n.segs),
         unatten.snp.mu[,2],
         Atten(unatten.snp.mu[,2], res[["em.fit"]][["theta"]][["at"]]))
   
   wes.conf = lapply(1:n.segs, function(i) CalcConf(mu = wes.major[i], sigma=res[["em.fit"]][["cap.e.mu"]][i,"sigma2"], conf=conf) )
   snp.conf = lapply(1:n.segs, function(i) {
            sigma = (res[["em.fit"]][["delta.tau.sd"]][i, "tau"]^2 + res[["em.fit"]][["delta.tau.sd"]][i, "delta"]^2)^(1/2) / 2
            mu = snp.major[i]
            CalcConf(mu=mu, sigma=sigma, conf=conf)
         } )
   all.conf[["major"]] = cbind(wes=sapply(wes.conf, function(x) abs(diff(x))), snp=sapply(snp.conf, function(x) abs(diff(x))))
   PlotEllipseScatter(x=snp.major, y=wes.major, a=all.conf[["major"]][,"snp"], b=all.conf[["major"]][,"wes"], xlim=lim, ylim=lim, xlab="SNP", 
         ylab="WES", main="", 
         add=T)   
   title(paste(type, "minor and major comparisons \n", round(100*conf), "% confidence interval", n.wes.na, "segs (", signif(frac.not.modeled * 100, 2), "%) not modeled by WES platform"))
   

   
}

PlotWesSnpConcordanceMinorMajorOld <- function(res, conf, type="unatten") {
   
   # conf = .7
   # res = iams.res
   # type = 'unatten'
   
   require(plotrix)
   lim = c(-0.5, 5)
   ScarterPar()
   par(mfcol=c(1,1))
   
   all.conf = list()
   unatten.snp.mu = res[["em.fit"]][["e.mu"]]
   
   
   #####
   #  WES minor v SNP minor
   #####
   
   n.segs = length(res[["em.fit"]][["cap.e.mu"]][,"mu1"])
   n.wes.na = sum(is.na(res[["em.fit"]][["cap.e.mu"]][,"mu1"]))
   not.modeled.idx = which(is.na(res[["em.fit"]][["cap.e.mu"]][,"mu1"]))
   frac.not.modeled = sum(sapply(not.modeled.idx, function(i) {
      res[["seg.dat"]][["seg.info"]][i, "End.bp"] - res[["seg.dat"]][["seg.info"]][i, "Start.bp"] + 1
      })) / 3.2e9
   wes.minor = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu1"], res[["em.fit"]][["theta"]][["at.capseg"]]), 
         res[["em.fit"]][["cap.e.mu"]][,"mu1"])
   snp.minor = ifelse(rep(type=="unatten", n.segs),
         unatten.snp.mu[,1],
         Atten(unatten.snp.mu[,1], res[["em.fit"]][["theta"]][["at"]]))
   
   wes.conf = lapply(1:n.segs, function(i) CalcConf(mu = wes.minor[i], sigma=res[["em.fit"]][["cap.e.mu"]][i,"sigma1"], conf=conf) )
   snp.conf = lapply(1:n.segs, function(i) {
            sigma = (res[["em.fit"]][["mu.post.sd"]][i, "tau"]^2 + res[["em.fit"]][["mu.post.sd"]][i, "delta"]^2)^(1/2) / 2
            mu = snp.minor[i]
            CalcConf(mu=mu, sigma=sigma, conf=conf)
         } )
   all.conf[["minor"]] = cbind(wes=sapply(wes.conf, function(x) abs(diff(x))), snp=sapply(snp.conf, function(x) abs(diff(x))))
   
   PlotEllipseScatter(x=snp.minor, y=wes.minor, a=all.conf[["minor"]][,"snp"], b=all.conf[["minor"]][,"wes"], xlim=lim, ylim=lim, xlab="SNP", 
         ylab="WES", main="", 
         add = F)
   
   
   #####
   # WES Major v SNP Major
   #####
   
   wes.major = ifelse(rep(type=="unatten", n.segs), 
         InvAtten(res[["em.fit"]][["cap.e.mu"]][,"mu2"], res[["em.fit"]][["theta"]][["at.capseg"]]), 
         res[["em.fit"]][["cap.e.mu"]][,"mu2"])
   snp.major = ifelse(rep(type=="unatten", n.segs),
         unatten.snp.mu[,2],
         Atten(unatten.snp.mu[,2], res[["em.fit"]][["theta"]][["at"]]))
   
   wes.conf = lapply(1:n.segs, function(i) CalcConf(mu = wes.major[i], sigma=res[["em.fit"]][["cap.e.mu"]][i,"sigma2"], conf=conf) )
   snp.conf = lapply(1:n.segs, function(i) {
            sigma = (res[["em.fit"]][["mu.post.sd"]][i, "tau"]^2 + res[["em.fit"]][["mu.post.sd"]][i, "delta"]^2)^(1/2) / 2
            mu = snp.major[i]
            CalcConf(mu=mu, sigma=sigma, conf=conf)
         } )
   all.conf[["major"]] = cbind(wes=sapply(wes.conf, function(x) abs(diff(x))), snp=sapply(snp.conf, function(x) abs(diff(x))))
   PlotEllipseScatter(x=snp.major, y=wes.major, a=all.conf[["major"]][,"snp"], b=all.conf[["major"]][,"wes"], xlim=lim, ylim=lim, xlab="SNP", 
         ylab="WES", main="", add=T)      
   title(paste(type, "minor and major comparisons \n", round(100*conf), "% confidence interval \n", n.wes.na, "segs (", signif(frac.not.modeled * 100, 2), "%) not modeled by WES platform"))
}

PlotSnpAscn <- function(i, res, y.lab="Allelic copy-ratio", xcrds=NA, chr=NA,
      seg.wd=2.3, add=FALSE, vertical=FALSE,
      plot.markers=TRUE, col.markers=TRUE,
      plot.segs=TRUE, het.seg.col.1=NA, het.seg.col.2=NA,
      loci.annot=NA, merge.prob=NA, seg.log.ev=NA,
      plot.tcn.seg=TRUE, d.quantiles=c(.1, .9),
      marker.gt.pal=NA, pch=".", cex=1.2, xaxt="s") {
   
   
# y.lab="Allelic copy-ratio"
# xcrds=NA; chr=NA;
# seg.wd=1; add=FALSE; vertical=FALSE;
# plot.markers=TRUE; col.markers=TRUE;
# plot.segs=TRUE; het.seg.col.1=NA; het.seg.col.2=NA;
# loci.annot=NA; merge.prob=NA; seg.log.ev=NA;
# plot.tcn.seg=TRUE; d.quantiles=c(.1, .9);
# marker.gt.pal=NA; pch=19; cex=.5; xaxt="s"
# res = iams.res
# i = 1
   
   d = res$as.res$h.seg.dat$h.snp.d[[i]] 
   theta = res$em.fit$theta
   mu = GetMeans(res[["em.fit"]][["delta.tau"]][i,1], res[["em.fit"]][["delta.tau"]][i,2])
   e.mu = Atten(mu, theta[["at"]])
   sigma.epsilon = theta[["sigma.epsilon"]]
   sigma.eta = theta[["sigma.eta"]]
   snp.clust.p = res[["em.fit"]][["snp.clust.p"]][[i]]
   snp.annot = res$as.res$h.seg.dat$h.snp.annot[[i]]
   wes.minor.af = res[["em.fit"]][["wes.f"]][i]
   e.mu.wes = Atten(c(wes.minor.af * mu[3], (1 - wes.minor.af) * mu[3]), theta[["at"]])
   array.name = basename(results.dir)
   
   if (length(d) == 0 ) {
      plot(1, type="n", main=paste("No SNP data in this seg."))
      return(NULL)
   }
   
   
   quant = quantile(c(d[1,], d[2, ]), probs=seq(from=.01, 1, length.out=100))
   min = floor(quant[[round(100*d.quantiles[1])]])
   max = ceiling(quant[[round(100*d.quantiles[2])]])
   
   
   
   
   if (is.na(xcrds)) xcrds = snp.annot[["pos"]] / 1e6
   if (is.na(chr)) chr = snp.annot[["chr"]] 
   
   
   max.c <- apply(snp.clust.p, 1, which.max)
   
   het.af.col = "orange"
   out.color <- "white"
   het.col <- "green"
   hom.col <- "black"
   if (any(is.na(marker.gt.pal))) {
      pal <- GetHapsegMarkerPal()
   } else {
      pal <- marker.gt.pal
   }
   
   cols <- pal(1000)
   state.mat <- matrix(c(0,1,2,3), ncol=4, nrow=nrow(snp.clust.p), byrow=TRUE)
   
   col.ix.1 <- round(rowSums(state.mat * snp.clust.p[, c(3, 2, 4, 1)]) / 3 * 999, 0) + 1 
   col.ix.2 <- round(rowSums(state.mat * snp.clust.p[, c(1, 4, 2, 3)]) / 3 * 999, 0) + 1 
   
   if (col.markers) {
      pcols.1 <- cols[col.ix.1]
      pcols.2 <- cols[col.ix.2]
   } else {
      pcols.1 <- rep(1, length(col.ix.1))
      pcols.2 <- rep(1, length(col.ix.2))
   }
   
   if (is.na(het.seg.col.1)) {
      het.seg.col.1 <- het.col
   }
   if (is.na(het.seg.col.2)) {
      het.seg.col.2 <- het.col
   }
   
   hom.seg.col <- hom.col
   
   ## FIXME: double check these two &s, they look like they should be &&
   if (!add) {
      x.lab <- ifelse(xaxt=="n", "", paste("Chr", chr, "position (Mb)"))
      main = paste("Chr: ", chr, "seg:", i)
      plot(0, type="n", xlab=x.lab, ylab=y.lab, ylim=c(min, max), xlim=range(xcrds), xaxt=xaxt, main=main)
      if (any(is.na(xcrds))) {
         xcrds <- c(1:ncol(d))
         plot(0, type="n", xlab="Marker order", ylab=y.lab,
               ylim=c(min, max), xlim=range(xcrds), main=main)
      }
   }
   
   if (plot.markers) {
      ## Hom SNPs
      hom.ix <- which((max.c == 1) | (max.c == 3))
      
      points(xcrds[hom.ix], d[1, hom.ix], col=pcols.1[hom.ix], pch=pch,cex=cex)
      points(xcrds[hom.ix], d[2, hom.ix], col=pcols.2[hom.ix], pch=pch,cex=cex)
      
      ## Het SNPs
      if (sum(max.c == 2) > 0) {
         ## FIXME: double check that |, it looks like it should be ||
         het.ix <- which((max.c == 2) | (max.c == 4))
         
         points(xcrds[het.ix], d[1, het.ix], col=pcols.1[het.ix], pch=pch, cex=cex)
         points(xcrds[het.ix], d[2, het.ix], col=pcols.2[het.ix], pch=pch, cex=cex)
      }
   }
   
   wd <- 0.04
   
   if (plot.segs) {
      if (!vertical) {
         if (add) {
            rect(xleft=min(xcrds), ybottom=e.mu[1] - wd / 2, xright=max(xcrds),
                  ytop=e.mu[1] + wd / 2, border=NA, col=het.seg.col.1)
            rect(xleft=min(xcrds), ybottom=e.mu[2] - wd / 2,
                  xright=max(xcrds), ytop=e.mu[2] + wd / 2, border=NA,
                  col=het.seg.col.2)
         } else {
            lines(x=range(xcrds), y=c(e.mu[1], e.mu[1]), col=het.seg.col.1,
                  lwd=seg.wd)
            lines(x=range(xcrds), y=c(e.mu[2], e.mu[2]), col=het.seg.col.2,
                  lwd=seg.wd)
         }
         if (plot.tcn.seg) {
            lines(x=range(xcrds), y=c(e.mu[3], e.mu[3]), col=hom.seg.col,
                  lwd=seg.wd)
            lines(x=range(xcrds), y=c(0, 0), col=hom.seg.col, lwd=seg.wd)
         }
      }
   }
   
   if (! (is.na(wes.minor.af) | is.null(wes.minor.af) )) {
#    lines(x=range(xcrds), y=rep(e.mu.wes[1], 2), col=het.af.col, lwd=seg.wd)
#    lines(x=range(xcrds), y=rep(e.mu.wes[2], 2) * c(e.mu[3], e.mu[3]), col=het.af.col, lwd=seg.wd)
      
      abline(h=e.mu.wes[1], col=het.af.col, lwd=seg.wd)
      abline(h=e.mu.wes[2], col=het.af.col, lwd=seg.wd)
      
   }
   
   if (!is.na(loci.annot)) {
      n <- length(xcrds)
      msg <- paste("N = ", n, ", log_ev/N = ", round( seg.log.ev / n, 3 ),
            ", MPr = ", round(merge.prob[1], 5), ",\nlog10 MPr = ",
            round( merge.prob[2], 2), sep="" )
      mtext(msg, line=0, side=3, adj=0, cex=par("cex")*par("cex.axis"))
      ix <- (loci.annot[,"chr" ] == chr)
      
      if (sum(ix) == 0) {
         return()
      }
      
      loci.annot <- loci.annot[ix, ]
      gcrd <- loci.annot[, "pos"] / 1e6
      
      if (!vertical) {
         segments(x0=gcrd, y0=-1, x1=gcrd, y1=4, lty=2)
         text(gcrd, y=rep(4.5, length(gcrd)),
               labels=paste("Pr=", round(loci.annot[, "log10_prob"], 5),
                     sep=""), cex=0.5, srt=90)
      } else {
         lines(y=loci.annot[, "pos"] / 1e6, x=c(-1, 4), lty=2)
         text(y=gcrd, x=rep(5, length(gcrd)),
               labels=paste("Pr=", round(loci.annot[, "prob"], 5), sep=""), cex=0.5)
      }
   }
}
   

Plot2DAllele <- function(i, res, main="") {   
   # i = 163
   # res = iams.res
   
   d = res$as.res$h.seg.dat$h.snp.d[[i]] 
   snp.clust.p = res$em.fit$snp.clust.p[[i]]
   snp.gt.p = res[["as.res"]][["h.seg.dat"]][["h.snp.gt.p"]][[i]]
   theta = res$em.fit$theta
   e.mu = Atten(GetMeans(res[["em.fit"]][["delta.tau"]][i,1], res[["em.fit"]][["delta.tau"]][i,2]),
      theta[["at"]])
   


   
   sigma.epsilon <- theta[["sigma.epsilon"]]
   sigma.eta <- theta[["sigma.eta"]]
   nu <- theta[["nu"]]
   het.cov <- theta[["het.cov"]]
   het.cov <- min(het.cov, (het.cov * e.mu[1] * e.mu[2]))
   mu0 <- 0
   sigma.h <- GetSigmaH(sigma.epsilon, sigma.eta)
   
   
   x.tx <- HTx(x, sigma.epsilon, sigma.eta) 
   e.mu <- HTx(e.mu, sigma.epsilon, sigma.eta)
   mu0 <- 0
   mu0.tx <- HTx(mu0, sigma.epsilon, sigma.eta)
   
   ct <- colSums(snp.gt.p)
   ct <- ct / sum(ct)
   ngt.t <- c((ct[1] + ct[3]),  (ct[2] + ct[4]))
   
   ## 2d plots
   het.af.col = "orange"
   out.color <- "white"
   het.col <- "green"
   hom.col <- "black"
   t.min <- -0.5
   t.max <- 2.75
   cex <- 2
   min <- -0.5
   max <- 5.25
   gl <- 250
   x2d <- seq(t.min, t.max, length.out=gl)
   x1 <- rep(x2d, gl)
   x2 <- as.vector(sapply(x2d, rep, gl))
   x <- cbind(x1, x2)
   scale <- 1 / 2
   

   if (length(d) == 0) {
         plot(0, type="n", main="No Data", xlim=c(t.min, t.max), ylim=c(t.min, t.max), 
            xlab="A allele copy-ratio", ylab="B allele copy-ratio", bty="n")
         return(NULL)
      }

   sigma <- matrix(c(sigma.h^2, het.cov, het.cov, sigma.h^2), nrow=2, ncol=2)
   invert.sigma <- solve(sigma)
   hom.sigma <- matrix(c(sigma.h^2, 0, 0, sigma.h^2), nrow=2, ncol=2)
   invert.hom.sigma <- solve(hom.sigma)
   
   het.1 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[1], e.mu[2]), sigma, invert.sigma, nu), dim=c(gl, gl))
   het.2 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[2], e.mu[1]), sigma, invert.sigma, nu), dim=c(gl, gl))
   
   hom.1 <- array(scale * ngt.t[1] * DmvFunc(x, c(mu0.tx, e.mu[3]), hom.sigma, invert.hom.sigma, nu), dim=c(gl, gl))
   hom.2 <- array(scale * ngt.t[1] * DmvFunc(x, c(e.mu[3], mu0.tx), hom.sigma, invert.hom.sigma, nu), dim=c(gl, gl))
   csum2d <- het.1 + het.2 + hom.1 + hom.2
   
   ## compute densities at modes
   het.md <- DmvFunc(c(e.mu[1], e.mu[2]), c(e.mu[1], e.mu[2]), sigma, invert.sigma, nu) * scale * ngt.t[2] 
   hom.md <- DmvFunc(c(mu0.tx, e.mu[3]), c(mu0.tx, e.mu[3]), hom.sigma, invert.hom.sigma, nu) * scale * ngt.t[1]
   
   level_scales <- c(0.99, 0.95, 0.8, 0.5, 0.2)
   
   out.pch <- rep(".", ncol(d))
   pal <- GetHapsegMarkerPal()
   
   
   cols <- rev(pal(1000))
   state.mat <- matrix(c(0,1,2,3), ncol=4, nrow=nrow(snp.clust.p), byrow=TRUE)
   col.ix <- round(rowSums(state.mat * snp.clust.p[, c(3, 2, 4, 1)]) / 3 * 999, 0) + 1 
   pcols <- cols[col.ix]
   
   x <- d[1, ]
   y <- d[2, ]
   ix <- (x < max) & (x > min) & (y < max) & (y > min)
   
   
   ## tx data 2d fits
   d.tx <- HTx(d, sigma.epsilon, sigma.eta)
   plot(0, type="n", main=main, xlim=c(t.min, t.max), ylim=c(t.min, t.max), xlab="A allele copy-ratio", 
         ylab="B allele copy-ratio", bty="n")
   points(d.tx[1, ], d.tx[2, ], col=pcols, pch=out.pch, cex=cex)
   
   contour(x=x2d, y=x2d, z=het.1, add=TRUE, drawlabels=FALSE,
         levels=level_scales*het.md, col=het.col)
   contour(x=x2d, y=x2d, z=het.2, add=TRUE, drawlabels=FALSE,
         levels=level_scales*het.md, col=het.col)
   contour(x=x2d, y=x2d, z=hom.1, add=TRUE, drawlabels=FALSE,
         levels=level_scales*hom.md, col=hom.col)
   contour(x=x2d, y=x2d, z=hom.2, add=TRUE, drawlabels=FALSE,
         levels=level_scales*hom.md, col=hom.col)
   
   abline(v=0, lty=2)
   abline(h=0, lty=2)
   
}

PlotAll <- function(res, segs=NULL, plots=1:2, main="", save = F, subdir="") {

   print("Plotting...")
# res = iams.res
# segs = NULL
   if (is.null(segs)) {
      segs = 1:(length(res[["as.res"]][["h.seg.dat"]][["h.snp.d"]])-1)
   }
   segs = sapply(segs, function(i) {
            
            if (res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i]][["chr"]] != res[["as.res"]][["h.seg.dat"]][["h.snp.annot"]][[i+1]][["chr"]]) {
               return(NULL)
            } else return(i)
         })
   segs = unlist(segs)
   foreach(i = segs, .combine=c) %dopar% { loopStatus(i, step=1)
#    res = iams.res
#    i = 8
#    save = F
#    PROBE.TYPES = c("cn", "snp", "cap")
#    plots=1:3; main = ""
      
      seg1.col = c("red", "darkred")
      seg2.col = c("purple", "darkblue")
      
      if (1 %in% plots) {
         
         OpenDev(save=save, fn=file.path(results.dir, "plots", subdir, paste("Segs_", i,"_",i+1, "a.jpeg", sep="")))
         PlotAllIntensities(i, res)
         
         if (save) dev.off()
      }
      
      
      if (2 %in% plots) {
         OpenDev(save=save, fn=file.path(results.dir, "plots", subdir, paste("Segs_", i,"_",i+1, "b.jpeg", sep="")))
         
         # CN Model
         h0.cols = c("blue")
         h1.cols = c("green")
         ScarterPar()
         par(cex=4)
         par(mfrow=c(2, 3))
         xlim = c(-0.5, 5)
         ylim=c(0,1)
         
         # PDF
         PlotCnSegfit(i, res)
         PlotSnpSegfit(i, res)
         PlotCapsegSegfit(i, res)
         PlotCnSegfit(i+1, res)
         PlotSnpSegfit(i+1, res)
         PlotCapsegSegfit(i+1, res)
         title(array.name)
         if(save) dev.off()
      }      
      
      if (3 %in% plots ) {
         OpenDev(save=save, fn=file.path(results.dir, "plots", subdir, paste("Segs_", i,"_",i+1, "c.jpeg", sep="")))
         
         ScarterPar()
         par(mfrow=c(1, 2))
   
         PlotFFit(i, res, conf=.95, plot=T)
         
         PlotSnpAscn(i, res, y.lab="Allelic copy-ratio", xcrds=NA, chr=NA,
               add=FALSE, vertical=FALSE,
               plot.markers=TRUE, col.markers=TRUE,
               plot.segs=TRUE, het.seg.col.1=NA, het.seg.col.2=NA,
               loci.annot=NA, merge.prob=NA, seg.log.ev=NA,
               plot.tcn.seg=TRUE,
               marker.gt.pal=NA, pch=20, cex=.4, xaxt="s")
         
         if(save) dev.off()
      }
      
   }
   
   
}




FindClosestSegNum <- function(res, chr, pos) {
# res = iams.res
# chr = 1
# pos = 72759524
   
   segs = res$seg.dat$seg.info
   out = which(segs$Chromosome == chr & segs$Start.bp <= pos & segs$End.bp >= pos)
   if (length(out) != 1) stop (paste(length(out), "segs found matching chr:", chr, "pos:", pos))
   return(out)
}




########################################
## End Bryan's
########################################



PlotSegFit <- function(d, snp.gt.p, snp.clust.p, snp.annot, e.mu,
      theta, merged.loci=NA, merge.prob, seg.log.ev,
      marker.gt.pal=NA, verbose=FALSE) {
   
   i = 1
   d = iams.res$as.res$h.seg.dat$h.snp.d[[i]] 
   snp.gt.p = iams.res$as.res$h.seg.dat$h.snp.gt.p[[i]]
   snp.clust.p = iams.res$em.fit$snp.clust.p[[i]]
   e.mu = GetMeans(iams.res[["em.fit"]][["delta.tau"]][i,1], iams.res[["em.fit"]][["delta.tau"]][i,2])
   theta = iams.res$em.fit$theta
   y.lab = NULL
   snp.annot = iams.res$as.res$h.seg.dat$h.snp.annot[[i]]
   
   data.name <- "Calibrated Intensities"
   par(bty="n")
   if (verbose) {
      print(summary(snp.gt.p))
   }
   
   sigma.epsilon <- theta[["sigma.epsilon"]]
   sigma.eta <- theta[["sigma.eta"]]
   nu <- theta[["nu"]]
   het.cov <- theta[["het.cov"]]
   het.cov <- min(het.cov, (het.cov * e.mu[1] * e.mu[2]))
   mu0 <- 0
   sigma.h <- GetSigmaH(sigma.epsilon, sigma.eta)
   
   ct <- colSums(snp.gt.p)
   ct <- ct / sum(ct)
   ngt.t <- c((ct[1] + ct[3]),  (ct[2] + ct[4]))
   
   cex <- 1.2
   min <- -0.5
   max <- 5.25
   t.min <- -0.5
   t.max <- 2.75
   
   het.col <- "green"
   hom.col <- "black"
   
   df <- d[((d >= min) & (d < max))] 
   
   hist(df, breaks=100, freq=FALSE,
         main="", xlab="Allelic copy-ratio", xlim=c(min, max))
   
   xgl <- 1001
   xg <- array(NA, dim=c(4, xgl))
   x <- seq(min(df), max(df), length.out=xgl)
   
   ## transform grid/ params
   x.tx <- HTx(x, sigma.epsilon, sigma.eta) 
   e.mu <- HTx(e.mu, sigma.epsilon, sigma.eta)
   mu0.tx <- HTx(mu0, sigma.epsilon, sigma.eta)
   
   xg[1, ] <- ngt.t[1] * DFunc(x.tx, mu0.tx, sigma.h, nu) *
         HTxVar(x, sigma.epsilon, sigma.eta) 
   xg[2, ] <- ngt.t[2] * DFunc(x.tx, e.mu[1], sigma.h, nu) *
         HTxVar(x, sigma.epsilon, sigma.eta)
   xg[3, ] <- ngt.t[2] * DFunc(x.tx, e.mu[2], sigma.h, nu) *
         HTxVar(x, sigma.epsilon, sigma.eta)
   xg[4, ] <- ngt.t[1] * DFunc(x.tx, e.mu[3], sigma.h, nu) *
         HTxVar(x, sigma.epsilon, sigma.eta)
   
   csum <- apply(xg, 2, sum)
   
   crds <- par("usr")
   scale <- 1 / 2
   
   for (i in 1:4) {
      lines(x, scale * xg[i, ], lwd=2, col=hom.col)
   }
   lines(x, scale * csum, lwd=2, col="coral") 
   
   xg[1, ] <- ngt.t[1] * DFunc(x.tx, mu0.tx, sigma.h, nu)             
   xg[2, ] <- ngt.t[2] * DFunc(x.tx, e.mu[1], sigma.h, nu)
   xg[3, ] <- ngt.t[2] * DFunc(x.tx, e.mu[2], sigma.h, nu)
   xg[4, ] <- ngt.t[1] * DFunc(x.tx, e.mu[3], sigma.h, nu)
   csum <- apply(xg, 2, sum)
   
   d.tx <- HTx(d, sigma.epsilon, sigma.eta )
   d.tx.f <- d.tx[(d.tx > t.min) & (d.tx < t.max)] 
   hist(d.tx.f, breaks=100, freq=FALSE, main="", xlab="Allelic copy-ratio")
   mtext("Variance-stabilizing transformation", side=3, line=0, adj=0,
         cex=par("cex")*par("cex.axis"))
   
   for (i in 1:4) {
      lines(x.tx, scale * xg[i, ], lwd=2, col=hom.col)
   }
   lines(x.tx, scale * csum, lwd=2, col="coral") 
   
   ## 2d plots
   gl <- 250
   x2d <- seq(t.min, t.max, length.out=gl)
   x1 <- rep(x2d, gl)
   x2 <- as.vector(sapply(x2d, rep, gl))
   x <- cbind(x1, x2)
   hom.cov <- 0
   
   sigma <- matrix(c(sigma.h^2, het.cov, het.cov, sigma.h^2), nrow=2,
         ncol=2)
   invert.sigma <- solve(sigma)
   hom.sigma <- matrix(c(sigma.h^2, 0, 0, sigma.h^2), nrow=2, ncol=2)
   invert.hom.sigma <- solve(hom.sigma)
   
   het.1 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[1], e.mu[2]), sigma,
               invert.sigma, nu), dim=c(gl, gl))
   het.2 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[2], e.mu[1]), sigma,
               invert.sigma, nu), dim=c(gl, gl))
   
   hom.1 <- array(scale * ngt.t[1] * DmvFunc(x, c(mu0.tx, e.mu[3]), hom.sigma,
               invert.hom.sigma, nu), dim=c(gl, gl))
   hom.2 <- array(scale * ngt.t[1] * DmvFunc(x, c(e.mu[3], mu0.tx), hom.sigma,
               invert.hom.sigma, nu), dim=c(gl, gl))
   csum2d <- het.1 + het.2 + hom.1 + hom.2
   
   ## compute densities at modes
   het.md <- DmvFunc(c(e.mu[1], e.mu[2]), c(e.mu[1], e.mu[2]),
         sigma, invert.sigma, nu) * scale * ngt.t[2] 
   hom.md <- DmvFunc(c(mu0.tx, e.mu[3]), c(mu0.tx, e.mu[3]),
         hom.sigma, invert.hom.sigma, nu) * scale * ngt.t[1]
   
   level_scales <- c(0.99, 0.95, 0.8, 0.5, 0.2)
   
   out.pch <- rep(".", ncol(d))
   
   if (any(is.na(marker.gt.pal))) {
      pal <- GetHapsegMarkerPal()
   } else {
      pal <- marker.gt.pal
   }
   
   cols <- rev(pal(1000))
   state.mat <- matrix(c(0,1,2,3), ncol=4, nrow=nrow(snp.clust.p),
         byrow=TRUE)
   col.ix <- round(rowSums(state.mat * snp.clust.p[, c(3, 2, 4, 1)]) / 3 * 999, 0) + 1 
   pcols <- cols[col.ix]
   
   x <- d[1, ]
   y <- d[2, ]
   ix <- (x < max) & (x > min) & (y < max) & (y > min)
   if (sum(ix) > 5) {
      plot(0, type="n", main="", xlim=c(min, max), ylim=c(min, max),
            xlab="A allele copy-ratio", ylab="B allele copy-ratio", bty="n")
   } else {
      smoothScatter(x, y, main=data.name, xlab="A allele copy-ratio",
            ylab="B allele copy-ratio", bty="n")
   }
   
   points(d[1,], d[2,], col=pcols, pch=out.pch, cex=cex)
   cl.het = contourLines(x=x2d, y=x2d, z=het.1, levels=level_scales * het.md)
   cl.hom <- contourLines(x=x2d, y=x2d, z=hom.1, levels=level_scales * hom.md)
   
   for (i in seq_along(cl.het)) {
      txi.x <- HTxInv(cl.het[[i]][["x"]], sigma.epsilon, sigma.eta) 
      txi.y <- HTxInv(cl.het[[i]][["y"]], sigma.epsilon, sigma.eta)
      lines(txi.x, txi.y, col=het.col)
      
      txi.x <- HTxInv(cl.het[[i]][["y"]], sigma.epsilon, sigma.eta)
      txi.y <- HTxInv(cl.het[[i]][["x"]], sigma.epsilon, sigma.eta)
      lines(txi.x, txi.y, col=het.col)
   }
   
   for (i in seq_along(cl.hom)) {
      txi.x <- HTxInv(cl.hom[[i]][["x"]], sigma.epsilon, sigma.eta)
      txi.y <- HTxInv(cl.hom[[i]][["y"]], sigma.epsilon, sigma.eta)
      lines(txi.x, txi.y, col=hom.col)
      
      txi.x <- HTxInv(cl.hom[[i]][["y"]], sigma.epsilon, sigma.eta)
      txi.y <- HTxInv(cl.hom[[i]][["x"]], sigma.epsilon, sigma.eta)
      lines(txi.x, txi.y, col=hom.col)
   }
   
   abline(v=mu0, lty=2)
   abline(h=mu0, lty=2)
   
   ## tx data 2d fits
   d.tx <- HTx(d, sigma.epsilon, sigma.eta)
   plot(0, type="n", main="", xlim=c(t.min, t.max), ylim=c(t.min, t.max),
         xlab="A allele copy-ratio", ylab="B allele copy-ratio", bty="n")
   mtext("Variance-stabilizing transformation", side=3, line=0, adj=0,
         cex=par("cex")*par("cex.axis"))
   points(d.tx[1, ], d.tx[2, ], col=pcols, pch=out.pch, cex=cex)
   
   contour(x=x2d, y=x2d, z=het.1, add=TRUE, drawlabels=FALSE, levels=level_scales*het.md, col=het.col)
   contour(x=x2d, y=x2d, z=het.2, add=TRUE, drawlabels=FALSE, levels=level_scales*het.md, col=het.col)
   contour(x=x2d, y=x2d, z=hom.1, add=TRUE, drawlabels=FALSE, levels=level_scales*hom.md, col=hom.col)
   contour(x=x2d, y=x2d, z=hom.2, add=TRUE, drawlabels=FALSE, levels=level_scales*hom.md, col=hom.col)
   
   abline(v=0, lty=2)
   abline(h=0, lty=2)
   
   ## SNP signal by genome/GT/phase
   PlotSnpAscn(d, snp.clust.p, HTxInv(e.mu, sigma.epsilon, sigma.eta),
         "Allelic copy-ratio", xcrds=snp.annot[["pos"]] / 1e6,
         chr=snp.annot[["chr"]], seg.wd=1, loci.annot=merged.loci,
         merge.prob=merge.prob, seg.log.ev=seg.log.ev, min=min, max=max)
   
   PlotSnpAscn(HTx(d, sigma.epsilon, sigma.eta), snp.clust.p, e.mu,
         y.lab="Variance stabilized ", xcrds=snp.annot[["pos"]] / 1e6,
         chr=snp.annot[["chr"]], seg.wd=1, loci.annot=merged.loci,
         merge.prob=merge.prob, seg.log.ev=seg.log.ev, min=t.min,
         max=t.max)
}


## FIXME: needs heavy refactoring


GetHapsegMarkerPal <- function() {
   hom.color <- "darkgrey"
   het.1.color <- "red"
   het.2.color <- "blue"
   mid.color <- "darkviolet"
   
   pal <- colorRampPalette(c(hom.color, het.1.color, mid.color, het.2.color,
               hom.color))
   return(pal)
}
