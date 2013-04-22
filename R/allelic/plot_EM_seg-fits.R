## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

PlotEmSegFits <- function(plate.name, hapseg.files, pdf.fn,
                          verbose=FALSE) {
  pdf(pdf.fn, 16, 12)
  par(mfrow=c(3, 2))
  if (verbose) {
    print(paste("PDF =", pdf.fn))
  }
  for (i in seq_along(hapseg.files)) {
    cur.result <- hapseg.files[i]
    if (verbose) {
      print(cur.result)
    }
    load(cur.result)
    if (seg.dat[["normal"]]) {
      next
    }
    array.name <- gsub("\\.RData", "", basename(cur.result))
    PlotSampleComparativeSegs(seg.dat, array.name)
    if (verbose) {
      cat(".")
    }
  }
  dev.off()
}

PlotSampleComparativeSegs <- function(seg.dat, array.name) {
  xmax <- 2.5
  as.max <- 2.5
  
  ## setup pallete for seg std.errors
  pal <- colorRampPalette( c("blue", "purple", "gray") )
  se.pal <- pal( 1000 )
  se.range <- c( 1e-3, 0.1 )  ## range for plotting
  por <- seg.dat[["seg.info"]][,"copy_num"]
  por.ix <-  por < xmax 

  ## Kinda-FIXME: Split this out into 3 lines, there
  ## was some weird bug that I couldn't spot that wasn't
  ## really catting them, yet c&p the code into the REPL
  ## worked
  as1 <- seg.dat[["allele.segs"]][, "A1.Seg.CN"]
  as2 <- seg.dat[["allele.segs"]][, "A2.Seg.CN"]
  as <- c(as1, as2)  
  
  as.ix <- (as < as.max) & (as >= 0)
  binW <- 0.025
  W <- seg.dat[["seg.info"]][, "length"]
  W <- W / sum(W)
  chr <- seg.dat[["seg.info"]][, "Chromosome"]

  PlotSeglenHist(por[por.ix], W[por.ix], chr[por.ix], order.by=chr[por.ix],
                 x.max=xmax, bin.w=binW,
                 use.pal=rainbow(length(unique(chr[as.ix]))),
                 xlab="Total copy-ratio")
  title("HAPSEG Pipeline result")
  
  ## plot with sample stats summary.
  plot( 0, type="n", bty="n", xlab="", ylab="", main=array.name, axes=FALSE )
  hscn.params <- seg.dat[["em.res"]][["theta"]]

  mtext(substitute(paste( sigma[eta] == x1, sep=""),
                   list(x1=round(hscn.params[["sigma.eta"]], 3))), 
        line=-1, adj=0, cex=1.25)
  
  mtext(substitute(paste(sigma[epsilon]==x2,  sep=""),
                   list(x2=round(hscn.params[["sigma.epsilon"]], 3))),
        line=-2.5, adj=0, cex=1.25)
  mtext(substitute(paste( nu==x4,  sep=""), list(x4=round(hscn.params[["nu"]], 3))),
        line=-5.5, adj=0, cex=1.25)
  
  mtext(paste("het_cov = ", round(hscn.params[["het.cov"]], 4), sep=""),
        line=-7, adj=0, cex=1.25)
  mtext(paste("AT = ", round(hscn.params[["at"]], 3), sep="" ),
        line=-8.5, adj=0, cex=1.25)
  mtext(paste("BG = ", round(hscn.params[["bg"]], 3), sep="" ),
        line=-10, adj=0, cex=1.25)
  msg <- paste("sum_seg_log_ev = ",
               round(sum(seg.dat[["em.res"]][["seg.log.ev"]], na.rm=TRUE),
                     3), sep="" )
  mtext(msg, line=-11.5, adj=0, cex=1.25)
  W <- rep(seg.dat[["allele.segs"]][, 5], 2)
  W <- W / sum(W)
  chr <- rep(seg.dat[["allele.segs"]][, 1], 2)

  PlotSeglenHist(as[as.ix], W[as.ix], color.by=chr[as.ix],
                 order.by=chr[as.ix], x.max=as.max, bin.w=binW,
                 use.pal=rainbow(length(unique(chr[as.ix]))))
  title("HSCR estimates")
  
  AS.seg.sd <- rep(seg.dat[["allele.segs"]][, "AS.Seg.sd"], 2)
  AS.seg.sd[AS.seg.sd < se.range[1]] <- se.range[1]
  AS.seg.sd[AS.seg.sd > se.range[2]] <- se.range[2]
  PlotSeglenHist(as[as.ix], W[as.ix], color.by=log(AS.seg.sd[as.ix]),
                 color.range=log(se.range), use.pal=se.pal,
                 order.by=-W[as.ix], x.max=as.max, bin.w=binW)
  title("HSCR estimates")
  
  CN <- seg.dat[["em.res"]][["e.mu"]][, 3] / 2
  cn.ix <- CN < xmax & CN >= 0
  W <- seg.dat[["allele.segs"]][, 5]
  W <- W / sum(W)
  chr <- seg.dat[["allele.segs"]][, 1]
  PlotSeglenHist(CN[cn.ix], W[cn.ix], color.by=chr[cn.ix],
                 order.by=chr[cn.ix], x.max=xmax, bin.w=binW,
                 use.pal=rainbow(length(unique(chr[cn.ix]))),
                 xlab="Total copy-ratio")
  title("Total copy-ratio estimates") 
  
  tcn.seg.sd <- seg.dat[["allele.segs"]][, "tCN.Seg.sd"]
  tcn.seg.sd[tcn.seg.sd < se.range[1]] <- se.range[1]
  tcn.seg.sd[tcn.seg.sd > se.range[2]] <- se.range[2]

  PlotSeglenHist(CN[cn.ix], W[cn.ix], color.by=log(tcn.seg.sd[cn.ix]),
                 color.range=log(se.range), use.pal=se.pal,
                 order.by=-W[cn.ix], x.max=xmax, bin.w=binW,
                 xlab="Total copy-ratio")
  return(TRUE)
}




