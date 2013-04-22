## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

SummarizePlateStats <- function(plate.name, hapseg.files,
                                pdf.fn, seg.pdf.fn, save.fn,
                                calls.fn=NULL, verbose=FALSE) {
  if (verbose) {
    print("summarizing plate stats...")
  }

  if (is.null(calls.fn)) {
    call.rates <- NULL
  } else {
    call.rates <- read.delim(calls.fn, row.names=1, header=FALSE)
  }

  n.samples <- length(hapseg.files)
  
  seg.merge.prob <- vector("list", length=n.samples)

  sum.seg.log.ev <- rep(NA, n.samples)
  BG <- rep(NA, n.samples)
  sigma.eta <- rep(NA, n.samples)
  sigma.epsilon <- rep(NA, n.samples)
  sigma.h <- rep(NA, n.samples)
  het.cov <- rep(NA, n.samples)
  Nu <- rep(NA, n.samples)
  loglik <- rep(NA, n.samples)
  normal <- rep(NA, n.samples)
  array.name <- rep(NA, n.samples)
  AT <- rep(NA, n.samples)
  alpha <- rep(NA, n.samples)

  for (i in 1:n.samples) {
    cur.file <- hapseg.files[i]
    if (verbose) {
      print(cur.file)
    }
    load(cur.file)

    seg.merge.prob[[i]] <- seg.dat[["final.merge.prob"]]
    sigma.eta[i] <- seg.dat[["em.res"]][["theta"]][["sigma.eta"]]
    sigma.epsilon[i] <- seg.dat[["em.res"]][["theta"]][["sigma.epsilon"]]
    het.cov[i] <- seg.dat[["em.res"]][["theta"]][["het.cov"]]
    AT[i] <- seg.dat[["em.res"]][["theta"]][["at"]]
    sum.seg.log.ev[i] <- sum(seg.dat[["em.res"]][["seg.log.ev"]], na.rm=TRUE)
    BG[i] <- seg.dat[["em.res"]][["theta"]][["bg"]]
    alpha[i] <- seg.dat[["em.res"]][["theta"]][["bg"]]

    if (!is.null(seg.dat[["em.res"]][["theta"]][["nu"]])) {
      Nu[i] <- seg.dat[["em.res"]][["theta"]][["nu"]]
    }
    
    if (!is.null(seg.dat[["em.res"]][["loglik"]])) {
      loglik[i] <- seg.dat[["em.res"]][["loglik"]]
    }

    if (!is.null(seg.dat[["normal"]])) {
      normal[i] <- seg.dat[["normal"]] 
    } else{
      normal[i] <- FALSE
    }
    array.name[i] <- seg.dat[["array.name"]]

    if (verbose) {
      cat(".")
    }
  }
   
  ## fixes "TRUE" bug
  normal <- (normal == TRUE)

  if (!all(is.na(unlist(seg.merge.prob)))) {
    PlotSampleSegMergeProbs(seg.merge.prob, normal, array.name, seg.pdf.fn)
  }

  ## fit chi-sq dist
  if (sum(normal) > 5) {
    X <- sigma.epsilon[normal == TRUE]^2
    n.samples <- sum(normal == TRUE)
    res <- SInvChisqMl(X)

    pi.sigma.epsilon <- list(sigma=res[["sigma.hat"]], nu=res[["nu.hat"]],
                        sigma.n=n.samples, nu.n=n.samples)
    if (verbose) {
      print(pi.sigma.epsilon)
    }
  } else {
    pi.sigma.epsilon <- NA
  }

  ## fit het.cov vs. sigma.eta to a line 
  if (n.samples > 5 & sum(normal) > 5) {
    res <- lm(sqrt(het.cov[normal]) ~ sigma.eta[normal])
    coefs <- res[["coefficients"]]
  } else {
    coefs <- NA
  }

  stats <- cbind(sigma.epsilon, sigma.eta, sigma.h, het.cov, Nu, AT, alpha,
                 loglik, sum.seg.log.ev, BG)
  colnames(stats) <- c("sigma_nu", "sigma_eta", "sigma_h", "het_cov",
                       "Nu", "AT", "alpha", "loglik", "sum_seg_log_ev",
                       "BG")
  rownames(stats) <- array.name
  names(normal) <- array.name

  plate.stats <- list(stats=stats, normal=normal, pi.sigma.epsilon=pi.sigma.epsilon,
                      het.cov.coefs=coefs)

  save(plate.stats, file=save.fn)
      
   pdf(pdf.fn, 12, 12)
   par(mfrow=c(3,3))
   par(bty="n")
   
  plot(c(0,0), type='n', xlab="sigma.eta",
       ylab="Cumulative fraction of arrays", main="", xlim=range(sigma.eta),
       ylim=c(0,1))

  if (sum(normal==FALSE) > 0) {
    plot(ecdf(sigma.eta[normal==FALSE]), add=TRUE)
  }

  if (sum(normal==TRUE) > 0) {
    plot(ecdf(sigma.eta[normal==TRUE]), col.p=2, add=TRUE)
  }
   
  plot(c(0,0), type='n', xlab="sigma.epsilon",
       ylab="Cumulative fraction of arrays", main="",
       xlim=range(sigma.epsilon), ylim=c(0,1))

  if (sum(normal==FALSE) > 0) {
    plot(ecdf(sigma.epsilon[normal==FALSE]), add=TRUE)
  }

  ## chi-sq fit
  if (sum(normal) > 5) {
    dd <- RSInvChisq(1e3, pi.sigma.epsilon[["nu"]], pi.sigma.epsilon[["sigma"]])
    plot(ecdf(dd), col.p=3, add=TRUE)
    plot(ecdf(sigma.epsilon[normal==TRUE]), col.p=2, add=TRUE)
  }
  
  legend(0.2, 0.5, legend=c("Tumors", "Normals"), col=c(1,2), pch=1)
   
  use.col <- rep(1, nrow(stats))
  use.col[normal== TRUE] <- 2
  plot(sigma.eta, sigma.epsilon, xlab="sigma.eta", ylab="sigma.epsilon", col=use.col)
  plot(sigma.eta, sqrt(het.cov), xlab="sigma.eta", ylab="sqrt(het.cov)", col=use.col)
   
  ## fit het.cov vs. sigma.eta to a line 
  if (!is.na(coefs))  {
    abline(a=coefs[1], b=coefs[2], col=2, lty=2)
  }

  if (sum(normal) > 5) {
    if (!is.null(call.rates)) {
      plot(sigma.eta[normal==TRUE],
           call.rates[array.name[normal==TRUE], 1], xlab="sigma.eta",
           ylab="Call rate", col=2)
    }
  }

  plot(ecdf(Nu[normal==FALSE]), ylab="Cumulative fraction of arrays",
       xlab="Nu", main="", col.p=1, xlim=range(Nu))
  if (sum(normal == TRUE) > 0) {
    plot(ecdf(Nu[normal==TRUE]), add=TRUE, col.p=2)
  }

  plot(ecdf(loglik[normal==FALSE]),
       ylab="Cumulative fraction of arrays",
       xlab="Sample-fit loglik", main="", col.p=1, xlim=range(loglik))
  if (sum(normal == TRUE) > 0) {
    plot(ecdf(loglik[normal==TRUE]), add=TRUE, col.p=2)
  }
   
  plot(ecdf(AT[normal==FALSE]), ylab="Cumulative fraction of arrays",
       xlab="Sample-fit Attenuation", main="", col.p=1, xlim=range(AT))
  if (sum(normal == TRUE) > 0) {
    plot(ecdf(AT[normal==TRUE]), add=TRUE, col.p=2)
  }

  at.col <- rep(1, n.samples)
  at.col[normal==TRUE] <- 2
  plot(AT, sigma.eta, col=at.col, main="", xlab="AT", ylab="sigma.eta")
  plot(AT, sigma.epsilon, col=at.col, main="", xlab="AT", ylab="sigma.epsilon")
  
  plot(ecdf(BG[normal==FALSE]), ylab="Cumulative fraction of arrays",
       xlab="Sample-fit BackGround", main="", col.p=1, xlim=range(BG))
  if (sum(normal == TRUE) > 0) {
    plot(ecdf(BG[normal==TRUE]), add=TRUE, col.p=2)
  }
  
  dev.off()
}

PlotSampleSegMergeProbs <- function(seg.merge.prob, normal,
                                    array.names, pdf.fn) {
  pd <- 4
  pdf(pdf.fn, pd*3, pd*3)
  par(mfcol=c(pd, pd))

  for (i in seq_along(normal)) {
    if (normal[i]) {
      next
    }
    if (all(is.na(seg.merge.prob))) {
      next
    }

    mp <- seg.merge.prob[[i]][, 2]
    mp <- mp[is.finite(mp)]

    if (length(mp) == 0) {
      next
    }

    ## convert to log-log
    mp <- -log(-mp, 10)

    plot(ecdf(mp), col.p=1, xlab="-log10(-log10 merge prob)",
         ylab="Fraction of seg-junctions", bty="n", main="")
    mtext(text=array.names[i], cex=0.5, line=2, side=3, adj=0)
    abline(v=-1, lty=2)
  }
       
  dev.off()
}
