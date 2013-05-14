## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

SegInitTheta <- function() {
  list(eps=1e-3, p.snp.cond.out=1, nu=0.5, out.p=0.01)
}

SegThetaOpt <- function(h.d, h.snp.gt.p, h.e.mu, out.p, theta) {
  theta[["nu"]] <- HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta,
                             "nu", limits=list("lower"= 0.1, "upper"=1.0),
                             symbol="%", opttol=1e-3)
  theta[["out.p"]] <- HThetaOpt(h.d, h.snp.gt.p, h.e.mu, out.p, theta,
                                "out.p",
                                limits=list("lower"= 0.001, "upper"=0.5),
                                symbol="#", opttol=1e-3)

  print(paste("theta nu =", round(theta[["nu"]], 3)))
  print(paste("theta out.p =", round(theta[["out.p"]], 3)))

  theta
}

SegAtten <- function(r, at) {
  return(r)
}

SegTauPlatformSpecificInitialization <- function(h.d, seg.info) {
  ## platform specific initialization for optim.
  ## This function returns a three part data frame, with the
  ## following values for tau:
  ## ( 0.25 point, 0.75 point, 1.0 point)
  # print(summary(seg.info))

  twoTimesSegSix <- 2 * seg.info[, 6]

  cbind(twoTimesSegSix * 0.25,
        twoTimesSegSix * 0.75,
        twoTimesSegSix)
}

SegPlatformSpecificOptimization <- function(cur.par, d, out.p, snp.gt.p, theta) {
  ## platform specific likelihood function for the WES dataset
  PrivateLL <- function(f, d, out.p, snp.gt.p, theta, tau) {
    ## FIXME: This looks like it should be ||
    if ((f < 0) | (f > 1)) {
      return (Inf)
    }
    e.mu <- c(f * tau, (1 - f) * tau, tau)
    SegGetLLInternal(e.mu, d, out.p, snp.gt.p, theta)
  }

  res <- optimize(interval=c(0, 0.5), maximum=TRUE, f=PrivateLL,
                  d=d, out.p=out.p, snp.gt.p=snp.gt.p, theta=theta,
                  tau=cur.par[3])
  GetMeans(f=res[["maximum"]], cur.par[3])
}

SegGetMeans <- function(f, t) {
  c(f * t, (1 - f) * t, t)
}

SegGetLL <- function(e.mu, d, out.p, snp.gt.p, theta, tau) {
  ## Platform specific likelihood function for the WES dataset
  if ((e.mu < 0) || (e.mu > 1)) {
    Inf
  } else {
    SegGetLLInternal(e.mu, d, out.p, snp.gt.p, theta)
  }
}

SegGetLLInternal <- function(e.mu, d, out.p, snp.gt.p, theta) {
  ## using the SNP log likelihood function, calculate the likelihood of
  ## the data given the e.mu
  ll <- sum(CalcSnpLogLik(d, e.mu, out.p, snp.gt.p, theta))
  ## if we got a non-numeric likelihood, return negative infinity
  ifelse(is.nan(ll), -Inf, ll)
}

SegCalcSnpLogLik <- function(d, e.mu, dummy, snp.gt.p, theta) {
  ## This function calculates the SNP log likelihood, given
  ## the following parameters:
  ## d - the input data, containing a row for each allele (A & B)
  ## e.mu - the error model
  DoCalcSnpLogLik(d, e.mu, theta[["out.p"]], snp.gt.p, theta)
}

SegGetSnpClustLik <- function(d, e.mu, theta) {
   n.col <- ncol(d)
   clust.lik <- array(0, dim=c(n.col, 5))
   d <- t(d)

   if (n.col == 1) {
     ## for the one case, we can't calculate the beta density
     return(NA)
   }

   eps <- theta[["eps"]]
   f <- e.mu[1] / e.mu[3]
   clust.lik[, 1] <- DmvFunc(d, c(eps, 1 - eps), theta)
   clust.lik[, 2] <- DmvFunc(d, c(1.0 - f, f), theta)
   clust.lik[, 3] <- DmvFunc(d, c(1 - eps, eps), theta)
   clust.lik[, 4] <- DmvFunc(d, c(f, 1.0 - f), theta)
   clust.lik[, 5] <- theta[["p.snp.cond.out"]]
   clust.lik
}

SegDmvFunc <- function(x, f, theta, log=FALSE) {
  ## calculates the beta density based likelihood
  cov <- rowSums(x)
  allele.fract <- x[, 1] / cov

  ##  beta density is not defined at 0 or 1.
  allele.fract[allele.fract == 0] <- 1e-3
  allele.fract[allele.fract == 1] <- 1-1e-3

  dbeta(allele.fract, f * (cov * theta[["nu"]]), (1.0 - f) * (cov * theta[["nu"]]))
}

PlotWesSegFit <- function(h.d, h.e.mu, theta, segment.info,
                          signal, seg.number) {
  print(paste("Processing Segment:", seg.number, "SNP =",
              length(h.d[1, ]), "start = ",
              segment.info[["start.bp"]], "end = ",
              segment.info[["end.bp"]], "Chromosome =",
              segment.info[["hromosome"]]))

  data <- h.d[1, ] / (h.d[2, ] + h.d[1, ])
  data.total <- h.d[2, ] + h.d[1, ]
  f.value <- h.e.mu[1] / h.e.mu[3]
  n <- median(data.total) * theta[["nu"]]
  if (length(data) > 1) {
    sq <- seq(0, 1, by=0.001)
    ## get the delta
    my.plot2 <- hist(data, breaks=15, ylab="allele fraction",
                     col=rgb(.8, .8, .8, .3),
                     main=paste("all data - segment", h.e.mu),
                     xlim=c(0, 1))
    one.beta <- dbeta(sq, (1 - f.value) * n, (f.value) * n)
    two.beta <- dbeta(sq, (f.value) * n, (1 - f.value) * n)

    one.max <- max(one.beta[is.finite(one.beta)])
    two.max <- max(two.beta[is.finite(two.beta)])
    lines(x=sq, y=(one.beta / one.max * (max(my.plot2[["counts"]]))),
          col="red", lwd=4)
    lines(x=sq, y=(two.beta / two.max * (max(my.plot2[["counts"]]))),
          col="blue", lwd=4)
  }
  dt <- t(h.d)
  vals <- (dt[, 1] / (dt[, 1] + dt[, 2]))

  plot(vals, xlab="Position", ylab="Allelic Fraction",
       col=rgb(0, 0, 0, 0.5), ylim=c(0, 1), pch=18,
       main=paste(segment.info[["Chromosome"]],
         ":", segment.info[["start.bp"]], "-", segment.info[["end.bp"]]))
  abline(h=(1 - f.value), col="red")
  abline(h=(f.value), col="blue")

  pl.data <- exp(signal[(signal[["contig"]] == segment.info[["chromosome"]]) &
                        (signal[["start"]] >= segment.info[["start.bp"]]) &
                        (signal[["stop"]] <= segment.info[["end.bp"]]), 5])
  color <- max(1 / length(pl.data), 0.25)
  if (length(pl.data) > 0) {
    plot(pl.data,pch=18, col=rgb(1, 0, 0, color), ylim=c(0, 2),
         main=paste("region", segment.info[["chromosome"]],
           ":", segment.info[["start.bp"]], "-", segment.info[["end.bp"]]))
  } else {
    plot(1, main="No points")
  }
}
