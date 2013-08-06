## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

FastDmvt <- function(x, delta, sigma, invert.sigma=NULL, df=1, log=TRUE) {
#  if (is.null(invert.sigma)) {
#    invert.sigma <- solve(sigma)
#  }
  
  if (df == 0) {
    return(DmvNorm(x, delta, sigma, log=log))
  }

  sigma.ncol <- NCOL(sigma)

 # dist.val <- mahalanobis(x, delta, invert.sigma, inverted=TRUE)
  dist.val <- mahalanobis(x, delta, sigma)
  log.det <- sum(log(eigen(sigma, TRUE, only.values=TRUE)[["values"]]))
  log.retval <- lgamma((sigma.ncol + df) / 2) - (lgamma(df / 2) + 0.5 *
                                                 (log.det + sigma.ncol * logb(pi * df))) -
                                                   0.5 * (df + sigma.ncol) *
                                                     logb(1 + dist.val / df)
  if (!log) {
    log.retval <- exp(log.retval)
  }
  
  return(log.retval)
}


LogAdd <- function(X) {
     ##  Calculates log( sum(exp(x)) )  without 'leaving' log space
    if (is.vector(X)) {
        mix <- which.max(X)
        max <- X[mix]
        res <- max + log(sum(exp(X - max)))
    }
    
    if (is.matrix(X)) {
        mv <- apply(X, 1, max)
        res <- mv + log(rowSums(exp(X - mv)))
    }
    
    return(res)
}

dDirichlet <- function(X, alpha, log.p = FALSE) {
    if (isTRUE(all.equal(sum(X), 1))) {
        ix <- is.finite(log(X))
        X[!ix] <- .Machine[["double.xmin"]]
        
        log.lik <- sum(((alpha - 1) * log(X)))
        log.z <- sum(lgamma(alpha)) - lgamma(sum(alpha))
        log.dens <- sum(log.lik) - log.z
    } else {
        log.dens <- -Inf
    }
    
    if (log.p) {
        return(log.dens)
    } else {
        return(exp(log.dens))
    }
}



dPolya = function( N, alpha, log.p=FALSE )
{
  ## aka. the Dirichlet-multinomial distribution

  ## supports vectorization if N is a matrix 
   if( is.matrix(N) )
   {
      alpha_sum = sum(alpha)
      alpha = matrix( alpha, nrow=nrow(N), ncol=K, byrow=TRUE )
   
      logZ = lgamma(alpha_sum) - lgamma( rowSums(N + alpha) )
      LL = rowSums( lgamma(N+alpha) - lgamma(alpha) )
   }
   else
   {
      logZ = lgamma(sum(alpha)) - lgamma( sum(N + alpha) )
      LL = sum( lgamma(N+alpha) - lgamma(alpha) )
   }

   if(log.p) { return(LL-logZ) }
   else
   {
      return( exp(LL-logZ) )
   }
}


weighted_sample_sd = function( c_data, c_mean, c_var )
{
   c_SD = sqrt( sum( ((c_data-c_mean)^2 + c_var) / c_var) / sum( 1/c_var ) )
   return(c_SD)
}



d_scaled_t =  function( x, mu, sigma, nu, log=FALSE )
{
   log_norm = log( gamma((nu+1)/2) / ( gamma(nu/2) * sqrt(nu*pi) * sigma ) )

   log_d = log(1+ ((x-mu)/sigma)^2 * 1/nu ) * -(nu+1)/2 + log_norm

   if( log==FALSE )
   {
      res = exp(log_d)
   }
   else
   {
      res = log_d 
   }

   return(res) 
}


d_beta_binom = function( k, a, b, n, log=FALSE )
{
  ## k must be an integer for this to be correct!!
  ## a and b are prior sample sizes (the prior beta params)
   log_d = (lbeta( k + a, n - k + b) + lchoose(n,k)) - lbeta(a,b)
   # if (is.nan(log_d)) log_d <- -Inf
   if( log==TRUE ) { return( log_d ) }
   else { return( exp(log_d) ) }
}



p_beta_binom = function( k, a, b, n, lower.tail=TRUE )
{
  p = rep(NA, length(k) )
 
  for( i in 1:length(k) )
  {
     p[i] = sum( d_beta_binom( c(0:k[i]), a, b, n[i] ) )
  }

  if( !lower.tail ) { p = 1-p }
  p[p<0]=0

  return(p)
}



d_igamma = function( x, alpha, beta )
{
   d = (beta^alpha)/gamma(alpha) *  x^(-alpha-1) * exp(-beta/x)

   return(d)
}

d_S_inv_chisq = function( x, nu, sigma )
{
   alpha = nu/2
   beta = alpha * sigma^2

   d = d_igamma( x, alpha, beta )

   return(d)
}



S_inv_chisq_joint_MAP = function( X, pi_nu, pi_nu_N, pi_sigma, pi_sigma_N )
{
   pi_alpha = pi_nu / 2                 ## init alpha (gamma shape)
   pi_beta = ( pi_alpha * pi_sigma^2 )^-1    ## init beta (inverse-scale of gamma dist)

   pi_alpha_N = pi_nu_N
   pi_beta_N = pi_sigma_N

   res = gamma_joint_MAP( 1/X, pi_alpha, pi_alpha_N, pi_beta, pi_beta_N )

   alpha = res$alpha
   scale = res$scale
   beta = scale^-1

   post_nu = alpha*2
   post_sigma = sqrt(beta / alpha)

   return( list("sigma"=post_sigma, "nu"=post_nu ))
}



softmax = function( a, log=FALSE )
{
   if( log == FALSE )
   {
#      val = exp(a)/ sum(exp(a))
      val = exp( a - LogAdd(a) )
   }
   else
   {
      val = a - LogAdd(a)
   }

   return( val )
}


inv_softmax = function( p, PC=0 )
{
   p = p+PC  ## pseudo-count
   p = p / sum(p)
  
   lv = log(p) 
   lv[ !is.finite(lv) ] = log(.Machine$double.xmin )
   mv = mean(lv)

   a = lv - mv

   return( a )
}


Laplace_approx = function( Hess_mat, mode_value )
{
   d = nrow(Hess_mat)
   inv_neg_H = solve(Hess_mat, diag(d) )
   logI = (d/2)*log(2*pi) + 0.5*determinant(inv_neg_H, log=TRUE)$modulus[1] + mode_value

   if ((!is.finite(logI)) ) 
   {
     es = eigen(inv_neg_H, only.values=TRUE, symmetric=TRUE )$values
     d_approx = product(es[es>0])
     logI = (d/2)*log(2*pi) + 0.5 * log(d_approx) + mode_value

     print( paste("WARNING: NON-FINITE Laplace approx in SCNA_model - approx: ", round(logI,4), sep="") )
   }
   
   return(logI)


# SLC orig:   
#   curvature <- abs(det(hess.mat / (2 * pi)))
#   mode.curv <- log((curvature)^(-1/2) ) - log(2)

# Kemp 2008:
#% Laplace approximation to p(D|S)
#% We approximate the integral over branch lengths
#% minus sign because datal computes -ll
#H=-hessiangrad(datal, X, 1e-5);
#
#includeind = find(X<upper_bound - 5);
#H = H(includeind, includeind);
#d = length(includeind);
#if includeind(1) ~= 1
#  disp('WARNING: sigma blows up');
#end
#
#logI = (d/2)*log(2*pi)+0.5*mylogdet(inv(-H))+ll;
#if ~isreal(logI)
#  disp('WARNING: laplacian approx gone awry');
#  es = eig(inv(-H));
#  dapprox = real(prod(es(es>0)));
#  logI = (d/2)*log(2*pi)+0.5*log(dapprox)+ll;
#end
}

