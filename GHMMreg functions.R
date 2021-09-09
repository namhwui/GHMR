# wrapper function for ipar
# inserts equal weights if weight is not given
rgpar <- function(y=NULL, p=NULL, G=NULL, w=NULL) {
  if (is.null(w)) {
    w <- matrix(1/G, nrow=length(y), ncol=G)
  }
  
  val    <- list()
  for (g in 1:G) {
    val[[g]] <- ipar(y=y, p=p, wt=w[,g])
  }
  val$pi <- rep(1/G, G)
  
  return(val)
}	

# larger wrapper function for initialisation
# based on kmeans clustering
kmeansgpar <- function(y=NULL, x=NULL, G=NULL, n=10, label=NULL, mtol=NULL, mmax=NULL, skewness=FALSE) {
  
  p    <- ncol(x)
  lw   <- kmeans(x=y, centers=G, iter.max=100, nstart=10)$cluster
  w    <- combinewk(weights=matrix(0, nrow=length(y), ncol=G), label=lw)
  gpar <- rgpar(y=y, p=p, G=G, w=w)
  
  #gpar$gamma <- as.matrix(rep(0.1, ncol(x)), ncol=1)
  
  for (j in 1:n) try({gpar <- EMgrstep(y=y, x=x, gpar=gpar, label=label, w=w, mtol=mtol, mmax=mmax, skewness=skewness)}, TRUE)
  try({gpar <- EMn(y=y, x=x, gpar=gpar, label=label, mtol=mtol, mmax=mmax, skewness=skewness )$gpar}, TRUE)
  
  return(gpar)
}

# uniform initialisation
# because regmixEM is a dumbass and doesn't use kmeans
unifgpar <- function(y=NULL, x=NULL, G=NULL, n=10, label=NULL, mtol=NULL, mmax=NULL, skewness=FALSE) {
  
  p    <- ncol(x)
  #lw   <- kmeans(x=y, centers=G, iter.max=100, nstart=10)$cluster
  lw <- sample(1:G, length(y), replace=TRUE)
  w    <- combinewk(weights=matrix(0, nrow=length(y), ncol=G), label=lw)
  gpar <- rgpar(y=y, p=p, G=G, w=w)
  
  #gpar$gamma <- as.matrix(rep(0.1, ncol(x)), ncol=1)
  
  for (j in 1:n) try({gpar <- EMgrstep(y=y, x=x, gpar=gpar, label=label, w=w, mtol=mtol, mmax=mmax, skewness=skewness)}, TRUE)
  try({gpar <- EMn(y=y, x=x, gpar=gpar, label=label, mtol=mtol, mmax=mmax, skewness=skewness )$gpar}, TRUE)
  
  return(gpar)
}


logbesselKv <- function(x, y) {
  log(besselK(x=y, nu=x, expon.scaled=TRUE)) - log(y)
}


besselKv <- function(x, y) {
  besselK(x=y, nu=x)
}


gig <- function(y=NULL, x=NULL, par=NULL) {
  # returns a matrix with dim length(a) x 3
  alpha  <- par$alpha
  sigma  <- par$sigma
  omega  <- par$ol[1]
  lambda <- par$ol[2]
  mu     <- x %*% par$gamma
  #print('here 2')
  a <- omega + alpha^2/sigma
  b <- omega + (y - mu)^2/sigma
  v <- lambda - 1/2
  
  sab  <-  sqrt(a * b) 
  kv1  <- besselK(sab, nu=v+1, expon.scaled=T)
  kv   <- besselK(sab, nu=v, expon.scaled=T)
  kv12 <- kv1/kv
  
  sb.a <- sqrt(b/a)
  w    <- kv12 * sb.a
  invw <- kv12/sb.a - 2 * v/b
  logw <- log(sb.a) + grad(logbesselKv, x=rep(v,length(sab)), y=sab, 
                           method="Richardson",  
                           method.args=list(eps=1e-8, d=0.0001, 
                                            zero.tol=sqrt(.Machine$double.eps/7e-7), 
                                            r=6, v=2, show.details=FALSE))
  #print(kv12)
  #print(cbind(kv1,kv))
  #print(kv1[is.nan(kv12), ])
  val <- cbind(w, invw, logw)	
  return(val)
}


weighted.sum <- function(z=NULL, wt=NULL, ...) {
  return(sum(z*wt, ...))
}

update.ggam <- function(y=NULL, x=NULL, abc=NULL, w=NULL, par=NULL) {
  
  M <- diag(w * abc[, 2])
  
  A <- solve(t(x) %*% M %*% x)
  B <- t(x) %*% as.matrix(M%*%y + w*par$alpha, ncol=1)
  
  
  val <- A %*% B
  return(val)
}

ggam.unbiased <- function(gamma=NULL, x=NULL, par=NULL, w=NULL) {
  
  M <- diag(w * par$abc[, 2])
  A <- solve(t(x) %*% M %*% x)
  B <- t(x) %*% as.matrix(par$alpha * (M%*%par$abc[, 1] - w), ncol=1)
  
  val <- gamma - A %*% B
  return(val)
}


update.maRol <- function(y=NULL, x=NULL, par=NULL, weights=NULL, alpha.known=NULL) {
  
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  }
  
  # expectations of w, 1/w, log(w) given x	
  abc <- gig(y=y, x=x, par=par)
  #print(abc)
  
  sumw <- sum(weights)
  ABC  <- apply(abc, 2, weighted.sum, wt=weights)/sumw
  
  gam.new <- update.ggam(y=y, x=x, abc=abc, w=weights, par=par)
  
  if (is.null(alpha.known)) {
    A <- ABC[1]
    B <- ABC[2]
    u <- (B - abc[, 2]) * weights
    t <- (A*abc[, 2] - 1) * weights
    t2 <- sum(t)
    #t2 <- sum(weights * abc[, 1])
    #mu.new    = apply(x, 2, weighted.sum, wt=t)/t2
    #alpha.new <- apply(x, 2, weighted.sum, wt=u)/t2
    alpha.new <- sum((y - x %*% gam.new) * weights)/sum(weights * abc[, 1])
  } 
  
  else { 
    alpha.new <- alpha.known
    #mu.new   = apply(x, 2, weighted.mean, w=abc[,2]*weights) - alpha.new/ABC[2]
  }	
  
  #alpha.new <- alpha.new * v
  ol.new <- update.ol(ol=par$ol, ABC=ABC, n=2)
  
  A <- sum((abc[, 2]*weights) * (y - x %*% gam.new)^2)/sumw
  r <- sum((y - x %*% gam.new) * weights)/sumw
  R <- A - 2*r*alpha.new + ABC[1]*alpha.new^2
  
  new.par <- list(alpha=alpha.new, sigma=R, ol=ol.new, gamma=gam.new, wg=abc[, 1], abc=abc)
  return(new.par)
}


RRlamz <- function(x,lam=NULL,z=1) {
  val =Rlam(x, lam=-lam)+ Rlam(x, lam= lam)
  zval = val - z
  return(zval)
}

Rlam <- function(x, lam=NULL) {
  v1 = besselK(x, nu= lam+1, expon.scaled=FALSE) 
  v0 = besselK(x, nu= lam, expon.scaled=FALSE) 
  val = v1/v0
  return(val)
}

update.ol <- function(ol=NULL, ABC=NULL, n=1) {
  
  for (i in 1:n) {
    # lambda 
    #print(ol)
    if (ABC[3] == 0) {
      ol[2] = 0
    } else {
      bv = grad( logbesselKv, x=ol[2], y=ol[1], method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      ol[2] = ABC[3]*(ol[2]/bv)
    }
    
    lam0 = ol[2]
    omg0 = ol[1]
    
    Rp = Rlam(omg0,lam=+lam0)
    Rn = Rlam(omg0,lam=-lam0)
    
    f1 = Rp + Rn - (ABC[1]+ABC[2])
    f2 = ( Rp^2 - (2*lam0+1)/omg0 *Rp -1 ) + ( Rn^2 - (2*(-1*lam0)+1)/omg0 *Rn -1 ) 
    #print(c(Rp,Rn,ABC[1], ABC[2]))
    if ( ol[1] > f1/f2 ) ol[1] = ol[1] - f1/f2
  }
  return(ol)
}

EMgrstep <- function(y=NULL, x=NULL, gpar=NULL, label=NULL, w=NULL, mtol=NULL, mmax=NULL, skewness=TRUE) {
  # data = design matrix (centralised, so no columns of 1)
  # x = residuals
  
  if (is.null(w)) {
    w <- weights(y=y, x=x, gpar=gpar)
  }
  if (!is.null(label)) {
    w <- combinewk(weights=w, label=label)
  }
  #print('weight')
  
  G  <- length(gpar$pi);
  Sk <- numeric(G)
  
  for (g in 1:G) {
    if (skewness) {
      gpar[[g]] <- update.maRol(y=y, x=x, par=gpar[[g]], weights=w[, g], alpha.known=NULL)
    } else {
      gpar[[g]] <- update.maRol(y=y, x=x, par=gpar[[g]], weights=w[, g], alpha.known=0)
    }
  }
  
  gpar$pi <- apply(w,2,mean)
  
  
  return(gpar)
}


getall <- function(loglik=NULL, eps=NULL) {
  
  if (length(loglik) < 3) {
    stop("must have at least 3 likelihood values")
  }
  
  n       <- length(loglik)
  lm1     <- loglik[n]
  lm      <- loglik[(n - 1)]
  lm_1    <- loglik[(n - 2)]
  am      <- (lm1 - lm)/(lm - lm_1)
  lm1.Inf <- lm + (lm1 - lm)/(1 - am)
  val     <- lm1.Inf - lm	
  
  continue <- TRUE
  if (is.nan(val)) {
    val <- 0
  }
  if (val < eps & val >= 0) {
    continue <- FALSE
  }
  
  return(continue)
}


EM <- function(y=NULL, x=NULL, gpar0=NULL, G=2, max.iter=100, 
               epsilon=1e-2, print.loglik=FALSE, label=NULL, 
               nstart=0, mtol=1e-8, mmax=10, skewness=TRUE,
               gam.unbiased=FALSE) {
  
  
  gpar   <- gpar0
  loglik <- numeric(max.iter)
  
  for (i in 1:3) {
    #print(i)
    loglik[i] <- llik(y=y, x=x, gpar=gpar)
    #print(dim(x))
    gpar      <- EMgrstep(y=y, x=x, 
                          gpar=gpar, label=label, 
                          mtol=mtol, mmax=mmax, skewness=skewness)
  }
  
  
  while (getall(loglik=loglik[1:i], eps=epsilon) & (i < max.iter)) {
    i         <- i + 1
    gpar      <- EMgrstep(y=y, x=x,
                          gpar=gpar, label=label, 
                          mtol=mtol, mmax=mmax, skewness=skewness)
    loglik[i] <- llik(y=y, x=x, gpar=gpar)
    
  }
  
  val <- list(loglik=loglik[1:i], 
              gpar=gpar, 
              z=weights(y=y, x=x, gpar=gpar), 
              map=MAP(y=y, x=x, gpar=gpar, label=label), 
              skewness=skewness)
  
  if (gam.unbiased) {
    w <- val$z
    for (g in 1:G) {
      val$gpar[[g]]$gamma.unb <- ggam.unbiased(gamma=val$gpar[[g]]$gamma, x=x, par=val$gpar[[g]], w=w[, g])
    }
  }
  
  return(val)
}

llik <- function(y=NULL, x=NULL, gpar=NULL) {
  logz <- matrix(0, nrow=length(y), ncol=length(gpar$pi))
  
  for (g in 1:length(gpar$pi)) {
    logz[, g] <- dghyp(y=y, x=x, par=gpar[[g]], logd=TRUE)
  }
  
  val <- sum(log(apply(logz, 1, function(z, wt=NULL) {
    return(sum(exp(z) * wt))
  },wt=gpar$pi)))
  
  return(val)
}


# Computes density for univariate GHD
dghyp <- function(y=NULL, x=NULL, par=NULL, logd=FALSE) {
  # univariate GHD density
  sigma  <- par$sigma
  alpha  <- par$alpha
  omega  <- par$ol[1]
  lambda <- par$ol[2]
  #print(x)
  mu     <- x %*% par$gamma
  #print('here')
  
  pa  <- omega + alpha^2/sigma
  mx  <- omega + (y - mu)^2/sigma
  kx  <- sqrt(mx * pa)
  
  lvx     <- matrix(0, nrow=length(y), 3)
  lvx[,1] <- (lambda - 1/2) * log(kx)
  lvx[,2] <- log(besselK(kx, nu=lambda-1/2, expon.scaled=T)) - kx
  lvx[,3] <- (y - mu) * alpha/sigma 
  
  lv    <- numeric(3)
  lv[1] <- -1/2*(log(sigma) + log(2*pi))
  lv[2] <-  -log(besselK(omega, nu=lambda, expon.scaled=FALSE)) 
  lv[3] <- (1/2 - lambda) * log(pa)
  
  val <- apply(lvx, 1, sum) + sum(lv)
  
  if (!logd) {
    val <- exp(val)
  }
  #print(length(val))
  return(val)
}


kmeansgpar <- function(y=NULL, x=NULL, G=NULL, n=10, label=NULL, mtol=NULL, mmax=NULL, skewness=FALSE) {
  
  p    <- ncol(x)
  lw   <- kmeans(x=y, centers=G, iter.max=100, nstart=5)$cluster
  w    <- combinewk(weights=matrix(0, nrow=length(y), ncol=G), label=lw)
  gpar <- rgpar(y=y, p=p, G=G, w=w)
  
  #gpar$gamma <- as.matrix(rep(0.1, ncol(x)), ncol=1)
  
  for (j in 1:n) try({gpar <- EMgrstep(y=y, x=x, gpar=gpar, label=label, w=w, mtol=mtol, mmax=mmax, skewness=skewness)}, TRUE)
  try({gpar <- EMn(y=y, x=x, gpar=gpar, label=label, mtol=mtol, mmax=mmax, skewness=skewness )$gpar}, TRUE)
  
  return(gpar)
}


EMn <- function(y=NULL, x=NULL, gpar0=NULL, G=2, n=10, label=NULL, nstart=0, mtol=1e-8, mmax=10, skewness= TRUE) {
  #if (is.null(gpar0)) gpar = igpar(data=data, g=G, nstart=nstart, mtol= mtol, mmax= mmax, covtype= covtype, skewness= skewness )
  #else gpar  = gpar0
  gpar   <- gpar0
  loglik <- numeric(n)
  
  for (i in 1:n) {
    #print(i)
    gpar      <- EMgrstep(y=y, x=x,
                          gpar=gpar, label=label, 
                          mtol=mtol, mmax=mmax, skewness=skewness)		
    #print('here 1')
    loglik[i] <- llik(y=y, x=x, gpar=gpar)
    #print('here 2')
  }
  #print(gpar[[1]]$gamma)
  #print(gpar[[2]]$gamma)
  #print(gpar[[3]]$gamma)
  #print(MAP(y=y, x=x, gpar=gpar, label=label))
  #print('here 3')
  val <- list(loglik=loglik, 
              gpar=gpar, 
              z=weights(y=y, x=x, gpar=gpar), 
              map=MAP(y=y, x=x, gpar=gpar, label=label), 
              skewness=skewness)
  
  return(val)
}


MAP <- function(y=NULL, x=NULL, gpar=NULL, label=NULL) {
  w <- weights(y=y, x=x, gpar=gpar)
  if (!is.null(label)) {
    w <- combinewk(weights=w, label=label)
  }
  
  z <- apply(w, 1, function(z) {
    z <- (1:length(z))[z == max(z)]
    return(z[1]) 
  })
  z <- as.numeric(z)
  
  return(z)	
}

# combines weight with label info
combinewk <- function(weights=NULL, label=NULL)	{
  # known is a numeric with 
  # 0 if unknown group membership 
  # 1,2,3,.. for label of known group
  if (is.null(label)) {
    stop('label is null')
  }
  
  kw <- label != 0
  for (j in 1:ncol(weights)) {
    weights[kw, j] <- (label == j)[kw]
  }
  
  return(weights)	
}

ipar <- function(y=NULL, p=NULL, wt=NULL) {
  
  n <- length(y)
  if (is.null(wt)) {
    wt <- rep(1, n)
  }
  
  val       <- list()
  val$alpha <- rnorm(1, 0, sqrt(1/n))
  val$sigma <- (n-1)/n * sd(sqrt(wt) * y)^2 + sd(y)^2 + val$alpha^2
  val$ol    <- c(1, -2)
  val$wg    <- rep(1, n)
  val$gamma <- matrix(0.1, nrow=p, ncol=1)
  
  return(val)
}


# computes the weight of each data point based on
# estimated component parameters
weights <- function(y=NULL, x=NULL, gpar=NULL) {
  
  G <- length(gpar$pi)	
  n <- length(y)
  
  if (G > 1) {
    zlog <- matrix(0, nrow=n, ncol=length(gpar$pi))
    for (g in 1:G) {
      zlog[, g] <- dghyp(y=y, x=x, par=gpar[[g]], logd=TRUE)
      #print('dghyp')
    }
    
    w <- t(apply(zlog, 1, function(z, wt) { 
      x <- exp(z + log(wt));
      
      if (sum(x) == 0) {
        x <- rep(1, length(x))
      }
      
      x <- x/sum(x)
      return(x) 
    }, wt=gpar$pi))
  } else {
    w <- matrix(1, nrow=n, ncol=G)
  }
  
  return(w)
}

bic <- function(model) {
  loglik <- tail(model$loglik, 1)
  G      <- length(unique(model$map))
  npar   <- G * (nrow(model$gpar[[1]]$gamma) + 4) + (G - 1)
  val    <- 2 * loglik - npar * log(length(model$map))
  return(val)
}
