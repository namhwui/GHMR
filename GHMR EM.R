# Generalised Hyperbolic Mixture Regression
# Functions for: Parameter estimation via EM algorithm

# needs GH density functions
source("GHMR density.R")

# required packages
library(numDeriv)


# logged Bessel K function as a function of index parameter (s)
# needed for numerical derivatives
log_besselK_index <- function(x) {
  function(s) {
    log(besselK(x, s, expon.scaled = T)) - x
  }
}

update_lambda <- function(abc_bar, param) {
  abc_bar[3] * param$lambda / grad(log_besselK_index(param$omega), param$lambda)
}

update_omega <- function(abc_bar, param) {
  
  q_function <- function(s) {
    val <- numeric(3)
    val[1] <- (param$lambda - 1) * abc_bar[3]
    val[2] <- -s * sum(abc_bar[-3]) / 2
    val[3] <- -log(besselK(s, param$lambda))
    sum(val)
  }
  
  param$omega - grad(q_function, param$omega) / c(hessian(q_function, param$omega))

} 

update_gamma <- function(y, x, abc, wt, beta) {
  M <- diag(wt * abc[, 2])
  solve(t(x) %*% M %*% x) %*% t(x) %*% (M %*% y + beta * wt)
}

# computes aig, big, cig (conditional expectations of GIG random variable)
# component-wise
component_Estep <- function(y, x, param) {
  
  with(param, {
    mu <- x %*% gamma
    A  <- omega + (y - mu)^2 / sigma2
    B  <- omega + beta^2 / sigma2
    v  <- lambda - 0.5
    
    sqrt_AB <- sqrt(A * B)
    sqrt_AdB <- sqrt(A / B)
    bK_ratio <- besselK(sqrt_AB, v + 1) / besselK(sqrt_AB, v)
    
    abc <- matrix(nrow = length(y), ncol = 3)
    abc[, 1] <- sqrt_AdB * bK_ratio
    abc[, 2] <- -2 * v / A + bK_ratio / sqrt_AdB
    abc[, 3] <- sapply(sqrt_AB, function(x, v) {
      grad(log_besselK_index(x), v)
    }, v = v)
    abc[, 3] <- abc[, 3] + log(sqrt_AdB)
    
    abc
  })
}

component_Mstep <- function(y, x, abc, wt, param) {
  
  n <- length(wt)
  ng <- sum(wt)
  abc_bar <- colSums(sweep(abc, 1, wt, "*")) / ng
  
  r <- y - x %*% param$gamma
  param$beta <- abc_bar[2] * sum(wt * r) / sum(wt * abc[, 1])
  temp <- numeric(2)
  temp[1] <- sum(wt * abc[, 2] * r^2) / ng
  temp[2] <- -2 * param$beta * sum(wt * r) / ng 
  param$sigma2 <- sum(temp) + abc_bar[1] * param$beta^2
  param$lambda <- update_lambda(abc_bar, param)
  param$omega  <- update_omega(abc_bar, param)
  param$gamma  <- update_gamma(y, x, abc, wt, param$beta)
  
  param
}

component_EM_once <- function(y, x, wt, param) {
  abc <- component_Estep(y, x, param)
  component_Mstep(y, x, abc, wt, param)
}

MAP_label <- function(obj) {
  
  c(apply(obj$wt, 1, function(x) {
    possible_labels <- (1:length(x))[x == max(x)]
    possible_labels[1]
  }))
  
}



# Aitken's acceleration for convergence
# if converged, returns FALSE
# otherwise, returns TRUE
aitken <- function(loglik, eps = 0.01) {
  if (missing(loglik)) {
    stop("Log-likelihood values (loglik) are missing.")
  }
  if (length(loglik) < 3) {
    stop("At least 3 log-likelihood values are needed for Aitken's acceleration.")
  }
  
  three <- tail(loglik, 3)
  accel <- (three[3] - three[2]) / (three[2] - three[1])
  asymp_loglik <- three[2] + (three[3] - three[2]) / (1 - accel)
  difference   <- asymp_loglik - three[2]
  
  if (is.nan(difference)) {
    return(F)
  }
  else if (difference >= 0 & difference < eps) {
    return(F)
  } else {
    return(T)
  }
}

EM_fixed_iter <- function(obj, iter = 100) {
  
  if (missing(obj)) {
    stop("Data and parameter (obj) are missing.")
  }
  
  G_count <- 1:length(obj$prop)
  logl <- numeric(iter)
  
  for (ii in 1:iter) {
    obj$parameter <- lapply(G_count, function(g, obj) {
      component_EM_once(obj$y, obj$x, obj$wt[, g], obj$param[[g]])
    }, obj = obj) 
    
    temp     <- density_mixture(obj, F)
    obj$wt   <- sweep(temp, 1, rowSums(temp), "/")
    logl[ii] <- loglik(obj)
  }
  
  obj$loglik <- logl
  obj$label <- MAP_label(obj)
  obj$prop <- colSums(obj$wt) / length(obj$label)
  obj
}

EM_until_converge <- function(obj, eps = 0.01, max_iter = 5000) {
  
  if (missing(obj)) {
    stop("Data and parameter (obj) are missing.")
  }
  
  G_count <- 1:length(obj$prop)
  logl <- numeric(max_iter)
  
  ii <- 1
  continue <- T
  while ((ii <= max_iter) & continue)  {
    obj$parameter <- lapply(G_count, function(g, obj) {
      component_EM_once(obj$y, obj$x, obj$wt[, g], obj$param[[g]])
    }, obj = obj) 
    
    temp     <- density_mixture(obj, F)
    obj$wt   <- sweep(temp, 1, rowSums(temp), "/")
    logl[ii] <- loglik(obj)
    
    if (ii >= 3) {
      continue <- aitken(logl[1:ii], eps)
    }
    ii <- ii + 1
  }
  
  obj$loglik <- logl[1:(ii - 1)]
  obj$label <- MAP_label(obj)
  obj$prop <- colSums(obj$wt) / length(obj$label)
  obj
}


bic <- function(model) {
  loglik <- tail(model$loglik, 1)
  G      <- length(unique(model$map))
  npar   <- G * (nrow(model$gpar[[1]]$gamma) + 4) + (G - 1)
  val    <- 2 * loglik - npar * log(length(model$map))
  return(val)
}


BIC <- function(obj) {
  
  if (missing(obj$loglik)) {
    stop("Log-likelihood value is missing. EM iterations are required for log-likelihood to be added to the model object.")
  }
  
  logl <- tail(obj$loglik, 1)
  G <- length(obj$prop)
  # number of parameters for in each component:
  # ncol(obj$x) for regression coefficients
  # 4 for sigma^2, beta, lambda, omega
  
  # prop:  G - 1 parameters
  # total: G * (ncol(obj$x) + 4) + G - 1
  num_param <- G * (ncol(obj$x) + 4) + G - 1
  
  # BIC formula based on Schwarz (maximisation)
  # 2 * log-likelihood - number of free parameters * log(n)
  2 * logl - num_param * log(length(obj$y))
}
