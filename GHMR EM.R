# Generalised Hyperbolic Mixture Regression
# Functions for: Parameter estimation via EM algorithm

# needs GH density functions
source("GHMR density.R")

# required packages
library(numDeriv)
library(Bessel)

# logged Bessel K function as a function of index parameter (s)
# needed for numerical derivatives
log_besselK_index <- function(x) {
  function(s) {
    log(BesselK(x, s, expon.scaled = T)) - x
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
    val[3] <- -log(BesselK(s, param$lambda))
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
    
    
    bK_ratio <- BesselK(sqrt_AB, v + 1, expon.scaled = T) / BesselK(sqrt_AB, v, expon.scaled = T)
    #bK_ratio <- besselK(sqrt_AB, v + 1, expon.scaled = T) / besselK(sqrt_AB, v, expon.scaled = T)
    
    #print(any(is.nan(bK_ratio)))
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
  #print(abc[, 1])
  param$beta <- sum(wt * r) / sum(wt * abc[, 1])
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
  
  if (is.null(obj$loglik)) {
    stop("Log-likelihood value is missing. EM iterations are required for log-likelihood to be added to the model object.")
  }
  
  logl <- tail(obj$loglik, 1)
  G <- length(obj$prop)
  # number of parameters for in each component:
  # ncol(obj$x) for regression coefficients
  # 4 for sigma^2, beta, lambda, omega
  
  # prop:  G - 1 parameters
  # total: G * (ncol(obj$x) + 4) + G - 1.
  num_param <- G * (ncol(obj$x) + 4) + G - 1
  
  # BIC formula based on Schwarz (maximisation)
  # 2 * log-likelihood - number of free parameters * log(n).
  2 * logl - num_param * log(length(obj$y))
}


model_selection <- function(y, x, G_range, sel_iter, init_method, criterion, add_intercept, centre) {
  sapply(G_range, function(g, y, x, method, add_intercept, centre, sel_iter, criterion) {
    tryCatch({
      
      # using tryCatch in case of errors from numerical issues
      # caused by numDeriv.
      
      obj <- object_mixture(y, x, g, label = NULL, method, add_intercept, centre)
      obj <- EM_fixed_iter(obj, sel_iter)
      criterion(obj)
    },
    error = function(cond) {
      message("\n Criterion value of -Inf will be returned due to error. The error could be numerical, or due to a component size being too small.")
      message("Here is the original error message:")
      message(cond)
      return(-Inf)
    },
    finally = {}) 
  }, y = y, x = x, method = init_method, add_intercept = add_intercept, centre = centre, sel_iter = sel_iter, criterion = criterion)
}

EM <- function(y, x, G_range = 1:8, label = NULL, init_method = "kmeans", 
               criterion = BIC, sel_iter = 100, iter = NULL, eps = 0.01, max_iter = 5000, 
               add_intercept = T, centre = T, seed = NULL, max_attempt = 5) {
  
  
  # IF seed is given, set the seed
  if (!is.numeric(seed)) {
    set.seed(seed)
  }
  
  # IF label is NULL:
  #   initialise with object_mixture for each G.
  #   run EM_fixed_iter on all values in G_range.
  #   select the best G based on sel_method.
  # ELSE:
  #   set G = number of classes in label.
  
  # run EM_fixed_iter or EM_until_converge on best G, based on NULLness of iter.
  # report the best model object (including its BIC), and the BIC value of others from selection stage.
  selection <- NULL
  if (is.null(label)) {
    selection <- model_selection(y, x, G_range, sel_iter, init_method, criterion, add_intercept, centre)
    names(selection) <- G_range
    G <- G_range[which.max(selection)]
  } else {
    G <- length(unique(label))
  }
  
 # print(selection)
  
  for (kk in 1:max_attempt) {
    tryCatch({
      obj <- object_mixture(y, x, G, label, init_method, add_intercept, centre)
      if (!is.null(selection)) {
        obj$BIC_all_G <- selection
      }
      
      if (is.null(iter)) {
        obj <- EM_until_converge(obj, eps, max_iter)
      } else {
        obj <- EM_fixed_iter(obj, iter)
      }
      
      # minimum threshold for scale parameter to avoid degenerate estimates
      all_sigma2 <- sapply(obj$parameter, function(comp) {comp$sigma2})
      if (any(all_sigma2 < 0.0001)) {
        stop("The scale parameters are too small.")
      }
      
      break
    },
    error = function(cond) {
      message("\n An error occurred. It's likely due to degenerate estimates.")
      message(paste("Here is the original error:", cond))
      message(paste0("Attempts used: ", kk, "/", max_attempt))
    }, 
    finally = {
      #if (kk == max_attempt) {
      #  stop("Max number of attempts reached with no convergence.")
      #}
    })
  }
  obj$BIC <- BIC(obj)
  obj
}




