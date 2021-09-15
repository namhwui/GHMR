# Generalised Hyperbolic Mixture Regression
# Functions for: component-wise and mixture density
library(Bessel)

# computes density for univariate GHD
density_GH <- function(y, x, param, logged = F, summed = F) {
  
  if (missing(param)) {
    stop("Parameters and data (param) are missing.")
  }
  
  val <- with(param, {
    mu <- x %*% gamma
    A <- omega + (y - mu)^2 / sigma2
    B <- omega + beta^2 / sigma2
    C <- exp((y - mu) * beta / sigma2)
    
    log_use_data <- matrix(nrow = length(y), ncol = 3)
    log_use_data[, 1] <- (lambda - 0.5) / 2 * (log(A) - log(B))
    log_use_data[, 2] <- Bessel::besselK(sqrt(A * B), nu = lambda - 0.5, expon.scaled = T)
    log_use_data[, 2] <- log(log_use_data[, 2]) - sqrt(A * B)
    log_use_data[, 3] <- (mu - y) * beta / sigma2
    
    # error check for quantities involving y
    if (any(is.na(log_use_data)) | any(is.nan(log_use_data))) {
      stop("A quantity involving response variable returned NA or NaN.")
    }
    
    log_no_data <- numeric(2)
    log_no_data[1] <- log(sqrt(2 * pi * sigma2))
    log_no_data[2] <- Bessel::besselK(omega, nu = lambda, expon.scaled = T)
    log_no_data[2] <- log(log_no_data[2]) - omega
    
    apply(log_use_data, 1, sum) - sum(log_no_data)
  })
  
  if (!logged) {
    val <- exp(val)
  }
  if (summed) {
    return(sum(val))
  }
  val
}


# computes mixture density
density_mixture <- function(obj, summed = T) {
  
  # number of components
  G <- length(obj$prop)
  
  # each column belongs to one component
  val <- sapply(1:G, function(g) {
    with(obj, prop[g] * density_GH(y, x, parameter[[g]]))
  })
  
  # sum across row to get mixture density
  if (summed) {
    val <- apply(val, 1, sum)
  }
  
  val
}

# computes log-likelihood of observed data
loglik <- function(obj) {
  sum(log(density_mixture(obj)))
}
