# Generalised Hyperbolic Mixture Regression
# Functions for: model comparison, other than BIC

# required libraries
library(e1071)
library(mixtools)
library(RobMixReg)
library(flexmix)

MAP <- function(wt) {
  c(apply(wt, 1, function(x) {
    possible_labels <- (1:length(x))[x == max(x)]
    possible_labels[1]
  }))
  
}

BIC_regmix <- function(model) {
  2 * tail(model$loglik, 1) - log(length(model$y)) * (prod(dim(model$beta)) + length(model$sigma))
}


RandIndex <- function(lab1, lab2, adjusted = T) {
  val <- classAgreement(table(lab1, lab2))
  if (adjusted) {
    return(val$crand)
  }
  val$rand
}

Dist <- function(gamma_model, gamma_true) {
  if (missing(gamma_true)) {
    stop("The true regression coefficient (gamma_true) are missing.")
  }
  G <- nrow(gamma_true)
  #gamma_model <- t(sapply(model$parameter, function(x) x$gamma))
  val <- as.matrix(dist(rbind(gamma_true, gamma_model)))
  sum(val[(G + 1):nrow(val), 1:G]) / G
}


regmix_over_G <- function(y, x, G_range = 1:8, criterion = BIC_regmix, sel_iter = 50, max_iter = 2000, max_attempt = 5) {
  criterion_over_G <- sapply(G_range, function(g) {
    model <- regmixEM(y, x, k = g, maxit = sel_iter)
    criterion(model)
  })
  G <- G_range[which.max(criterion_over_G)]
  for (ii in 1:max_attempt) {
    tryCatch({
      model <- regmixEM(y, x, k = G, maxit = max_iter)
      if (missing(model)) {
        stop()
      }
      break
    },
    error = function(cond) {
      message(paste0("Re-initialising due to error. Attempts used: ", ii, "/", max_attempt))
    }, 
    finally = {})
  }
  model$criterion <- criterion(model)
  model
}


robmixreg_bic <- function(compcoef, data) {
  
  individual_density <- function(x, y) {
    z <- sapply(1:ncol(compcoef), function(g) {
      mu <- sum(cbind(1,x) * head(compcoef[, g], -2))
      return(tail(compcoef[, g], 1) * dnorm(y, mean = mu, sd = tail(compcoef[, g], 2)[1]))
    })
    return(sum(z))
  }
  
  val <- sapply(1:nrow(data), function(ii) {
    individual_density(data[ii, -1], data[ii, 1])
  })
  
  num_par <- prod(dim(compcoef)) - 1
  
  2 * sum(log(val)) - log(nrow(data)) * num_par 
}
