# Generalised Hyperbolic Mixture Regression
# Functions for: model comparison, other than BIC

# required libraries
library(e1071)
library(mixtools)

MAP <- function(wt) {
  c(apply(wt, 1, function(x) {
    possible_labels <- (1:length(x))[x == max(x)]
    possible_labels[1]
  }))
  
}

BIC_regmix <- function(model) {
  2 * tail(model$loglik, 1) - log(length(model$y)) * (prod(dim(model$beta)) + length(model_GMM$sigma))
}


RandIndex <- function(lab1, lab2, adjusted = T) {
  val <- classAgreement(table(lab1, lab2))
  if (adjusted) {
    return(val$crand)
  }
  val$rand
}

Dist <- function(model, gamma_true) {
  if (missing(gamma_true)) {
    stop("The true regression coefficient (gamma_true) are missing.")
  }
  G <- nrow(gamma_true)
  gamma_model <- t(sapply(model$parameter, function(x) x$gamma))
  val <- as.matrix(dist(rbind(gamma_true, gamma_model)))
  sum(val[(G + 1):nrow(val), 1:G])
}


regmix_over_G <- function(y, x, G_range = 1:8, criterion = BIC_regmix, sel_iter = 50, max_iter = 2000) {
  criterion_over_G <- sapply(G_range, function(g) {
    model <- regmixEM(y, x, k = g, maxit = sel_iter)
    criterion(model)
  })
  G <- G_range[which.max(criterion_over_G)]
  regmixEM(y, x, k = G, maxit = max_iter)
}
