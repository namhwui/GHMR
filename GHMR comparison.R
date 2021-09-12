# Generalised Hyperbolic Mixture Regression
# Functions for: model comparison, other than BIC

# required libraries
library(e1071)

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
