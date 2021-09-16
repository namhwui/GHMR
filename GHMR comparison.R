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


combine_d <- function(gamma1, gamma2) {
  angle <- sum(gamma1 * gamma2) / sqrt(sum(gamma1^2) * sum(gamma2^2))
  distance <- sqrt(sum((gamma1 - gamma2)^2))
  angle + distance
}


combined_label <- function(model) {
  gamma_model <- t(sapply(1:length(model$prop), function(g) {
    model$parameter[[g]]$gamma
  }))
  mat <- as.matrix(proxy::dist(gamma_model, method = combine_d))
  mat[col(mat) >= row(mat)] <- Inf
  to_combine <- which(mat == min(mat), arr.ind = T)
  label0 <- model$label
  label0[model$label %in% to_combine] <- to_combine[1]
  if (min(label0) > 1) {
    label0 <- label0 - 1
  }
  return(label0)
}


combine_components <- function(model, criterion) {
  if (length(model$prop) == 1) {
    message("Combining not applicable: There is only one component.")
    return(model)
  }
  
  if (length(model$prop) == 2) {
    model1 <- EM(model$y, model$x, G_range = 1, max_iter = 1000, centre = F)
    if (model$criterion > model1$criterion) {
      return(model)
    } else {
      return(model1)
    }
  }
  
  continue <- T
  criterion_new <- model$criterion
  criterion_old <- model$criterion
  while (continue) {
    new_label <- combined_label(model)
    new_model <- EM(model$y, model$x[, -1], label = new_label, max_iter = 1000, centre = F)
    criterion_new <- new_model$criterion
    
    if (length(new_model$prop) > 2) {
      if (criterion_new > criterion_old) {
        model <- new_model
      } else {
        continue <- F
        return(model)
      }
    } else {
      combine_components(model, criterion)
    }
  }
}


