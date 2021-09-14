# Generalised Hyperbolic Mixture Regression
# Functions for: Parameter initialisation


# check if all columns in x are numeric
all_numeric <- function(x) {
  all(apply(x, 2, is.numeric))
}

# initialise parameter for a single component
param_component <- function(y, x, member = NULL, beta = NULL, omega = NULL, lambda = NULL, gamma = NULL) {
  
  if (missing(y)) {
    stop("Response variable (y) is missing.")
  }
  if (missing(x)) {
    stop("Covariate (x) is missing.")
  }
  if (missing(member)) {
    stop("Component membership (member) is missing.")
  }
  
  n <- length(y)
  val <- list()
  
  val$beta   <- ifelse(is.null(beta), 0, beta) 
  val$omega  <- ifelse(is.null(omega), 1, omega)
  val$lambda <- ifelse(is.null(lambda), 1, lambda)
  val$sigma2 <- var(y[member])
  val$gamma  <- gamma
  
  if (is.null(gamma)) {
    val$gamma <- lm(y[member] ~ x[member, ] - 1)$coefficients
    names(val$gamma) <- NULL
  } 
  
  val
}

# initialise parameters for all components
# the returned list contains the data set as well, 
# so that only one list needs pushing through the functions.
object_mixture <- function(y, x, G, label = NULL, method = "kmeans", add_intercept = T, centre = T) {
  
  if (missing(y)) {
    stop("Response variable (y) is missing.")
  }
  if (missing(G)) {
    stop("Component count (G) is missing.")
  }
  if (missing(x)) {
    stop("Covariate matrix (x) is missing.")
  }
  
  if (is.null(ncol(x))) {
    x <- matrix(x, ncol = 1)
  } else {
    x <- as.matrix(x)
  }
  
  if (!all_numeric(x)) {
    stop("Not all covariates are numeric. Please convert non-numeric variables.")
  }
  
  if (centre) {
    x <- scale(x, scale = F)
  }
  
  lab <- label
  if (is.null(lab)) {
    if (method == "kmeans") {
      lab <- kmeans(cbind(y, x), centers = G, nstart = 5, iter.max = 5)$cluster
    } else {
      stop("Other methods for label initialisation are not implemented yet.")
    }
  }
  
  if (add_intercept) {
    x <- cbind(1, x)
  }
  
  val <- list()
  val$parameter <- list()
  val$wt <- sapply(1:G, function(g) {
    (lab == g) * 1
  })
  
  for (g in 1:G) {
    val$parameter[[g]] <- param_component(y, x, member = (lab == g))
  }
  
  val$prop <- colSums(val$wt) / length(y)
  val$y    <- y
  val$x    <- x
  val
}	




