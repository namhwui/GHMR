source("GHMR init.R")
source("GHMR density.R")
source("GHMR EM.R")
source("GHMR comparison.R")

ghmr <- EM_all(data[,2], data[,1], G_range = 1:6, max_iter = 1000, centre = F)

gamma_model <- t(sapply(ghmr$parameter, function(comp) {
  comp$gamma
}))

result <- list(model = ghmr,
               BIC = ghmr$criterion,
               ARI = RandIndex(label, ghmr$label),
               dist = Dist(gamma_model, gamma))
