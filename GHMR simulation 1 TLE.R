source("GHMR comparison.R")

model_list <- lapply(1:6, function(g) {
  model <- rmr(lr.method = "TLE", 
               formula   = as.formula("y ~ x"),
               data      = data.frame(data),
               MaxIt     = 200, 
               nc        = g)
  return(model)
})


bic_list <- sapply(model_list, function(model) {
  robmixreg_bic(model@compcoef, data)
})


tle <- model_list[[which.max(bic_list)]]
gamma_model <- t(tle@compcoef[1:2, ])

result <- list(model = tle, 
               BIC = max(bic_list), 
               ARI = RandIndex(label, tle@ctleclusters), 
               dist = Dist(gamma_model = gamma_model, gamma))
