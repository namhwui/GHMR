source("GHMR comparison.R")

model_list <- lapply(1:6, function(g) {
  flexmix(y ~ ., data = data.frame(data), k = g, model = FLXMRrobglm(), control = list(iter = 500))
})

bic_list <- sapply(model_list, function(x) -stats::BIC(x))
rgmm <- model_list[[which.max(bic_list)]]
gamma_model <- t(sapply(rgmm@components, function(comp) {
  comp[[1]]@parameters$coef
}))

result <- list(model = rgmm,
               BIC = -stats::BIC(rgmm),
               ARI = RandIndex(label, rgmm@cluster),
               dist = Dist(gamma_model, gamma))
