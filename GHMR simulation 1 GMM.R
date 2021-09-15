source("GHMR comparison.R")

gmm1 <- lm(y ~ x, data = data.frame(data))
gmm1$criterion <- -AIC(gmm1, k = nrow(data))
gmm <- regmix_over_G(data[,1], data[,2], G_range = 2:6)

ARI <- RandIndex(label, MAP(gmm$posterior))
D <- Dist(t(gmm$beta), gamma)
criterion <- BIC_regmix(gmm)
result <- list(model = gmm,
               BIC = gmm$criterion,
               ARI = ARI,
               dist = D)


if (gmm1$criterion > criterion) {
  result$model
  result$BIC <- gmm1$criterion
  result$ARI <- RandIndex(label, rep(1, length(label)))
  result$dist <- Dist(gmm1$coefficients, gamma)
}