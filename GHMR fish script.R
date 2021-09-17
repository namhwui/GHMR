# GHMR Fish Market data analysis

source("GHMR init.R")
source("GHMR density.R")
source("GHMR EM.R")
source("GHMR comparison.R")

set.seed(50)

fish <- read.csv("Fish.csv")
data <- scale(fish[, -1])
label <- fish$ï..Species
label[!(label %in% c("Perch", "Bream"))] <- "Other"
y <- data[, 1]
x <- data[, -1]
names(data) <- colnames(fish)[-1]


ghmr_list <- lapply(1:50, function(ii) {
  print(ii)
  EM(y, x, G_range = 1:6, sel_iter = 50, max_iter = 800, centre = F)
})

gmm_list <- lapply(1:50, function(ii) {
  print(ii)
  regmix_over_G(y, as.matrix(x), G_range = 2:6)
})

rgmm_list <- lapply(1:50, function(ii) {
  print(ii)
  val <- sapply(1:6, function(g) {
    model <- flexmix(y ~ x, data = data.frame(data), k = g, model = FLXMRrobglm(), control = list(iter = 50))
    -stats::BIC(model)
  })
  G <- which.max(val)
  flexmix(y ~ x, data = data.frame(data), k = G, model = FLXMRrobglm(), control = list(iter = 500))
})

tle_list <- lapply(1:50, function(ii) {
  print(ii)
  val <- sapply(1:6, function(g) {
    model <- rmr(lr.method = "TLE", 
                 formula   = Weight ~.,
                 data      = data.frame(data),
                 MaxIt     = 50, 
                 nc        = g)
    robmixreg_bic(model@compcoef, data)
  })
  G <- which.max(val)
  rmr(lr.method = "TLE", 
      formula   = Weight ~.,
      data      = data.frame(data),
      MaxIt     = 500, 
      nc        = G)
})

#saveRDS(ghmr_list, "GHMR fish GHMR.rds")
#saveRDS(gmm_list, "GHMR fish GMM.rds")
#saveRDS(rgmm_list, "GHMR fish RGMM.rds")
#saveRDS(tle_list, "GHMR fish TLE.rds")
ghmr_best <- which.max(sapply(ghmr_list, function(x) x$criterion))
gmm_best <- which.max(sapply(gmm_list, function(x) x$criterion))
rgmm_best <- which.max(sapply(rgmm_list, function(x) -stats::BIC(x)))
tle_best <- which.max(sapply(tle_list, function(x) robmixreg_bic(x@compcoef, data)))

c(ghmr_list[[ghmr_best]]$criterion, gmm_list[[gmm_best]]$criterion, -stats::BIC(rgmm_list[[rgmm_best]]), robmixreg_bic(tle_list[[tle_best]]@compcoef, data))
c(RandIndex(label, ghmr_list[[ghmr_best]]$label), RandIndex(label, MAP(gmm_list[[gmm_best]]$posterior)), RandIndex(label, rgmm_list[[rgmm_best]]@cluster), RandIndex(label, tle_list[[tle_best]]@ctleclusters))

proxy::dist(t(sapply(ghmr_list[[ghmr_best]]$parameter, function(x) x$gamma)), method = combine_d)


# plots
plot(density(data[, 1]),
     main = "Fish Weight",
     xlab = "Scaled weight")
boxplot(data[,1] ~ label, main = "Fish Weight by Label", ylab = "Scaled weight")
par(mfcol = c(2,2))
boxplot(data[,1] ~ ghmr_list[[ghmr_best]]$label, main = "GHMR", xlab = "Label", ylab = "Scaled weight")
boxplot(data[,1] ~ MAP(gmm_list[[gmm_best]]$posterior), main = "GMR", xlab = "Label", ylab = "Scaled weight")
boxplot(data[,1] ~ rgmm_list[[rgmm_best]]@cluster, main = "RGMR", xlab = "Label", ylab = "Scaled weight")
boxplot(data[,1] ~ tle_list[[tle_best]]@ctleclusters, main = "TLE", xlab = "Label", ylab = "Scaled weight")

table(label, ghmr_list[[ghmr_best]]$label)
table(label, MAP(gmm_list[[gmm_best]]$posterior))
table(label, rgmm_list[[rgmm_best]]@cluster)
table(label, tle_list[[tle_best]]@ctleclusters)


