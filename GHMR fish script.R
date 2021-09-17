# GHMR Fish Market data analysis

source("GHMR init.R")
source("GHMR density.R")
source("GHMR EM.R")
source("GHMR comparison.R")

set.seed(50)

fish <- read.csv("Fish.csv")
data <- scale(fish[, -1])
label <- fish$ï..Species
y <- data[, 1]
x <- data[, -1]
names(data) <- colnames(fish)[-1]


ghmr_list <- lapply(1:50, function(ii) {
  EM(y, x, G_range = 1:6, sel_iter = 50, max_iter = 800, centre = F)
})

gmm_list <- lapply(1:50, function(ii) {
  regmix_over_G(y, as.matrix(x), G_range = 2:6)
})

rgmm_list <- lapply(1:50, function(ii) {
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
