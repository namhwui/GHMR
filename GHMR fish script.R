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


ghmr_list <- lapply(1:50, function(ii) {
  EM(y, x, G_range = 1:6, sel_iter = 50, max_iter = 800, centre = F)
})

gmm_list <- lapply(1:50, function(ii) {
  regmix_over_G(y, as.matrix(x), G_range = 2:6)
})

rgmm_list <- lapply(1:50, function(ii) {
  flexmix(y ~ x, data = data.frame(data), k = g, model = FLXMRrobglm(), control = list(iter = 500))
})

tle_list <- lapply(1:50, function(ii) {
  rmr(lr.method = "TLE", 
      formula   = as.formula("y ~ x"),
      data      = data.frame(data),
      MaxIt     = 200, 
      nc        = g)
})