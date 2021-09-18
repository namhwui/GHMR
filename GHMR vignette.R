# GHMR vignette

set.seed(50)
source("GHMR init.R")
source("GHMR density.R")
source("GHMR EM.R")
source("GHMR comparison.R")

data(cars)
y <- cars$speed
x <- cars$dist

# there are built-in error messages during initialization,
# so please don't panick if you see an error.
# it's most likely due to having too many components at initialization,
# so that component count won't be selected.
# it's kinda slow currently because it involve several error checkpoints 
# and numerical approximations.
# we plan on speeding it up somehow before releasing it as a package.
model <- EM(y, x, G_range = 1:6, max_iter = 1000)

# you can let the G_range be a single number as well.
model <- EM(y, x, G_range = 2, max_iter = 1000)

# model selection criterion value
model$criterion

# adjusted rand index
# if the data set has a label set
#RandIndex()

# plot the clustered data
plot(x, y, col = model$label, pch = 16, main = "vignette with cars")

# model parameters
model$parameter
model$prop

