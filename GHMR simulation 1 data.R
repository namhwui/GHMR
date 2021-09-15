# GHMR
# simulation 1: 2 dimensions, 3 components, gaussian error
# data generation

# set n (sample size)
#n <- 50

ng <- round(n * c(0.3, 0.3, 0.4))
label <- c(rep(1, ng[1]), rep(2, ng[2]), rep(3, ng[3]))
x <- list(runif(round(0.3*n), -2, 2), 
          runif(round(0.3*n), -2, 2), 
          runif(round(0.4*n), -2, 2))
err <- list(rnorm(round(0.3*n), sd = 1),
            rnorm(round(0.3*n), sd = 0.3),
            rnorm(round(0.4*n), sd = 1.2))
gamma <- rbind(c(1, 2),
               c(1, 0.5),
               c(-1, -2))
y <- lapply(1:3,
            function(g) {
              gamma[g, 1] + gamma[g, 2] * x[[g]] + err[[g]]
            })
data <- cbind(unlist(y), unlist(x))
colnames(data) <- c("y", "x")