# GHMR
# simulation 4: 8 dimensions, 3 components, skew normal error
# data generation

library(sn)

# set n (sample size)
#n <- 50

ng <- round(n * c(0.3, 0.3, 0.4))
label <- c(rep(1, ng[1]), rep(2, ng[2]), rep(3, ng[3]))
x <- list(matrix(runif(round(0.3*n)*7, -2, 2), ncol = 7), 
          matrix(runif(round(0.3*n)*7, -2, 2), ncol = 7) , 
          matrix(runif(round(0.4*n)*7, -2, 2), ncol = 7))
err <- list(rsn(round(0.3*n), dp = c(0, 1, 5)),
            rsn(round(0.3*n), dp = c(0, 2, -2)),
            rsn(round(0.4*n), dp = c(0, 1.5, -10)))
gamma <- rbind(c(1, 2, 1, 0, 0, 0, 0, 0),
               c(0, 0, -1, 3, 2, 0, 0, 0),
               c(0, 0, 0, 0, 0, 1, 0, -2))
y <- lapply(1:3,
            function(g) {
              cbind(1, x[[g]]) %*% gamma[g, ]  + err[[g]]
            })
data <- cbind(unlist(y), do.call("rbind", x))
colnames(data) <- c("y", paste0("x", 1:7))
