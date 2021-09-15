# GHMM regression simulation parent file
library(parallel)
library(pbapply)

cl <- makeCluster(20)
reps <- 300
n <- 50
clusterExport(cl, 'n')


result1 <- pblapply(1:reps, function(ii) {
  set.seed(ii)
  #write.table(paste("Current seed is", ii), paste0("seed", ii, ".txt"), row.names = F, col.names = F)
  source('GHMR simulation 1 data.R', local = TRUE)
  source("GHMR simulation 1 GHMR.R", local = TRUE)
  return(result)
}, cl = cl)
saveRDS(result1, "GHMR simulation 1 GHMR.rds")

result2 <- pblapply(1:reps, function(ii) {
  set.seed(ii)
  #write.table(paste("Current seed is", ii), paste0("seed", ii, ".txt"), row.names = F, col.names = F)
  source('GHMR simulation 1 data.R', local = TRUE)
  source("GHMR simulation 1 GMM.R", local = TRUE)
  return(result)
}, cl = cl)
saveRDS(result2, "GHMR simulation 1 GMM.rds")

result3 <- pblapply(1:reps, function(ii) {
  set.seed(ii)
  #write.table(paste("Current seed is", ii), paste0("seed", ii, ".txt"), row.names = F, col.names = F)
  source('GHMR simulation 1 data.R', local = TRUE)
  source("GHMR simulation 1 RGMM.R", local = TRUE)
  return(result)
}, cl = cl)
saveRDS(result3, "GHMR simulation 1 RGMM.rds")

result4 <- pblapply(1:reps, function(ii) {
  set.seed(ii)
  #write.table(paste("Current seed is", ii), paste0("seed", ii, ".txt"), row.names = F, col.names = F)
  source('GHMR simulation 1 data.R', local = TRUE)
  source("GHMR simulation 1 TLE.R", local = TRUE)
  return(result)
}, cl = cl)
saveRDS(result4, "GHMR simulation 1 TLE.rds")




