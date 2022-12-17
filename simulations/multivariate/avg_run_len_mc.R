source("simulations/multivariate/helper_functions.R")

CORES <- 16

#####################################################
############ training average run length ############
#####################################################

N <- 1e4

target_arl <- 4000

# data with no change
Y_nc <- lapply(1:100, function(i) generate_sequence(n = N, cp = 500, magnitude = 0, dens = 0, seed = i))                # to evaluate the average run length
Y_train <- lapply(1:100, function(i) generate_sequence(n = 500, cp = 199, magnitude = 0, dens = 0, seed = 600 + i))     # to train the mu0 value
Y_monte_carlo <- lapply(1:100, function(i) generate_sequence(n = target_arl + 100, cp = 500, magnitude = 0, dens = 0, seed = i))       # to train the monte carlo treshold


##################################
###### FOCUS0 oracle #############
##################################

focus0_mc <- mclapply(Y_monte_carlo, function(y) {
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = rep(0, 100))
  rbind(res$maxs, res$sums)
}, mc.cores = CORES)

focus0_res <- mclapply(Y_nc, function(y) {
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = rep(0, 100))
  rbind(res$maxs, res$sums)
}, mc.cores = CORES)

max_stats <- lapply(focus0_mc, function (r) apply(r, 1, max)) %>%
  Reduce(f = rbind)
# this thresholds achieves an average run length of at least 4000
foc0_thres <- apply(max_stats, 2, quantile, prob = .69)

runs <- sapply(focus0_res, function (stat) {
  for(t in seq_len(ncol(stat))) {
    over_the_treshold <- which(stat[ , t] >= foc0_thres)
    if(sum(over_the_treshold > 0)) return(t)
  }
  return(t)
})
mean(runs)

# hitted <- sapply(focus0_res, function (stat) {
#   for(t in seq_len(ncol(stat))) {
#     over_the_treshold <- which(stat[ , t] >= foc0_thres)
#     if(sum(over_the_treshold > 0)) return(over_the_treshold)
#   }
#   return("none")
# })
# table(hitted)



######################################
########## FOCuS0 est ################
######################################

focus0_est_mc <- mclapply(1:100, function(i) {
  y_tr <- Y_train[[i]]
  y <- Y_monte_carlo[[i]]
  mu0hat <- apply(y_tr, 1, mean)
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = mu0hat)
  rbind(res$maxs, res$sums)
                         # applying cumulative sums
  }, mc.cores = CORES)

focus0_est_res <- mclapply(1:100, function(i) {
  y_tr <- Y_train[[i]]
  y <- Y_nc[[i]]
  mu0hat <- apply(y_tr, 1, mean)
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf), mu0 = mu0hat)
  rbind(res$maxs, res$sums)
                         # applying cumulative sums
  }, mc.cores = CORES)

max_stats <- lapply(focus0_est_mc, function (r) apply(r, 1, max)) %>%
  Reduce(f = rbind)

foc0_est_thres <- apply(max_stats, 2, quantile, prob = .61)

runs <- sapply(focus0_est_res, function (stat) {
  for(t in seq_len(ncol(stat))) {
    over_the_treshold <- which(stat[ , t] >= foc0_est_thres)
    if(sum(over_the_treshold > 0)) return(t)
  }
  return(t)
})
mean(runs)

#################################
########## FOCuS ################
#################################

focus_mc <- mclapply(1:100, function(i) {
  y <- Y_monte_carlo[[i]]
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf))
  rbind(res$maxs, res$sums)
}, mc.cores = CORES)

focus_res <- mclapply(1:100, function(i) {
  y <- Y_nc[[i]]
  res <- FOCuS_multi_JMLR(y, c(Inf, Inf))
  rbind(res$maxs, res$sums)
}, mc.cores = CORES)


max_stats <- lapply(focus_mc, function (r) apply(r, 1, max)) %>%
  Reduce(f = rbind)

foc_thres <- apply(max_stats, 2, quantile, prob = .65)

runs <- sapply(focus_res, function (stat) {
  for(t in seq_len(ncol(stat))) {
    over_the_treshold <- which(stat[ , t] >= foc_thres)
    if(sum(over_the_treshold > 0)) return(t)
  }
  return(t)
})
mean(runs)

###################################
######### ocd oracle ##############
###################################

ocd_stat <- MC_ocd_v6(Y_monte_carlo, 1, "auto")

p <- .75
avg_run_len <- 0
while (avg_run_len < target_arl) {
  ocd_thres <- apply(ocd_stat, 2, quantile, prob = p)

  names(ocd_thres) <- colnames(ocd_stat)

  res <- mclapply(1:100, function(i) {
      y <- Y_nc[[i]]
      ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
      r <- ocd_detecting(y, ocd_det)
      r$t
    }, mc.cores = CORES)
  avg_run_len <- mean(unlist(res))
  print(avg_run_len)
  p <- min(1, p + 0.1)
}

################################
######### ocd est ##############
################################

ocd_est_stat <- MC_ocd_v6(Y_monte_carlo, 1, "auto", training_data = Y_train)

p <- .7
avg_run_len <- 0
while(avg_run_len < target_arl) {
  ocd_est_thres <- apply(ocd_est_stat, 2, quantile, prob = p)
  names(ocd_est_thres) <- colnames(ocd_est_stat)

  res <- mclapply(1:100, function(i) {
      y <- Y_nc[[i]]
      y_tr <- Y_train[[i]]

      ocd_det <- ocd_training(y_tr, ocd_est_thres)
      r <- ocd_detecting(y, ocd_det)
      r$t
    }, mc.cores = CORES)
  avg_run_len <- mean(unlist(res))
  p <- min(1, p + 0.1)
  print(avg_run_len)
}

save(foc_thres, foc0_thres, foc0_est_thres, ocd_thres, ocd_est_thres, file = "simulations/multivariate/thres.RData")

load("simulations/multivariate/thres.RData")