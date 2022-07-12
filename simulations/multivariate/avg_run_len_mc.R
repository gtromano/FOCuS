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



compute_mvfocus_thres <- function (thres, a, nu, p) {
  P1 <- function(j) 2 * thres + (p * nu + 2 * sqrt(p * nu * thres))
  
  P2 <- function(j) 2 * thres + a *  2 * j * log(p)
  sapply(1:p, function (j) min(P1(j), P2(j)))
}



##################################
###### FOCUS0 oracle #############
##################################

focus0_mc <- mclapply(Y_nc, function(y) {
  res_focus0 <- FOCuS(y, Inf, mu0 = rep(0, 100))
  sortd <- apply(res_focus0$maxs, 2, sort, decreasing = T)   # sorting the maxs of the statistics
  apply(sortd, 2, cumsum)                            # applying cumulative sums
  #apply(cums, 1, max)
}, mc.cores = 8)

foc0_thres <- 5.36
threshold <- compute_mvfocus_thres(foc0_thres, a = .7, nu = 2.1, p = 100)
runs <- sapply(focus0_mc, function (stat) {
  for(t in seq_len(ncol(stat))) {
    over_the_treshold <- which(stat[ , t] >= threshold)
    if(sum(over_the_treshold > 0)) return(t)
  }
  return(t)
})

hitted <- sapply(focus0_mc, function (stat) {
  for(t in seq_len(ncol(stat))) {
    over_the_treshold <- which(stat[ , t] >= threshold)
    if(sum(over_the_treshold > 0)) return(ifelse(max(over_the_treshold) > 50, "dense", "sparse"))
  }
  return("none")
})
table(hitted)

(avg_run_len <- mean(runs))
increment <- .01

while (avg_run_len < target_arl) {
  foc0_thres <- foc0_thres + increment
  
  runs <- mclapply(focus0_mc, function (stat) {
    for(t in seq_len(ncol(stat))) {
      over_the_treshold <- which(stat[ , t] >= compute_mvfocus_thres(foc0_thres, a = .7, nu = 2.1, p = 100))
      if(sum(over_the_treshold > 0)) return(t)
    }
    return(t)
  }, mc.cores = CORES)
  avg_run_len <- mean(unlist(runs), na.rm = T)
  print(avg_run_len)
  print(foc0_est_thres)
}



######################################
########## FOCuS0 est ################
######################################

focus0_est_mc <- mclapply(1:100, function(i) {
  y_tr <- Y_train[[i]]
  y <- Y_nc[[i]]
  mu0hat <- apply(y_tr, 1, mean)
  res_focus0 <- FOCuS(y, Inf, mu0 = mu0hat)
  sortd <- apply(res_focus0$maxs, 2, sort, decreasing = T)   # sorting the maxs of the statistics
  apply(sortd, 2, cumsum)                            # applying cumulative sums
  }, mc.cores = CORES)

foc0_est_thres <- 96.9
increment <- .1
avg_run_len <- 0

while (avg_run_len < target_arl) {
  foc0_est_thres <- foc0_est_thres + increment
  
  runs <- mclapply(focus0_est_mc, function (stat) {
    for(t in seq_len(ncol(stat))) {
      over_the_treshold <- which(stat[ , t] >= compute_mvfocus_thres(foc0_est_thres, a = .7, nu = 2.1, p = 100))
      if(sum(over_the_treshold > 0)) return(t)
    }
    return(t)
  }, mc.cores = CORES)
  avg_run_len <- mean(unlist(runs), na.rm = T)
  print(avg_run_len)
  print(foc0_est_thres)
}


#################################
########## FOCuS ################
#################################

focus_mc <- mclapply(1:100, function(i) {
  y <- Y_nc[[i]]
  res_focus <- FOCuS(y, Inf)
  sortd <- apply(res_focus$maxs, 2, sort, decreasing = T)   # sorting the maxs of the statistics
  apply(sortd, 2, cumsum)                            # applying cumulative sums
}, mc.cores = CORES)

foc_thres <- 12.28
increment <- .01
avg_run_len <- 0

while (avg_run_len < target_arl) {
  foc_thres <- foc_thres + increment
  
  runs <- mclapply(focus_mc, function (stat) {
    for(t in seq_len(ncol(stat))) {
      over_the_treshold <- which(stat[ , t] >= compute_mvfocus_thres(foc_thres, a = .7, nu = 2.1, p = 100))
      if(sum(over_the_treshold > 0)) return(t)
    }
    return(t)
  }, mc.cores = CORES)
  avg_run_len <- mean(unlist(runs), na.rm = T)
  print(avg_run_len)
  print(foc_thres)
}




Y_monte_carlo <- lapply(1:100, function(i) generate_sequence(n = target_arl + 100, cp = 500, magnitude = 0, dens = 0, seed = i))       # to train the monte carlo treshold

###################################
######### ocd oracle ##############
###################################

ocd_thres <- MC_ocd_v5(Y_monte_carlo, 1, "auto")

################################
######### ocd est ##############
################################

ocd_est_thres <- MC_ocd_v5(Y_monte_carlo, 1, "auto", training_data = Y_train)



save(foc_thres, foc0_thres, foc0_est_thres, ocd_thres, ocd_est_thres, file = "simulations/multivariate/thres.RData")

load("simulations/multivariate/thres.RData")