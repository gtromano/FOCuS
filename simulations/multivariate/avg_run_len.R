source("simulations/multivariate/helper_functions.R")

CORES <- 16

#####################################################
############ training average run length ############
#####################################################

# data with no change
Y_nc <- lapply(1:100, function(i) generate_sequence(n = 5000, cp = 500, magnitude = 0, dens = 0, seed = i))
Y_train <- lapply(1:100, function(i) generate_sequence(n = 500, cp = 200, magnitude = 0, dens = 0, seed = 600 + i))

# we keep increasing the threshold until we either hit 1% false positives and avg run leng > 1000

### FOCuS0 - pre-change mean oracle ###

foc0_thres <- 4.1
increment <- .02

avg_run_len <- 0
while (avg_run_len < 2000) {

  foc0_thres <- foc0_thres + increment
  
  focus_res <- mclapply(Y_nc, function(y) {
    res_focus0 <- FOCuS(y, foc0_thres, mu0 = rep(0, 100))
    ifelse(res_focus0$t == -1, 2000, res_focus0$t)
  }, mc.cores = CORES)
  
  avg_run_len <- mean(unlist(focus_res), na.rm = T)
  print(avg_run_len)
  
  
}


### FOCuS0 - pre-change mean estimated ###

foc0_est_thres <- 30
increment <- .5

avg_run_len <- 0
while (avg_run_len < 2000) {
  foc0_est_thres <- foc0_est_thres + increment
  
  focus_res <- mclapply(1:100, function(i) {
    y_tr <- Y_train[[i]]
    y <- Y_nc[[i]]
    mu0hat <- apply(y_tr, 1, mean)
    res_focus0 <- FOCuS(y, foc0_est_thres, mu0 = mu0hat)
    ifelse(res_focus0$t == -1, 2000, res_focus0$t)
  }, mc.cores = CORES)
  
  avg_run_len <- mean(unlist(focus_res), na.rm = T)
  print(avg_run_len)
}


### FOCuS - pre change mean unknown ###

foc_thres <- 8.5
increment <- .05

while(avg_run_len < 2000) {
  foc_thres <- foc_thres + increment
  
  focus_res <- mclapply(Y_to_check, function(y) {
    res_focus <- FOCuS(y, foc_thres)  # here we do not provide any information about mu0 
    ifelse(res_focus$t == -1, 2000, res_focus$t)
  }, mc.cores = CORES)
  
  avg_run_len <- mean(unlist(focus_res), na.rm = T)
  print(avg_run_len)
  
}


### ocd - pre change mean oracle ###

# let's get an initial estimate of the threshold
ocd_thres <- MC_ocd_v2(100, 2000, 1, "auto", 50)

ocd_res <- mclapply(Y_nc, function(y) {
  ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
  res_ocd <- ocd_detecting(y, ocd_det)
  res_ocd$t
}, mc.cores = CORES)

avg_run_len <-  mean(unlist(ocd_res), na.rm = T)


# fine tuning
increment <- c(.5, 2, 1)

while (avg_run_len < 2000) {
  
  ocd_thres <- ocd_thres + increment

  ocd_res <- mclapply(Y_nc, function(y) {
    ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
    res_ocd <- ocd_detecting(y, ocd_det)
    res_ocd$t
  }, mc.cores = CORES)
  
  avg_run_len <- mean(unlist(ocd_res), na.rm = T)
  print(avg_run_len)
}


save(foc_thres, foc0_thres, foc0_est_thres, ocd_thres, file = "simulations/multivariate/thres.RData")

# 
# # a small increment for all 3 thresholds
# increment <- c(.1, .5, .3)
# 
# Y_to_check <- Y_nc
# while (length(Y_to_check) > 1) {
#   
#   ocd_thres <- ocd_thres + increment
#   
#   ocd_res <- mclapply(Y_to_check, function(y) {
#     ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
#     res_ocd <- ocd_detecting(y, ocd_det)
#     res_ocd$t
#   }, mc.cores = CORES)
#   
#   
#   fp <- which(unlist(ocd_res) < 500)
#   cat("False positives:", fp, "\n")
#   Y_to_check <- Y_to_check[fp]
#   
# }


#### ocd - pre change mean estimated ####

ocd_est_thres <-  c(89, 568, 249)

# to_check <- 1:100
# while (length(to_check) > 1) {
#   
#   ocd_est_thres <- ocd_est_thres + increment
#   
#   ocd_res <- mclapply(to_check, function(i) {
#     y_tr <- Y_train[[i]]
#     y <- Y_nc[[i]]
#     
#     ocd_det <- ocd_training(y_tr, ocd_est_thres) # this function trains ocd
#     res_ocd <- ocd_detecting(y, ocd_det)
#     res_ocd$t
# 
#   }, mc.cores = CORES)
#   
#   fp <- which(unlist(ocd_res) < 500)
#   cat("False positives:", fp, "\n")
#   to_check <- to_check[fp]
#   
# }

# we tune now up to guaranteed average run length of minimum 1000

avg_run_len <- 0

increment <- c(1, 10, 5)

while (avg_run_len < 2000) {
  
  ocd_est_thres <- ocd_est_thres + increment
  
  
  ocd_res <- mclapply(1:100, function(i) {
    y_tr <- Y_train[[i]]
    y <- Y_nc[[i]]
    
    ocd_det <- ocd_training(y_tr, ocd_est_thres) # this function trains ocd
    res_ocd <- ocd_detecting(y, ocd_det)
    res_ocd$t
    
  }, mc.cores = CORES)
  
  avg_run_len <- mean(unlist(ocd_res), na.rm = T)
  print(avg_run_len)
}



save(foc_thres, foc0_thres, foc0_est_thres, ocd_thres, ocd_est_thres, file = "simulations/multivariate/thres.RData")



load("simulations/multivariate/thres.RData")

#
# ################### sparse change ########################
#
# Y <- lapply(1:100, function(i) generate_sequence(n = 2000, cp = 500, magnitude = 1, dens = .05, seed = i))
#
# focus_res <- mclapply(Y, function(y) {
#   res_focus <- FOCuS(y, foc_thres)
#   ifelse(res_focus$t == -1, 2000, res_focus$t)
# }, mc.cores = CORES)
# ddfoc <- sapply(focus_res, function(r) ifelse(r-500 > 0, r-500, NA))
#
#
# focus0_res <- mclapply(Y, function(y) {
#   res_focus0 <- FOCuS(y, foc0_thres, mu0 = rep(0, 100))
#   ifelse(res_focus0$t == -1, 2000, res_focus0$t)
# }, mc.cores = CORES)
# ddfoc0 <- sapply(focus0_res, function(r) ifelse(r-500 > 0, r-500, NA))
#
#
# ocd_res <- mclapply(Y, function(y) {
#   ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
#   res_ocd <- ocd_detecting(y, ocd_det)
#   res_ocd$t
# }, mc.cores = CORES)
# ddocd <- sapply(ocd_res, function(r) ifelse(r-500 > 0, r-500, NA))
#
# mean(ddfoc, na.rm = T)
# mean(ddocd, na.rm = T)
#
# cbind(ddfoc, ddocd)
# mean(log(ddfoc / ddocd), na.rm = T)
