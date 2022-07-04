source("simulations/multivariate/helper_functions.R")

CORES <- 16

##################################################
############ training false positives ############
##################################################
# 
# set.seed(1)
# Y <- lapply(1:100, function(i) generate_sequence(n = 2000, cp = 500, magnitude = 0, dens = 0))
# 
# monte_carlo_focus <- function(Y) {
#   focus_MC_res <- mclapply(Y, function(y) {
#     res_focus0 <- FOCuS(y, Inf, a = 1, mu0 = rep(0, 100))
#     sortd <- apply(res_focus0$maxs, 2, sort, decreasing = T)
#     cums <- apply(sortd, 2, cumsum)
#     apply(cums, 1, max)
#   }, mc.cores = 8)
#   Reduce(rbind, focus_MC_res)
# }
# 
# focus_thres <- monte_carlo_focus(Y)
# thres <- apply(focus_thres, 2, quantile, probs = .93)


# data with no change
Y_nc <- lapply(1:100, function(i) generate_sequence(n = 2000, cp = 500, magnitude = 0, dens = 0, seed = i))
Y_train <- lapply(1:100, function(i) generate_sequence(n = 250, cp = 200, magnitude = 0, dens = 0, seed = 600 + i))

# we keep increasing the threshold until we either hit 1% false positives and avg run leng > 1000

### FOCuS0 - pre-change mean oracle ###

foc0_thres <- 3
increment <- .03

Y_to_check <- Y_nc
while (length(Y_to_check) > 1) {

  foc0_thres <- foc0_thres + increment
  
  focus_res <- mclapply(Y_to_check, function(y) {
    res_focus0 <- FOCuS(y, foc0_thres, mu0 = rep(0, 100))
    ifelse(res_focus0$t == -1, 2000, res_focus0$t)
  }, mc.cores = CORES)
  
  fp <- which(unlist(focus_res) < 500)
  cat("False positives:", fp, "\n")
  Y_to_check <- Y_to_check[fp]
  
}


### FOCuS0 - pre-change mean estimated ###

foc0_est_thres <- 8.5
increment <- .05

to_check <- 1:100
while (length(to_check) > 1) {
  
  foc0_est_thres <- foc0_est_thres + increment
  
  focus_res <- mclapply(to_check, function(i) {
    y_tr <- Y_train[[i]]
    y <- Y_nc[[i]]
    mu0hat <- apply(y_tr, 1, mean)
    res_focus0 <- FOCuS(y, foc0_est_thres, mu0 = mu0hat)
    ifelse(res_focus0$t == -1, 2000, res_focus0$t)
  }, mc.cores = CORES)
  
  fp <- which(unlist(focus_res) < 500)
  cat("False positives:", fp, "\n")
  to_check <- to_check[fp]
  
}



### FOCuS - pre change mean unkown ###

foc_thres <- 8.5
increment <- .05

Y_to_check <- Y_nc
while (length(Y_to_check) > 1) {
  
  foc_thres <- foc_thres + increment
  
  focus_res <- mclapply(Y_to_check, function(y) {
    res_focus <- FOCuS(y, foc_thres)  # here we do not provide any information about mu0 
    ifelse(res_focus$t == -1, 2000, res_focus$t)
  }, mc.cores = CORES)
  
  fp <- which(unlist(focus_res) < 500)
  cat("False positives:", fp, "\n")
  Y_to_check <- Y_to_check[fp]
  
}



### ocd - pre change mean oracle ###

# let's get an initial estimate of the threshold
ocd_thres <- MC_ocd_v2(100, 2500, 1, "auto", 5)

# a small increment for all 3 thresholds
increment <- c(.1, .5, .3)

Y_to_check <- Y_nc
while (length(Y_to_check) > 1) {
  
  ocd_thres <- ocd_thres + increment
  
  ocd_res <- mclapply(Y_to_check, function(y) {
    ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
    res_ocd <- ocd_detecting(y, ocd_det)
    res_ocd$t
  }, mc.cores = CORES)
  
  
  fp <- which(unlist(ocd_res) < 500)
  cat("False positives:", fp, "\n")
  Y_to_check <- Y_to_check[fp]
  
}


#### ocd - pre change mean estimated ####

ocd_est_thres <-  c(19.95634, 221.05078, 75.78118)
increment <- c(.2, 1, .5)

to_check <- 1:100
while (length(to_check) > 1) {
  
  ocd_est_thres <- ocd_est_thres + increment
  
  ocd_res <- mclapply(to_check, function(i) {
    y_tr <- Y_train[[i]]
    y <- Y_nc[[i]]
    
    ocd_det <- ocd_training(y_tr, ocd_est_thres) # this function trains ocd
    res_ocd <- ocd_detecting(y, ocd_det)
    res_ocd$t

  }, mc.cores = CORES)
  
  fp <- which(unlist(ocd_res) < 500)
  cat("False positives:", fp, "\n")
  to_check <- to_check[fp]
  
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
