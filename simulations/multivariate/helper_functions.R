library(ocd)
library(FOCuS)
library(parallel)

generate_sequence <- function(n = 1000, p = 100, cp = 200, sd = 1, magnitude = 1, dens = .1, seed = 42){
  
  set.seed(seed)
  noise <- matrix(rnorm(p * n, sd = sd), nr = p, nc = n)

  pre_change <- rep(0, p)#rnorm(p, mean = 0)
  #pre_change <- rnorm(p, mean = 0)
  
  post_change <- pre_change + c(sample(c(magnitude, -magnitude), floor(p * dens), replace = T), rep(0, floor(p * (1-dens))))
  noise + cbind(matrix(pre_change, nr = p, nc = cp), matrix(post_change, nr = p, nc = n - cp))

}


# this function takes treshold and known values for pre and post change mean and return an instance of
# an ocd detector ready to do some changepoint monitoring.
ocd_known <- function (thresh, mu0, sd0){
  detector <- ChangepointDetector(dim=length(mu0), method='ocd', beta=1, thresh=thresh)
  detector <- setBaselineMean(detector, mu0)
  detector <- setBaselineSD(detector, sd0)
  
  setStatus(detector, 'monitoring')
}


# this function takes some training observations and return an instance of
# an ocd detector ready to do some changepoint monitoring.
ocd_training <- function (Y_train, thresh){
  detector <- ChangepointDetector(dim=nrow(Y_train), method='ocd', beta=1, thresh=thresh)
  detector <- setStatus(detector, 'estimating')
  for (t in seq_len(ncol(Y_train))) {
    detector <- getData(detector, Y_train[,t])
  }
  setStatus(detector, 'monitoring')
}

ocd_detecting <- function (Y, detector) {
  for (t in seq_len(ncol(Y))){
    detector <- getData(detector, Y[,t])
    if (is.numeric(attr(detector, "status")))
      break
  }
  list(t = t, det = detector)
}


MC_ocd_v2 <- function (dim, patience, beta, sparsity, MC_reps) 
{
  peak_stat <- matrix(0, MC_reps, 3)
  colnames(peak_stat) <- c("diag", "off_d", "off_s")
  if (sparsity == "sparse") 
    peak_stat <- peak_stat[, -2]
  if (sparsity == "dense") 
    peak_stat <- peak_stat[, -3]
  for (rep in 1:MC_reps) {
    cat(rep, " ")
    A <- matrix(0, dim, 1)
    tail <- matrix(0, dim, floor(log2(dim)) * 2 + 4)
    for (i in 1:patience) {
      x_new <- rnorm(dim)
      ret <- ocd_update(x_new, A, tail, beta, sparsity)
      A <- ret$A
      tail <- ret$tail
      peak_stat[rep, ] <- pmax(peak_stat[rep, ], ret$stat)
    }
  }
  cat("\n")
  thresh_est <- function(v) quantile(sort(v), exp(-1))
  th_individual <- apply(peak_stat, 2, thresh_est)
  th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 
                                    2, max))
  th <- th_individual * th_multiplier
  names(th) <- colnames(peak_stat)
  return(th)
}


# takes a list of matrices rather then generating the observations at random
MC_ocd_v3 <- function (Y, beta, sparsity) 
{
  peak_stat <- matrix(0, length(Y), 3)
  colnames(peak_stat) <- c("diag", "off_d", "off_s")
  if (sparsity == "sparse") 
    peak_stat <- peak_stat[, -2]
  if (sparsity == "dense") 
    peak_stat <- peak_stat[, -3]
  for (rep in 1:length(Y)) {
    cat(rep, " ")
    y <- Y[[rep]]
    dim <- nrow(y)
    A <- matrix(0, dim, 1)
    tail <- matrix(0, dim, floor(log2(dim)) * 2 + 4)
    for (i in 1:ncol(y)) {
      x_new <- y[, i]
      ret <- ocd_update(x_new, A, tail, beta, sparsity)
      A <- ret$A
      tail <- ret$tail
      peak_stat[rep, ] <- pmax(peak_stat[rep, ], ret$stat)
    }
  }
  cat("\n")
  thresh_est <- function(v) quantile(sort(v), exp(-1))
  th_individual <- apply(peak_stat, 2, thresh_est)
  th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 
                                    2, max))
  th <- th_individual * th_multiplier
  names(th) <- colnames(peak_stat)
  return(th)
}

MC_ocd_v4 <- function (Y, beta, sparsity, training_data = NA, CORES = 16) 
{

  peak_stat <- mclapply(1:length(Y), function(rep) {
    cat(rep, " ")
    ps <- c(0, 0, 0)
    y <- Y[[rep]]
    dim <- nrow(y)
    A <- matrix(0, dim, 1)
    
    if(is.na(training_data[1])){
      mu0 <- rep(0, dim)
    } else {
      mu0 <- apply(training_data[[rep]], 1, mean)
    }
    
    tail <- matrix(0, dim, floor(log2(dim)) * 2 + 4)
    for (i in 1:ncol(y)) {
      x_new <- y[, i] - mu0
      ret <- ocd_update(x_new, A, tail, beta, sparsity)
      A <- ret$A
      tail <- ret$tail
      ps <- pmax(ps, ret$stat)
    }
    return(ps)
  }, mc.cores = CORES)
  
  peak_stat <- Reduce(rbind, peak_stat)
  colnames(peak_stat) <- c("diag", "off_d", "off_s")
  
  cat("\n")
  thresh_est <- function(v) quantile(sort(v), exp(-1))
  th_individual <- apply(peak_stat, 2, thresh_est)
  th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 
                                    2, max))
  th <- th_individual * th_multiplier
  names(th) <- colnames(peak_stat)
  return(th)
}


MC_ocd_v5 <- function (Y, beta, sparsity, training_data = NA, CORES = 16) 
{
  
  peak_stat <- mclapply(1:length(Y), function(rep) {
    cat(rep, " ")
    ps <- c(0, 0, 0)
    y <- Y[[rep]]
    dim <- nrow(y)
    A <- matrix(0, dim, 1)
    
    if(is.na(training_data[1])){
      mu0 <- rep(0, dim)
    } else {
      mu0 <- apply(training_data[[rep]], 1, mean)
    }
    
    tail <- matrix(0, dim, floor(log2(dim)) * 2 + 4)
    for (i in 1:ncol(y)) {
      x_new <- y[, i] - mu0
      ret <- ocd_update(x_new, A, tail, beta, sparsity)
      A <- ret$A
      tail <- ret$tail
      ps <- pmax(ps, ret$stat)
    }
    return(ps)
  }, mc.cores = CORES)
  
  peak_stat <- Reduce(rbind, peak_stat)
  colnames(peak_stat) <- c("diag", "off_d", "off_s")
  
  cat("\n")
  thresh_est <- function(v) quantile(sort(v), .95)
  th_individual <- apply(peak_stat, 2, thresh_est)
  th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 
                                    2, max))
  th <- th_individual * th_multiplier
  names(th) <- colnames(peak_stat)
  return(th)
}




get_ocd_thres <- function(gamma, p) {
  psi <- function(t) p - 1 + t + sqrt(2 * (p-1) * t)
  diag <- log(24 * p * gamma * log(4 * p, 2))
  off_d <- psi(2 * log(24 * p * gamma * log(2 * p, 2)))
  off_s <- 8 * log(24 * p * gamma * log(2 * p, 2))
  setNames(c(diag, off_d, off_s), c('diag', 'off_d', 'off_s'))
}



