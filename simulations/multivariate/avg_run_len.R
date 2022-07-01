
source("simulations/multivariate/helper_functions.R")

set.seed(42)
Y <- generate_sequence(n = 2000, cp = 1500, magnitude = 0.5, dens = .05)


image(t(Y_test))
# mu0 estimation phase

mu0hat <- apply(Y_train, 1, mean)
res_focus0 <- FOCuS(Y_test, 2, a = .9, mu0 = rep(0, 100))
res_focus0$t

#ocd_det <- ocd_training(Y_train, thresh = setNames(c(20, 179.48, 54.87), c('diag', 'off_d', 'off_s')))

get_ocd_thres <- function(gamma, p) {
  psi <- function(t) p - 1 + t + sqrt(2 * (p-1) * t)
  diag <- log(24 * p * gamma * log(4 * p, 2))
  off_d <- psi(2 * log(24 * p * gamma * log(2 * p, 2)))
  off_s <- 8 * log(24 * p * gamma * log(2 * p, 2))
  setNames(c(diag, off_d, off_s), c('diag', 'off_d', 'off_s'))
}

ocd_det <- ocd_known(get_ocd_thres(1000, 100), rep(0, 100), rep(1, 100))
res_ocd <- ocd_detecting(Y_test, ocd_det)
res_ocd$t


##################################################
############ training false positives ############
##################################################

set.seed(1)
Y <- lapply(1:100, function(i) generate_sequence(n = 2000, cp = 500, magnitude = 0, dens = 0))

monte_carlo_focus <- function(Y) {
  focus_MC_res <- mclapply(Y, function(y) {
    res_focus0 <- FOCuS(y, Inf, a = 1, mu0 = rep(0, 100))
    sortd <- apply(res_focus0$maxs, 2, sort, decreasing = T)
    cums <- apply(sortd, 2, cumsum)
    apply(cums, 1, max)
  }, mc.cores = 8)
  Reduce(rbind, focus_MC_res)
}

focus_thres <- monte_carlo_focus(Y)
thres <- apply(focus_thres, 2, quantile, probs = .93)


set.seed(42)
Y <- lapply(1:100, function(i) generate_sequence(n = 2000, cp = 500, magnitude = 0, dens = 0))
focus_res <- mclapply(Y, function(y) {
  res_focus0 <- FOCuS(y, 5, a = 1, MC_thres = thres, mu0 = rep(0, 100))
  ifelse(res_focus0$t == -1, 2000, res_focus0$t)
}, mc.cores = 8)

mean(unlist(focus_res))
which(unlist(focus_res) < 500)



MC_ocd_v2 <- function (dim, patience, beta, sparsity, MC_reps) 
{
  peak_stat <- matrix(0, MC_reps, 3)
  colnames(peak_stat) <- c("diag", "off_d", "off_s")
  if (sparsity == "sparse") 
    peak_stat <- peak_stat[, -2]
  if (sparsity == "dense") 
    peak_stat <- peak_stat[, -3]
  for (rep in 1:MC_reps) {
    print(rep)
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
  thresh_est <- function(v) quantile(sort(v), exp(-1))
  th_individual <- apply(peak_stat, 2, thresh_est)
  th_multiplier <- thresh_est(apply(t(peak_stat)/th_individual, 
                                    2, max))
  th <- th_individual * th_multiplier
  names(th) <- colnames(peak_stat)
  return(th)
}



ocd_thres <- MC_ocd_v2(100, 2000, 1, "auto", 100)

ocd_res <- mclapply(Y, function(y) {
  ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
  res_ocd <- ocd_detecting(y, ocd_det)
  res_ocd$t
}, mc.cores = 8)

mean(unlist(ocd_res))
which(unlist(ocd_res) < 500)



#save(thres, file = "simulations/multivariate/foc_thres.RData")

load("simulations/multivariate/foc_thres.RData")
################### sparse change ########################

set.seed(42)
Y <- lapply(1:100, function(i) generate_sequence(n = 2000, cp = 500, magnitude = .5, dens = .05))

focus_res <- mclapply(Y, function(y) {
  res_focus0 <- FOCuS(y, 5, a = 1, MC_thres = thres, mu0 = rep(0, 100))
  ifelse(res_focus0$t == -1, 2000, res_focus0$t)
}, mc.cores = 8)
ddfoc <- sapply(focus_res, function(r) ifelse(r-500 > 0, r-500, NA))


ocd_thres <- c(11.6, 179.5, 54.9)
ocd_res <- mclapply(Y, function(y) {
  ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
  res_ocd <- ocd_detecting(y, ocd_det)
  res_ocd$t
}, mc.cores = 8)
ddocd <- sapply(ocd_res, function(r) ifelse(r-500 > 0, r-500, NA))

mean(log(ddfoc / ddocd), na.rm = T)

