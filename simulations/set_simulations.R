source("simulations/helper_functions.R")
run_simulation <- function(p, REPS, seed = 42, diff_thres = F) {
  print(p)

  set.seed(seed)
  data <- lapply(1:REPS, function (k) c(rnorm(p$changepoint,0), rnorm(p$N - p$changepoint, p$delta)))

  # FOCuS with no pruning costraint
  res <- lapply(data, function (y) FOCuS_offline_sim(y, p$threshold, grid = NA))
  cp <- sapply(res, function (r) r$cp)
  res_FOCuS <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("FOCus done")

  # FoCUS 5
  grid <- find_grid(0, 5, .1)
  res <- lapply(data, function (y) FOCuS_offline_sim(y, p$threshold, grid = grid))
  cp <- sapply(res, function (r) r$cp)
  res_FOCuS5 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS 5", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("page-CUSUM done")

  # Page CUSUM 5
  grid <- find_grid(0, 5, .1)
  cp <- unlist(mclapply(data, function (y) pageCUSUM_offline(y, p$threshold, grid = grid), mc.cores = 6))
  res_page25 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-CUSUM 5", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("page-CUSUM done")

  # Page CUSUM 50
  grid <- find_grid(0, 50, .01)
  cp <- unlist(mclapply(data, function (y) pageCUSUM_offline(y, p$threshold, grid = grid), mc.cores = 6))
  res_page50 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-CUSUM 50", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)



  # here put methods with different thresholds
  # CUSUM
  if (diff_thres)
    p$threshold <- 650
  cp <- unlist(mclapply(data, function (y) CUSUM_offline(y, p$threshold, 0), mc.cores = 6))
  res_CUSUM <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "CUSUM", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)

  # MOSUM
  if (diff_thres)
    p$threshold <- 27
  cp <- unlist(mclapply(data, function (y) MOSUM_offline(y, p$threshold, 30, 0), mc.cores = 6))
  res_MOSUM <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "MOSUM", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)


  return(rbind(res_FOCuS, res_FOCuS5, res_page25, res_page50, res_CUSUM, res_MOSUM))
}
