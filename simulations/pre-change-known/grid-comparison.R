source("simulations/pre-change-known/set_simulations.R")

CORES <- 16

run_simulation <- function(p, REPS, seed = 42, tlist) {
  print(p)
  grid <- find_grid(0, 50, .01, 1.3)
  set.seed(seed)
  data <- lapply(1:REPS, function (k) c(rnorm(p$changepoint,0), rnorm(p$N - p$changepoint, p$delta)))

  # FOCuS with no pruning costraint
  print("FOCus0")
  res <- mclapply(data, function (y) FOCuS_melk(y, tlist["FOCuS", 1], mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0", est = cp, real = p$changepoint, N = p$N)

  # FoCUS 10
  print("FOCus0 p10")
  res <- mclapply(data, function (y) FOCuS_melk(y,  tlist["FOCuS 10", 1], mu0 = 0, grid = grid[round(seq(1, 50, length.out = 10))], K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0-10p", est = cp, real = p$changepoint, N = p$N))
  #print("page-CUSUM done")

  # Page CUSUM 50
  print("Page 50p")
  res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist["Page-CUSUM 50", 1], mu0 = 0, grid = grid), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <-  rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-50p", est = cp, real = p$changepoint, N = p$N))

  # here put methods with different thresholds
  # CUSUM
  # print("CUSUM")
  # res <- mclapply(data, function (y) CUSUM_offline(y, tlist["CUSUM", 1], 0), mc.cores = CORES)
  # cp <- sapply(res, function (r) r$t)
  # output <- rbind(output,
  #                 data.frame(sim = 1:REPS, magnitude = p$delta, algo = "CUSUM", est = cp, real = p$changepoint, N = p$N))

  # # MOSUM
  # res <- mclapply(data, function (y) MOSUMwrapper(y, bandw = 50, thres = tlist["MOSUM", 1]), mc.cores = CORES)
  # cp <- sapply(res, function (r) r$cp)
  # res_MOSUM <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "MOSUM", est = cp, real = p$changepoint, N = p$N)


  return(output)
}



output_file = "./simulations/pre-change-known/results/dr_comp4.RData"

sim_grid <- expand.grid(
  N = 2e6,
  changepoint = 1e5,
  delta = seq(.3, .7, by = 0.01)
)

load("simulations/pre-change-known/thresholds.RData")

tlist <- thresholds %>%
  filter(run_len == 1e6) %>%
  group_by(algo) %>%
  summarise(tres = mean(threshold)) %>%
  column_to_rownames(var = "algo")

#run_simulation(sim_grid[10, ], NREP, tlist = tlist)

if (F) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP, tlist = tlist))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}
