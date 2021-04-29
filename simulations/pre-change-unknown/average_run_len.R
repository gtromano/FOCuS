# this script investigates the average run length to a false positive

source("simulations/helper_functions.R")

run_simulation <- function(p, REPS, seed = 42) {

  grid <- find_grid(0, 50, .01, 1.3)

  set.seed(seed)
  means <- runif(REPS, 1, 10)
  data <- lapply(1:REPS, function (k) rnorm(p$N, mean = means[k]))


  print("Running FOCuS")
  # FOCuS with no pruning costraint
  res <- mclapply(data, function (y) FOCuS_offline(y, Inf, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  
  res_FOCuS <- data.frame(sim = 1:REPS, algo = "FOCuS", est = cp, max1e3=max1e3, max1e4=max1e4,  max1e5=max1e5, real = p$changepoint, N = p$N, threshold = p$threshold)


  print("Running FOCuS pre-change 100")
  m <- 100
  res <- mclapply(data, function (y) FOCuS_offline(y[m:p$N], Inf, mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(c(rep(0, m), r$maxs)[1:1e3]))
  max1e4 <- sapply(res, function (r) max(c(rep(0, m), r$maxs)[1:1e4]))
  max1e5 <- sapply(res, function (r) max(c(rep(0, m), r$maxs)))
  
  res_FoCuSt100 <- data.frame(sim = 1:REPS, algo = "FOCuS t100",  est = cp, max1e3=max1e3, max1e4=max1e4, max1e5=max1e5, real = p$changepoint, N = p$N, threshold = p$threshold)


  
  print("Running FOCuS pre-change 1000")
  m <- 1000
  res <- mclapply(data, function (y) FOCuS_offline(y[m:p$N], Inf, mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(c(rep(0, m), r$maxs)[1:1e3]))
  max1e4 <- sapply(res, function (r) max(c(rep(0, m), r$maxs)[1:1e4]))
  max1e5 <- sapply(res, function (r) max(c(rep(0, m), r$maxs)))
  
  res_FoCuSt1000 <- data.frame(sim = 1:REPS, algo = "FOCuS t1000",  est = cp, max1e3=max1e3, max1e4=max1e4, max1e5=max1e5, real = p$changepoint, N = p$N, threshold = p$threshold)
  
  

  return(rbind(res_FOCuS, res_FoCuSt100, res_FoCuSt1000))
}


output_file = "./simulations/pre-change-unknown/results/avgl5.RData"

sim_grid <- expand.grid(
  N = 1e5,
  changepoint = -1,
  threshold = Inf
)

CORES <- 16


if (F) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)

}

load(output_file)


summary_df <-
  outDF[, c(2, 4:6)] %>% pivot_longer(
    -algo,
    names_to = "run_len",
    names_prefix = "max",
    names_transform = list(run_len = as.integer),
    values_to = "threshold"
  )


avg_run_len <- ggplot(summary_df %>% filter(algo != "MOSUM"),
       aes(x = run_len, y = threshold, group = run_len)) +
  stat_boxplot() +
  facet_grid(algo~., scales = "free") +
  scale_x_log10() +
  ylab("Threshold") +
  xlab("Average Run Length") +
  theme_idris()
avg_run_len

#ggsave("simulations/plots/avg_run_len.pdf", avg_run_len, width = 6, height = 10)


# plot of the grids
# grid <- find_grid(0, 50, .01, 1.3)
# par(mfrow = c(1,2))
# plot(NULL, xlim = c(-6, 6), ylim = c(-1, 1), ylab = " ", xlab = expression(mu), main = "Page-CUSUM 50 points grid")
# abline(v = grid)
# plot(NULL, xlim = c(-6, 6), ylim = c(-1, 1), ylab = " ", xlab = expression(mu), main = "FOCuS 10 points grid")
# abline(v = grid[round(seq(1, 50, length.out = 10))])
# par(mfrow = c(1,1))


thresholds <- summary_df
save(thresholds, file = "simulations/pre-change-unknown/thresholds-unknown.RData")
