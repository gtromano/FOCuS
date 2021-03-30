# this script investigates the average run length to a false positive

source("simulations/helper_functions.R")

run_simulation <- function(p, REPS, seed = 42) {

  grid <- find_grid(0, 50, .01, 1.3)

  set.seed(seed)
  means <- runif(REPS, 1, 10)
  data <- lapply(1:REPS, function (k) rnorm(p$N, mean = means[k]))
  m <- 100


  print("Running FOCuS")
  # FOCuS with no pruning costraint
  res <- mclapply(data, function (y) FOCuS_offline(y, Inf, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  res_FOCuS <- data.frame(sim = 1:REPS, algo = "FOCuS", est = cp, max1e3=max1e3,max1e4=max1e4, real = p$changepoint, N = p$N, threshold = p$threshold)


  print("Running FOCuS pre-change")
  # FoCUS 10
  res <- mclapply(data, function (y) FOCuS_melk(y, Inf, mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:(1e3 + m)]))
  max1e4 <- sapply(res, function (r) max(r$maxs))
  res_FoCuSpc <- data.frame(sim = 1:REPS, algo = "FOCuS est.",  est = cp, max1e3=max1e3,max1e4=max1e4, real = p$changepoint, N = p$N, threshold = p$threshold)


  print("Running yu method")
  res <- mclapply(data, function (y) yuCUSUM_v3(y, Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$cp)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs))
  res_yuCUSUM <- data.frame(sim = 1:REPS, algo = "Yu-CUSUM",  est = cp, max1e3=max1e3,max1e4=max1e4, real = p$changepoint, N = p$N, threshold = p$threshold)


  return(rbind(res_FOCuS, res_FoCuSpc, res_yuCUSUM))
}


output_file = "./simulations/pre-change-unkown/results/avgl3.RData"

sim_grid <- expand.grid(
  N = 1e4,
  changepoint = -1,
  threshold = Inf
)

CORES <- 16


if (T) {
  NREP <- 10
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)

}

load(output_file)


summary_df <-
  outDF[, c(2, 4:7)] %>% pivot_longer(
    -algo,
    names_to = "run_len",
    names_prefix = "max",
    names_transform = list(run_len = as.integer),
    values_to = "threshold"
  )


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
avg_run_len <- ggplot(summary_df %>% filter(algo != "MOSUM"),
       aes(x = run_len, y = threshold, group = run_len)) +
  stat_boxplot() +
  facet_grid(algo~., scales = "free") +
  scale_color_manual(values = cbPalette) +
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
save(thresholds, file = "simulations/thresholds-unkown.RData")
