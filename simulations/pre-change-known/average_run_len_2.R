source("simulations/helper_functions.R")

CORES <- 16

run_simulation <- function(p, REPS, seed = 42, tlist) {
  print(p)
  grid <- find_grid(0, 50, .005, 1.15)
  
  # set.seed(seed)
  # data <- mclapply(1:REPS, function (k) rnorm(p$N, mean = means[k]), mc.cores = CORES)


  # # FOCuS with no pruning constraint
  print("Running FOCuS")
  res <- mclapply(data, function (y) FOCuS_offline(y, p$threshold, mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  st <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, threshold = p$threshold, algo = "FOCuS0", est = st, real = p$changepoint, N = p$N)
  #print("FOCus done")

  # FOCuS with 10 points grid
  print("Running FOCuS10")

  res <- mclapply(data, function (y) FOCuS_offline(y, p$threshold, mu0 = 0, grid = grid[round(seq(1, 50, length.out = 10))], K = Inf), mc.cores = CORES)
  st <- sapply(res, function (r) r$t)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, threshold = p$threshold, algo = "FOCuS0-10p", est = st, real = p$changepoint, N = p$N))

  print("Running Page")
  # Page Cusum 50 points
  res <- mclapply(data, function (y) PageCUSUM_offline(y, p$threshold, mu0 = 0, grid = grid), mc.cores = CORES)
  st <- sapply(res, function (r) r$t)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, threshold = p$threshold, algo = "Page-50p", est = st, real = p$changepoint, N = p$N))

  # simple cusum
  # res <- mclapply(data, function (y) CUSUM_offline(y, p$threshold, 0), mc.cores = CORES)
  # st <- sapply(res, function (r) r$t)
  # output <- rbind(output,
  #                 data.frame(sim = 1:REPS, threshold = p$threshold, algo = "CUSUM", est = st, real = p$changepoint, N = p$N))

  return(output)
}


output_file <- "./simulations/pre-change-known/results/avgl2_1.RData"

sim_grid <- expand.grid(
  N = 5e6,
  changepoint = -1,
  threshold = seq(1, 22, by = .5)
)


if (T) {
  REPS <- 100

  set.seed(42)
  data <- mclapply(1:REPS, function (k)  rnorm(5e6), mc.cores = CORES) # generating the data outside the function

#  run_simulation(sim_grid[1, ], REPS, tlist = tlist)

  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, REPS, tlist = tlist))
  })
  
  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}


load(output_file)

summary_df <- outDF %>% mutate(stopt = if_else(est == -1, N, est))


cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 3, 4, 5, 6)]
avg_run_len_plot <- ggplot(summary_df, aes(x = threshold, y = stopt, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylab("Run Length") +
  scale_y_log10() +
  #scale_x_log10() +
  theme_idris()




