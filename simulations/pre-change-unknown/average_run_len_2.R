source("simulations/helper_functions.R")

CORES <- 16

run_simulation <- function(p, REPS, seed = 42, tlist) {
  print(p)
  grid <- find_grid(0, 50, .01, 1.3)
  
  # set.seed(seed)
  # means <- runif(REPS, 1, 10)
  # data <- mclapply(1:REPS, function (k) rnorm(p$N, mean = means[k]), mc.cores = CORES)


  # # FOCuS with no pruning constraint
  # res <- mclapply(data, function (y) FOCuS_offline(y, p$threshold, grid = NA, K = Inf), mc.cores = CORES)
  # st <- sapply(res, function (r) r$t)
  # output <- data.frame(sim = 1:REPS, threshold = p$threshold, algo = "FOCuS", est = st, real = p$changepoint, N = p$N)
  # #print("FOCus done")
  #
  # # FOCuS with estimate of mu0 (100 obs)
  # m <- 100
  # res <- mclapply(data, function (y) FOCuS_offline(y[m:length(y)], p$threshold, mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  # st <- sapply(res, function (r) r$t)
  # #st <- sapply(res, function (r) if_else(r$t != -1, r$t + m, r$t))
  # output <- rbind(output,
  #                 data.frame(sim = 1:REPS, threshold = p$threshold, algo = "FOCuS0 100", est = st, real = p$changepoint, N = p$N))
  #
  # # FOCuS with estimate of mu0 (1000 obs)
  # m <- 1000
  # res <- mclapply(data, function (y) FOCuS_offline(y[m:length(y)], p$threshold, mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  # st <- sapply(res, function (r) r$t)
  # #st <- sapply(res, function (r) if_else(r$t != -1, r$t + m, r$t))
  # output <- rbind(output,
  #                 data.frame(sim = 1:REPS, threshold = p$threshold, algo = "FOCuS0 1000", est = st, real = p$changepoint, N = p$N))
  #
  m <- 1e4
  res <- mclapply(data, function (y) FOCuS_offline(y[m:length(y)], p$threshold, mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  st <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, threshold = p$threshold, algo = "FOCuS0 10000", est = st, real = p$changepoint, N = p$N)

  m <- 1e5
  res <- mclapply(data, function (y) FOCuS_offline(y[m:length(y)], p$threshold, mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  st <- sapply(res, function (r) r$t)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, threshold = p$threshold, algo = "FOCuS0 100000", est = st, real = p$changepoint, N = p$N))

  return(output)
}


output_file <- "./simulations/pre-change-unknown/results/avgl2_2.RData"

sim_grid <- expand.grid(
  N = 5e6,
  changepoint = -1,
  threshold = seq(1, 30, by = .5)
)


if (F) {
  set.seed(42)
  means <- runif(100, 1, 10)
  data <- mclapply(1:100, function (k) rnorm(5e6, mean = means[k]), mc.cores = CORES)

  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP, tlist = tlist))
  })
  
  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}

#
# if (T) {
#   set.seed(42)
#   means <- runif(100, 1, 10)
#   data <- mclapply(1:100, function (k) rnorm(5e6, mean = means[k]), mc.cores = CORES)
#
#   NREP <- 100
#
#   # p <- sim_grid[2, ]
#   # run_simulation(p, NREP, tlist = tlist)
#
#   outDF2 <- lapply(seq_len(nrow(sim_grid)), function (i) {
#     p <- sim_grid[i, ]
#     return(run_simulation(p, NREP, tlist = tlist))
#   })
#
#   outDF2 <- Reduce(rbind, outDF2)
#
#   load("./simulations/pre-change-unknown/results/avgl2_1.RData")
#   outDF <- rbind(outDF, outDF2)
#   save(outDF, file = output_file)
# }

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




