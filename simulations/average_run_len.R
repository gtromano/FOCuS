# this script investigates the average run length to a false positive

source("simulations/helper_functions.R")

output_file = "./simulations/results/avgl1.RData"

sim_grid <- expand.grid(
  N = 1e6,
  changepoint = -1,
  threshold = Inf
)


CORES <- 16

run_simulation <- function(p, REPS, seed = 42, diff_thres = F) {

  set.seed(seed)
  data <- lapply(1:REPS, function (k) rnorm(p$N))

  print("Running FOCuS")
  # FOCuS with no pruning costraint
  res <- mclapply(data, function (y) FOCuS_melk(y, Inf, mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_FOCuS <- data.frame(sim = 1:REPS, algo = "FOCuS", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("FOCus done")

  print("Running FOCuS 5")
  # FoCUS 5
  grid <- find_grid(0, 5, .1)
  res <- mclapply(data, function (y) FOCuS_melk(y, Inf, mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_FOCuS5 <- data.frame(sim = 1:REPS, algo = "FOCuS 5",  est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("page-CUSUM done")

  print("Running Page-CUSUM")
  # Page CUSUM 100
  grid <- find_grid(0, 100, .01)
  res <- mclapply(data, function (y) PageCUSUM_offline(y, Inf, mu0 = 0, grid = grid), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_page100 <- data.frame(sim = 1:REPS, algo = "Page-CUSUM 100", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  
  print("Running CUSUM")
  res <- mclapply(data, function (y) CUSUM_offline(y, Inf, 0), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_CUSUM <- data.frame(sim = 1:REPS, algo = "CUSUM", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  
  # MOSUM
  print("Running MOSUM")
  res <- mclapply(data, function (y) MOSUMwrapper(y, 50, thres = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$cp)
  max1e3 <- sapply(res, function (r) max(r$out$stat[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$out$stat[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$out$stat[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$out$stat))
  res_MOSUM <- data.frame(sim = 1:REPS, algo = "MOSUM", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  
  return(rbind(res_FOCuS, res_FOCuS5, res_page100, res_CUSUM, res_MOSUM))
}


if (T) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)

}




load(output_file)


outDF
outDF %>% filter(algo == "FOCuS") %>% summary
outDF %>% filter(algo == "Page-CUSUM 100") %>% summary
outDF %>% filter(algo == "CUSUM") %>% summary


summary_df <-
  outDF[, c(2, 4:7)] %>% pivot_longer(
    -algo,
    names_to = "run_len",
    names_prefix = "max",
    names_transform = list(run_len = as.integer),
    values_to = "threshold"
  )


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
avg_run_len <- ggplot(summary_df,
       aes(x = run_len, y = threshold, group = run_len)) +
  stat_boxplot() +
  facet_grid(algo~., scales = "free") +
  scale_color_manual(values = cbPalette) +
  scale_x_log10() +
  ylab("Threshold") +
  xlab("Average Run Length") +
  theme_idris()
avg_run_len

ggsave("simulations/results/avg_run_len.pdf", avg_run_len, width = 5, height = 8)
