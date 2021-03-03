# this script investigates the average run length to a false positive

source("simulations/helper_functions.R")

output_file = "./simulations/results/avgl1.RData"

sim_grid <- expand.grid(
  N = 1e6,
  changepoint = -1,
  threshold = Inf
)


run_simulation <- function(p, REPS, seed = 42, diff_thres = F) {

  set.seed(seed)
  data <- lapply(1:REPS, function (k) rnorm(p$N))
  
  # FOCuS with no pruning costraint
  res <- mclapply(data, function (y) FOCuS_melk(y, Inf, mu0 = 0, grid = NA, K = Inf), mc.cores = 6)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_FOCuS <- data.frame(sim = 1:REPS, algo = "FOCuS", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("FOCus done")
  
  # FoCUS 5
  grid <- find_grid(0, 5, .1)
  res <- mclapply(data, function (y) FOCuS_melk(y, Inf, mu0 = 0, grid = NA, K = Inf), mc.cores = 6)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_FOCuS5 <- data.frame(sim = 1:REPS, algo = "FOCuS 5",  est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("page-CUSUM done")
  
  # Page CUSUM 100
  grid <- find_grid(0, 100, .01)
  res <- mclapply(data, function (y) PageCUSUM_offline(y, Inf, mu0 = 0, grid = grid), mc.cores = 6)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_page100 <- data.frame(sim = 1:REPS, algo = "Page-CUSUM 100", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  

  res <- mclapply(data, function (y) CUSUM_offline(y, Inf, 0), mc.cores = 6)
  cp <- sapply(res, function (r) r$t)
  max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  max1e6 <- sapply(res, function (r) max(r$maxs))
  res_CUSUM <- data.frame(sim = 1:REPS, algo = "CUSUM", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  
  # # MOSUM
  # res <- mclapply(data, function (y) MOSUM_offline(y, Inf, 100, 0), mc.cores = 6)
  # cp <- sapply(res, function (r) r$cp)
  # max1e3 <- sapply(res, function (r) max(r$maxs[1:1e3]))
  # max1e4 <- sapply(res, function (r) max(r$maxs[1:1e4]))
  # max1e5 <- sapply(res, function (r) max(r$maxs[1:1e5]))
  # max1e6 <- sapply(res, function (r) max(r$maxs))
  # res_MOSUM <- data.frame(sim = 1:REPS, algo = "MOSUM", est = cp, max1e3=max1e3,max1e4=max1e4,max1e5=max1e5,max1e6=max1e6, real = p$changepoint, N = p$N, threshold = p$threshold)
  #
  
  return(rbind(res_FOCuS, res_FOCuS5, res_page100, res_CUSUM))
}



if (T) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP))
  })

  outDF <- Reduce(rbind, outDF)
}

save(outDF, file = output_file)

load(output_file)


summary_df <- outDF
summary_df %>% filter(algo == "FOCuS") %>% summary
summary_df %>% filter(algo == "Page-CUSUM 100") %>% summary
summary_df %>% filter(algo == "CUSUM") %>% summary


cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
false_alarm_plot <- ggplot(summary_df,
       aes(x = threshold, y = false_alarm, color = algo, group = algo, by = threshold)) +
  stat_summary_bin(fun.data = "mean_se", geom = "line") +
  stat_summary_bin(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = .1, lty = 2, colour = "grey") +
  geom_vline(xintercept = 27, lty = 2, colour = "grey") +
  scale_color_manual(values = cbPalette) +
  ylab("False Alarm Rate") +
  xlab("threshold") +
  theme_idris()
false_alarm_plot



cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
avg_run_len_plot <- ggplot(summary_df,
       aes(x = threshold, y = run_len, color = algo, group = algo, by = threshold)) +
  stat_summary_bin(fun.data = "mean_se", geom = "line") +
  stat_summary_bin(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylab("Average Run Length") +
  xlab("threshold") +
  theme_idris()
avg_run_len_plot


tot_fp <- ggarrange(false_alarm_plot, avg_run_len_plot, labels = "AUTO", nrow = 2, common.legend = T, legend = "right")
ggsave("simulations/results/fp.pdf", tot_fp, width = 10, height = 6)
