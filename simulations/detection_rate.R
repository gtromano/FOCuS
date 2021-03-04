source("simulations/set_simulations.R")

output_file = "./simulations/results/dr_new1.RData"

sim_grid <- expand.grid(
  N = 1e5,
  changepoint = 9e4,
  delta = seq(.05, 2, by = 0.05)
)

load("simulations/thresholds.RData")

tlist <- thresholds %>% filter(run_len == 1e5) %>% group_by(algo) %>% summarise(tres = mean(threshold))
tlist <- tlist %>% column_to_rownames(var = "algo")


run_simulation <- function(p, REPS, seed = 42, tlist) {
  print(p)
  grid <- find_grid(0, 100, .005, 1.15)
  set.seed(seed)
  data <- lapply(1:REPS, function (k) c(rnorm(p$changepoint,0), rnorm(p$N - p$changepoint, p$delta)))

  # FOCuS with no pruning costraint
  res <- mclapply(data, function (y) FOCuS_melk(y, tlist["FOCuS", 1], mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  res_FOCuS <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS", est = cp, real = p$changepoint, N = p$N)
  #print("FOCus done")

  # FoCUS 10
  grid <- grid[seq(1, 100, length.out = 10)]
  res <- mclapply(data, function (y) FOCuS_melk(y,  tlist["FOCuS 10", 1], mu0 = 0, grid = grid[seq(1, 100, length.out = 10)], K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  res_FOCuS10 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS 10", est = cp, real = p$changepoint, N = p$N)
  #print("page-CUSUM done")


  # Page CUSUM 50
  res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist["Page-CUSUM 100", 1], mu0 = 0, grid = grid), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  res_page100 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-CUSUM 100", est = cp, real = p$changepoint, N = p$N)

  # here put methods with different thresholds
  # CUSUM
   res <- mclapply(data, function (y) CUSUM_offline(y, tlist["CUSUM", 1], 0), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  res_CUSUM <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "CUSUM", est = cp, real = p$changepoint, N = p$N)

  # MOSUM
  res <- mclapply(data, function (y) MOSUMwrapper(y, tlist["MOSUM", 1], thres = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$cp)
  res_MOSUM <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "MOSUM", est = cp, real = p$changepoint, N = p$N)


  return(rbind(res_FOCuS, res_FOCuS10, res_page100, res_CUSUM, res_MOSUM))
}



if (F) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP, diff_thres = T))
  })

  outDF <- Reduce(rbind, outDF)
}

save(outDF, file = output_file)

load(output_file)

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    det_delay = ifelse(est - real > 0, est - real, NA),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (is.na(det_delay)), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positiv
  )


grouped <- summary_df %>% group_by(algo, magnitude) %>%
  summarise(no_detection = mean(no_detection), false_alarm = mean(false_alarm), tp_rate = mean(true_positive), det_del = mean(det_delay, na.rm = T))
print(grouped, n = 80)


### true positive rate ####

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
tp_rate <- ggplot(summary_df, aes(x = magnitude, y = true_positive, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  facet_grid(~ threshold) +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  ylab("True Positive Rate") +
  theme_idris()

tp_rate


### detection delay ####
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
detection_delay <- ggplot(summary_df %>% filter(true_positive == 1),
                           aes(x = magnitude, y = det_delay, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  ylab("Detection Delay") +
  theme_idris()

detection_delay

# this is to show that we do always as good as if not better than the Page-CUSUM
# as proven in proposition 1

tot_dr <- ggarrange(detection_delay, detection_delay + scale_y_continuous(trans = "log10") + ylab("log detection delay"), labels = "AUTO", nrow = 2, common.legend = T, legend = "right")
ggsave("simulations/results/dr.pdf", tot_dr, width = 10, height = 6)




#########################
algo1 <- "FOCuS"
algo2 <- "Page-CUSUM 5"

detection_diff <- (summary_df %>% filter(algo == algo1))$est - (summary_df %>% filter(algo == algo2))$est

summary2 <- summary_df %>% filter(algo == algo1) %>%
  mutate(diff = detection_diff)


 ggplot(summary2 %>% filter(true_positive == 1),
                           aes(x = magnitude, y = diff)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
   geom_hline(yintercept = 0, lty = 2, colour = "grey") +
  xlab("magnitude") +
  ylab("Detection delay between FOCuS 5 and Page-CUSUM 5") +
  theme_idris()

thresholds <- summary_df
save(thresholds, file = "simulations/thresholds.RData")
