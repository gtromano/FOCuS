source("simulations/helper_functions.R")

CORES <- 16
SEED <- 45

run_simulation <- function(p, REPS, noise, tlist) {
  print(p)
  grid <- find_grid(0, 26, .01, 1.74)
  data <- mclapply(noise, function (epsilon) c(rep(0, p$changepoint), rep(p$delta, p$N - p$changepoint)) + epsilon, mc.cores = CORES)

  # FOCuS with no pruning costraint
  print("FOCus0")
  res <- mclapply(data, function (y) FOCuS_offline(y, tlist["FOCuS"], mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0", est = cp, real = p$changepoint, N = p$N)

  # FoCUS melk
  print("FOCus0 melk")
  res <- mclapply(data, function (y) simpleMelkman(y,  F, F), mc.cores = CORES)
  cp <- sapply(res, function (r) which(r$maxs >= tlist["FOCuSMelk"])[1])
  output <- rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0Melk", est = cp, real = p$changepoint, N = p$N))
  #print("page-CUSUM done")

  # FoCUS 10
  print("FOCus0 p10")
  res <- mclapply(data, function (y) FOCuS_offline(y,  tlist["FOCuS 10"], mu0 = 0, grid = grid[c(3, 6, 8, 11, 13, 14, 16, 19, 21, 24)], K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0-10p", est = cp, real = p$changepoint, N = p$N))
  #print("page-CUSUM done")

  # Page CUSUM 50
  print("Page 25p")
  res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist["Page-CUSUM 25"], mu0 = 0, grid = grid), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <-  rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-25p", est = cp, real = p$changepoint, N = p$N))

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




output_file <- "./simulations/pre-change-known/results/dr_new9.RData"


gg <- find_grid(0, 26, .01, 1.74)[14:24]

sim_grid <- expand.grid(
  N = 2e6,
  changepoint = 1e5,
  delta = c(.05, .07, seq(.1, 2, by = 0.1), .25, gg)
)


load("simulations/pre-change-known/results/avg_run_len_NEW2.RData")
tlist <- apply(avg_run_len, 2, function (len) row.names(avg_run_len)[which(len >= 1e6-1)][1] %>% as.numeric)



if (T) {
  NREP <- 100
  set.seed(SEED)
  noise <- lapply(1:NREP, function (i) rnorm(2e6))
  #run_simulation(sim_grid[10, ], NREP, noise, tlist = tlist)
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP, noise = noise, tlist = tlist))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}


load(output_file)

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    det_delay = ifelse(est - real > 0, est - real, NA),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (is.na(det_delay)), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positive
  )


det_del_table <- summary_df %>% filter(magnitude > 0, magnitude < 2) %>% group_by(magnitude, algo) %>% summarise(dd = mean(det_delay, na.rm = T), no_det = mean(no_detection, na.rm = T), fa = mean(false_alarm, na.rm = T))
print(det_del_table, n = 100)

pivot_wider(det_del_table[1:3], names_from = algo, values_from = dd) %>% mutate(diff = FOCuS0 - `Page-25p`) %>%  print(n = 100)

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 5, 6, 4, 3)]
detection_delay <-
  ggplot(summary_df %>% filter(true_positive == 1, algo != "FOCuS0m",  algo != "FOCuS0-10p"),
    aes(
      x = magnitude,
      y = log10(det_delay),
      group = algo,
      col = algo
    )
  ) +
  geom_vline(xintercept = gg, color = "grey") +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  ylab("Detection Delay") +
  scale_y_continuous(trans='log10') +
  xlim(0, .6) +
  theme_idris()

detection_delay

algolist <- c("FOCuS0", "Page-25p")
summary1 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = FOCuS0 - `Page-25p`,
         comparison = "FOCuS0 against Page-50p") %>%
    select(sim, magnitude, diff, comparison)

algolist <- c("FOCuS0", "FOCuS0-10p")
summary2 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = FOCuS0 - `FOCuS0-10p`,
         comparison = "FOCuS0 against FOCuS0-10p") %>%
    select(sim, magnitude, diff, comparison)

algolist <- c("FOCuS0-10p", "Page-25p")
summary3 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = `FOCuS0-10p` - `Page-25p`,
         comparison = "FOCuS0-10p against Page-50p") %>%
  select(sim, magnitude, diff, comparison)

tot_summary_diff <- rbind(summary1, summary2, summary3)

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 5, 6, 4)]
dec_diff <- ggplot(tot_summary_diff, aes(x = magnitude, y = - diff)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  facet_grid(~comparison) +
  xlab("magnitude") +
  ylab("detection advantage") +
  theme_idris()

ggsave(dec_diff, filename = "simulations/pre-change-known/results/1-dec-diff-known.pdf", width = 15, height = 5)

tot_summary_diff %>%
  group_by(comparison, magnitude) %>%
  summarise(avg = mean(diff, na.rm = T)) %>%
  filter(magnitude == .05)
