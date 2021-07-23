source("simulations/helper_functions.R")

CORES <- 16

run_simulation <- function(p, REPS, seed = 42, tlist) {
  print(p)
  #grid <- find_grid(0, 50, .01, 1.3)

  set.seed(seed)
  means <- runif(REPS, 1, 10)
  data <- lapply(1:REPS, function (k) c(rnorm(1e5 + p$changepoint, means[k]), rnorm(p$N - p$changepoint, means[k] + p$delta)))

  # FOCuS with no pruning constraint
  res <- mclapply(data, function (y) FOCuS_offline(y[(1e5 + 1):length(y)], tlist["FOCuS", 1], grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS", est = cp, real = p$changepoint, N = p$N)
  #print("FOCus done")

  # FOCuS0 100 estimate
  m <- 100
  res <- mclapply(data, function (y) FOCuS_offline(y[(1e5 + 1):length(y)], tlist[paste("FOCuS0", m), 1], mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, 
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = paste("FOCuS0", m), est = cp, real = p$changepoint, N = p$N)) 

  # FOCuS0 100 estimate
  m <- 1000
  res <- mclapply(data, function (y) FOCuS_offline(y[(1e5 + 1):length(y)], tlist[paste("FOCuS0", m), 1], mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, 
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = paste("FOCuS0", m), est = cp, real = p$changepoint, N = p$N)) 
  
  m <- 10000
  res <- mclapply(data, function (y) FOCuS_offline(y[(1e5 + 1):length(y)], tlist[paste("FOCuS0", m), 1], mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, 
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = paste("FOCuS0", m), est = cp, real = p$changepoint, N = p$N)) 
  
  m <- 100000
  res <- mclapply(data, function (y) FOCuS_offline(y[(1e5 + 1):length(y)], tlist["FOCuS0 100000", 1], mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, 
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = paste("FOCuS0", m), est = cp, real = p$changepoint, N = p$N)) 
  
  # res <- mclapply(data, function (y) yuCUSUM_v3(y, tlist["Yu-CUSUM", 1]), mc.cores = CORES)
  # cp <- sapply(res, function (r) r$cp)
  # res_yuCUSUM <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Yu-CUSUM", est = cp, real = p$changepoint, N = p$N)


  return(output)
}



output_file = "./simulations/pre-change-unknown/results/dr_ukn6.RData"

sim_grid <- expand.grid(
  N = 2e6,
  changepoint = 1e5,
  delta = seq(.05, 2, by = 0.05)
)

load("simulations/pre-change-unknown/thresholds-unknown.RData")

tlist <- thresholds %>%
  filter(run_len == 1e6) %>%
  group_by(algo) %>%
  summarise(tres = mean(threshold)) %>%
  column_to_rownames(var = "algo")

#tlist[3,1] <- 12.5

#run_simulation(sim_grid[10, ], 10, tlist = tlist) # test run
if (F) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP, tlist = tlist))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}


load(output_file)

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    det_delay = ifelse(est - real >= 0, est - real, N - real),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (det_delay == N - real), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positive
  )


grouped <- summary_df %>% group_by(algo, magnitude) %>%
  summarise(no_detection = mean(no_detection), false_alarm = mean(false_alarm), tp_rate = mean(true_positive), det_del = mean(det_delay, na.rm = T))
print(grouped, n = 80)


### false alarm rate ####

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 3, 4, 5, 6)]
fa_rate <- ggplot(summary_df %>% filter(algo != "MOSUM"), aes(x = magnitude, y = false_alarm, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
#  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylim(0, 1) +
  xlab("magnitude") +
  ylab("False Alarm Rate") +
  theme_idris()

fa_rate

tp_rate <- ggplot(summary_df %>% filter(algo != "MOSUM"), aes(x = magnitude, y = true_positive, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
#  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylim(0, 1) +
  xlab("magnitude") +
  ylab("True Positive Rate") +
  theme_idris()

tp_rate


### detection delay ####
detection_delay <-
  ggplot(
#    summary_df %>% filter(magnitude >= .3, true_positive == 1),
    summary_df %>% filter(true_positive == 1),
    aes(
      x = magnitude,
      y = det_delay,
      group = algo,
      col = algo
    )
  ) +
  stat_summary(fun.data = "mean_se", geom = "line") +
#  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  ylab("Detection Delay") +
  scale_y_log10() +
  theme_idris() + theme(legend.position = "none")



ggsave(filename = "simulations/pre-change-unknown/results/2-det-delay.pdf", ggarrange(detection_delay, legend = F), width = 6, height = 6)

#tot_dr <- ggarrange(fa_rate, tp_rate, detection_delay, labels = "AUTO", nrow = 3, common.legend = T, legend = "right")
ggsave("simulations/pre-change-unknown/results/dr.pdf", tot_dr, width = 6, height = 10)

tot_sim <- ggarrange(avg_run_len_plot, detection_delay, nrow = 1, common.legend = T, legend = "right", labels = "AUTO")
ggsave("simulations/pre-change-unknown/results/tot_sim_ukn.pdf", tot_sim, width = 15, height = 7)

#########################
# algo1 <- "FOCuS"
# algo2 <- "Page-CUSUM 5"

# detection_diff <- (summary_df %>% filter(algo == algo1))$est - (summary_df %>% filter(algo == algo2))$est
#
# summary2 <- summary_df %>% filter(algo == algo1) %>%
#   mutate(diff = detection_diff)
#
#
#  ggplot(summary2 %>% filter(true_positive == 1),
#                            aes(x = magnitude, y = diff)) +
#   stat_summary(fun.data = "mean_se", geom = "line") +
#   stat_summary(fun.data = "mean_se", geom = "errorbar") +
#   scale_color_manual(values = cbPalette) +
#    geom_hline(yintercept = 0, lty = 2, colour = "grey") +
#   xlab("magnitude") +
#   ylab("Detection delay between FOCuS 5 and Page-CUSUM 5") +
#   theme_idris()


test <- (summary_df %>% filter(algo == "FOCuS"))$det_delay - (summary_df %>% filter(algo == "FOCuS t1000"))$det_delay
mean(test, na.rm = T)

test <- (summary_df %>% filter(magnitude > .1, algo == "FOCuS"))$det_delay - (summary_df %>% filter(magnitude > .1, algo == "FOCuS t1000"))$det_delay
mean(test, na.rm = T)

test <- (summary_df %>% filter(magnitude < .1, algo == "FOCuS"))$det_delay - (summary_df %>% filter(magnitude < .1, algo == "FOCuS t1000"))$det_delay
mean(test, na.rm = T)


