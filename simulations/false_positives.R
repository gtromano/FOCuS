# this script investigates the average run length to a false positive

source("simulations/set_simulations.R")

output_file = "./simulations/results/fp2.RData"

sim_grid <- expand.grid(
  N = 1e5,
  changepoint = 1e5,
  delta = 0,
  threshold = c(1:13, seq(15, 50, by = 10)),
  grid_points = c(100)
)

if (T) {
  NREP <- 10
  outDF <- mclapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP))
  }, mc.cores = 6)

  outDF <- Reduce(rbind, outDF)
}

outDF

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    det_delay = ifelse(est - real > 0, stop_time - real, NA),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (is.na(det_delay)), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0) # if it's not a missed detection nor it's a false alarm, then it's a true positive
  )

summary_df


cbPalette <- c("#0072B2", "#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7")
false_alarm_plot <- ggplot(summary_df,
       aes(x = threshold, y = false_alarm, color = algo, group = algo, by = threshold)) +
  stat_summary_bin(fun.data = "mean_se", geom = "line") +
  stat_summary_bin(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = .1, lty = 2, colour = "grey") +
  scale_color_manual(values = cbPalette) +
  ylab("False Alarm Rate") +
  xlab("beta") +
  theme_idris()
false_alarm_plot



cbPalette <- c("#0072B2", "#56B4E9", "#009E73", "#33cc00", "#E69F00", "#CC79A7")
avg_run_len_plot <- ggplot(summary_df,
       aes(x = threshold, y = run_len, color = algo, group = algo, by = threshold)) +
  stat_summary_bin(fun.data = "mean_se", geom = "line") +
  stat_summary_bin(fun.data = "mean_se", geom = "errorbar") +
  geom_hline(yintercept = .1, lty = 2, colour = "grey") +
  scale_color_manual(values = cbPalette) +
  ylab("False Alarm Rate") +
  xlab("beta") +
  theme_idris()
avg_run_len_plot
