# this script investigates the average run length to a false positive

source("simulations/set_simulations.R")

output_file = "./simulations/results/avgrlen7.RData"

sim_grid <- expand.grid(
  N = 1e5,
  changepoint = 1e5,
  delta = 0,
  threshold = c(1:10, seq(15, 35, by = 5))
  #threshold = seq(200, 1000, by = 10)
)


if (F) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP))
  })

  outDF <- Reduce(rbind, outDF)
}

#save(outDF, file = output_file)

load(output_file)

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    det_delay = ifelse(est - real > 0, est - real, NA),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (is.na(det_delay)), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positive
)

summary_df


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
