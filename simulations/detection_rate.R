source("simulations/set_simulations.R")


output_file = "./simulations/results/dr5.RData"

sim_grid <- expand.grid(
  N = 2e5,
  changepoint = 1e5,
  delta = seq(.05, .6, by = 0.05),
  threshold = c(16)
)


if (F) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP))
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


grouped = summary_df %>% group_by(algo, magnitude) %>%
  summarise(tp_rate = mean(true_positive), det_del = mean(det_delay, na.rm = T))
print(grouped, n = 50)


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

detection_diff <- (summary_df %>% filter(algo == "FOCuS 5"))$est - (summary_df %>% filter(algo == "Page-CUSUM 50"))$est

summary2 <- summary_df %>% filter(algo == "FOCuS 5") %>%
  mutate(diff = detection_diff)


 ggplot(summary2 %>% filter(true_positive == 1),
                           aes(x = magnitude, y = diff)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
   geom_hline(yintercept = 0, lty = 2, colour = "grey") +
  xlab("magnitude") +
  ylab("Detection delay between FOCuS 5 and Page-CUSUM") +
  theme_idris()
