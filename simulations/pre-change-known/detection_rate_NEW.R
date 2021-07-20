source("simulations/pre-change-known/set_simulations.R")

CORES <- 16


run_simulation <- function(p, REPS, noise, tlist) {
  print(p)
  grid <- find_grid(0, 50, .01, 1.3)
  data <- mclapply(noise, function (epsilon) c(rep(0, p$changepoint), rep(p$delta, p$N - p$changepoint)) + epsilon, mc.cores = CORES)

  # FOCuS with no pruning costraint
  print("FOCus0")
  res <- mclapply(data, function (y) FOCuS_melk(y, tlist["FOCuS"], mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0", est = cp, real = p$changepoint, N = p$N)

  # FoCUS 10
  print("FOCus0 p10")
  res <- mclapply(data, function (y) FOCuS_melk(y,  tlist["FOCuS 10"], mu0 = 0, grid = grid[round(seq(1, 50, length.out = 10))], K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0-10p", est = cp, real = p$changepoint, N = p$N))
  #print("page-CUSUM done")

  # Page CUSUM 50
  print("Page 50p")
  res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist["Page-CUSUM 25"], mu0 = 0, grid = grid[round(seq(1, 50, length.out = 25))]), mc.cores = CORES)
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




output_file <- "./simulations/pre-change-known/results/dr_new6.RData"


gg <- find_grid(0, 50, .01, 1.3)[round(seq(1, 50, length.out = 25))][13:25]

sim_grid <- expand.grid(
  N = 2e6,
  changepoint = 1e5,
  delta = c(seq(.05, 2, by = 0.05), gg)
)


load("simulations/pre-change-known/results/avg_run_len_NEW.RData")
tlist <- apply(avg_run_len, 2, function (len) thre_seq[which(len == 1e6)][1])



if (T) {
  NREP <- 100
  set.seed(42)
  noise <- mclapply(1:REP, function (i) rnorm(N), mc.cores = CORES)
  run_simulation(sim_grid[10, ], NREP, noise, tlist = tlist)
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




algolist <- c("FOCuS0", "Page-50p")
summary1 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = FOCuS0 - `Page-50p`,
         comparison = "FOCuS0 against Page-50p") %>%
    select(sim, magnitude, diff, comparison)

algolist <- c("FOCuS0", "FOCuS0-10p")
summary2 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = FOCuS0 - `FOCuS0-10p`,
         comparison = "FOCuS0 against FOCuS0-10p") %>%
    select(sim, magnitude, diff, comparison)

algolist <- c("FOCuS0-10p", "Page-50p")
summary3 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = `FOCuS0-10p` - `Page-50p`,
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

# grouped <- summary_df %>% group_by(algo, magnitude) %>%
#   summarise(no_detection = mean(no_detection), false_alarm = mean(false_alarm), tp_rate = mean(true_positive), det_del = mean(det_delay, na.rm = T))
# print(grouped, n = 80)


### false alarm rate ####

# cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 5, 6, 4)]
# fa_rate <- ggplot(summary_df %>% filter(algo != "MOSUM"), aes(x = magnitude, y = false_alarm, group = algo, col = algo)) +
#   stat_summary(fun.data = "mean_se", geom = "line") +
#   stat_summary(fun.data = "mean_se", geom = "errorbar") +
#   scale_color_manual(values = cbPalette) +
#   ylim(0, 1) +
#   xlab("magnitude") +
#   ylab("False Alarm Rate") +
#   theme_idris()
#
# fa_rate
#
# tp_rate <- ggplot(summary_df %>% filter(algo != "MOSUM"), aes(x = magnitude, y = true_positive, group = algo, col = algo)) +
#   stat_summary(fun.data = "mean_se", geom = "line") +
#   stat_summary(fun.data = "mean_se", geom = "errorbar") +
#   scale_color_manual(values = cbPalette) +
#   ylim(0, 1) +
#   xlab("magnitude") +
#   ylab("True Positive Rate") +
#   theme_idris()
#
# tp_rate


### detection delay ####
# cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 1, 3, 4, 5, 6)]
# detection_delay <- ggplot(summary_df %>% filter(true_positive == 1),
#                            aes(x = magnitude, y = det_delay, group = algo, col = algo)) +
#   stat_summary(fun.data = "mean_se", geom = "line") +
#   stat_summary(fun.data = "mean_se", geom = "errorbar") +
#   scale_color_manual(values = cbPalette) +
#   xlab("magnitude") +
#   scale_y_log10() +
#   ylab("Detection Delay") +
#   theme_idris()
#
# detection_delay
#
#
# tot_dr <- ggarrange(fa_rate, tp_rate, detection_delay, labels = "AUTO", nrow = 3, common.legend = T, legend = "right")
# ggsave("simulations/results/dr.pdf", tot_dr, width = 6, height = 10)


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


test <- (summary_df %>% filter(algo == "FOCuS0"))$det_delay - (summary_df %>% filter(algo == "Page-50p"))$det_delay
mean(test, na.rm = T)

test <- (summary_df %>% filter(magnitude < .1, algo == "FOCuS0"))$det_delay - (summary_df %>% filter(magnitude < .1, algo == "Page-50p"))$det_delay
mean(test, na.rm = T)

test = (summary_df %>% filter(algo == "FOCuS"))$det_delay - (summary_df %>% filter(algo == "FOCuS 10"))$det_delay
mean(test, na.rm = T)

test = (summary_df %>% filter(magnitude < .1, algo == "FOCuS"))$det_delay - (summary_df %>% filter(magnitude < .1, algo == "FOCuS 10"))$det_delay
mean(test, na.rm = T)




############################# mosum stuff


cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 5, 6, 4, 3)]
fa_rate <- ggplot(summary_df, aes(x = magnitude, y = false_alarm, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylim(0, 1) +
  xlab("magnitude") +
  ylab("False Alarm Rate") +
  theme_idris()

fa_rate

tp_rate <- ggplot(summary_df, aes(x = magnitude, y = true_positive, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylim(0, 1) +
  xlab("magnitude") +
  ylab("True Positive Rate") +
  theme_idris()

tp_rate


### detection delay ####
detection_delay <- ggplot(summary_df,
                           aes(x = magnitude, y = det_delay, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  ylab("Detection Delay") +
  theme_idris()

detection_delay








##################### pres plots ##################################à

source("simulations/pre-change-known/set_simulations.R")

CORES <- 16

# run_simulation <- function(p, REPS, seed = 42, tlist) {
#   print(p)
#   grid <- find_grid(0, 50, .01, 1.3)
#   set.seed(seed)
#   data <- lapply(1:REPS, function (k) c(rnorm(p$changepoint,0), rnorm(p$N - p$changepoint, p$delta)))
#
#   # FOCuS with no pruning costraint
#   res <- mclapply(data, function (y) FOCuS_offline(y, tlist["FOCuS", 1], mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
#   cp <- sapply(res, function (r) r$t)
#   res_FOCuS <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS", est = cp, real = p$changepoint, N = p$N)
#   #print("FOCus done")
#
#   # Page CUSUM 50
#   res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist["Page-CUSUM 50", 1], mu0 = 0, grid = grid), mc.cores = CORES)
#   cp <- sapply(res, function (r) r$t)
#   res_page50 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-CUSUM 50", est = cp, real = p$changepoint, N = p$N)
#
#
#   return(rbind(res_FOCuS, res_page50))
# }
#
#
#
# output_file = "./simulations/pre-change-known/results/dr_pres2.RData"
#
# sim_grid <- expand.grid(
#   N = 1e5,
#   changepoint = 9e4,
#   delta = seq(.05, 2, by = 0.05)
# )
#
# load("simulations/pre-change-known/thresholds.RData")
#
# tlist <- thresholds %>%
#   filter(run_len == 1e5) %>%
#   group_by(algo) %>%
#   summarise(tres = quantile(threshold, .75)) %>%
#   column_to_rownames(var = "algo")
#
# tlist[1, 1] <- 13.4
# tlist[3, 1] <- 13.4

#run_simulation(sim_grid[10, ], NREP, tlist = tlist)

# if (F) {
#   NREP <- 100
#   outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
#     p <- sim_grid[i, ]
#     return(run_simulation(p, NREP, tlist = tlist))
#   })
#
#   outDF <- Reduce(rbind, outDF)
#   save(outDF, file = output_file)
# }


load(output_file)

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    det_delay = ifelse(est - real > 0, est - real, NA),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (is.na(det_delay)), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positive
  )


summary_df %>%
  mutate(`change magnitude` = cut(magnitude, breaks = c(0, .3, .5, 2))) %>%
  group_by(`change magnitude`, algo) %>%
  summarise(`false positives rate` = mean(false_alarm, na.rm = T),
            `missed detection rate` = mean(no_detection, na.rm = T))

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2,  4, 3)]
fa_rate <- ggplot(summary_df %>% filter(algo %in% c("FOCuS", "Page-CUSUM 50")), aes(x = magnitude, y = false_alarm, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylim(0, 1) +
  xlab("magnitude") +
  ylab("False Alarm Rate") +
  theme_idris()

fa_rate

tp_rate <- ggplot(summary_df %>% filter(algo %in% c("FOCuS", "Page-CUSUM 50")), aes(x = magnitude, y = true_positive, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  ylim(0, 1) +
  xlab("magnitude") +
  ylab("True Positive Rate") +
  theme_idris()

tp_rate


### detection delay ####
detection_delay <- ggplot(summary_df %>% filter(algo %in% c("FOCuS", "Page-CUSUM 50"), true_positive == 1),
                          aes(x = magnitude, y = det_delay, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  ylab("Detection Delay") +
  theme_idris()
detection_delay + scale_y_log10()

tot_sim <- ggarrange(fa_rate, tp_rate, detection_delay, ncol = 3)



algolist <- c("FOCuS0", "Page-50p")
summary1 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = FOCuS0 - `Page-50p`,
         comparison = "FOCuS0 against Page-50p") %>%
    select(sim, magnitude, diff, comparison)


algolist <- c("FOCuS0", "FOCuS0-10p")
summary2 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = FOCuS0 - `FOCuS0-10p`,
         comparison = "FOCuS0 against FOCuS0-10p") %>%
    select(sim, magnitude, diff, comparison)



algolist <- c("FOCuS0-10p", "Page-50p")
summary3 <- summary_df %>% filter(algo %in% algolist) %>% select(sim, magnitude, real, N, algo, det_delay) %>%
  pivot_wider(names_from = "algo", values_from = "det_delay") %>%
  mutate(diff = `FOCuS0-10p` - `Page-50p`,
         comparison = "FOCuS0-10p against Page-50p") %>%
  select(sim, magnitude, diff, comparison)

tot_summary_diff <- rbind(summary1, summary2, summary3)

dec_diff <- ggplot(tot_summary_diff, aes(x = magnitude, y = - diff)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  facet_grid(~comparison) +
  xlab("magnitude") +
  ylab("detection advantage") +
  theme_idris()

dec_diff

tot_sim <- ggarrange(dec_diff, tp_rate, ncol = 2)


ggsave("simulations/pre-change-known/results/page-cusum-comp.pdf", dec_diff, width = 7, height = 7)

summary3 %>%
  mutate(`change magnitude` = cut(magnitude, breaks = c(0, .5, 1, 1.5, 2))) %>%
  group_by(`change magnitude`) %>%
  summarise(`false positives rate` = mean(diff, na.rm = T))


summary3 %>%
  mutate(`change magnitude` = cut(magnitude, breaks = c(0, .1, .3, .5, 1, 2))) %>%
  group_by(`change magnitude`) %>%
  summarise(`average difference` = mean(diff, na.rm = T)) %>% view
