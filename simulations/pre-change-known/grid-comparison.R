source("simulations/pre-change-known/set_simulations.R")

CORES <- 16

run_simulation <- function(p, REPS, seed = 42, tlist) {
  print(p)
  grid <- find_grid(0, 50, .01, 1.3)
  set.seed(seed)
  data <- mclapply(1:REPS, function (k) c(rnorm(p$changepoint,0), rnorm(p$N - p$changepoint, p$delta)), mc.cores = CORES)

  # FOCuS with no pruning costraint
  print("FOCus0")
  res <- mclapply(data, function (y) FOCuS_melk(y, tlist["FOCuS", 1], mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0", est = cp, real = p$changepoint, N = p$N)

  # FoCUS 10
  print("FOCus0 p10")
  res <- mclapply(data, function (y) FOCuS_melk(y,  tlist["FOCuS 10", 1], mu0 = 0, grid = grid[round(seq(1, 50, length.out = 10))], K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0-10p", est = cp, real = p$changepoint, N = p$N))
  #print("page-CUSUM done")

  # Page CUSUM 50
  print("Page 50p")
  res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist["Page-CUSUM 50", 1], mu0 = 0, grid = grid), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <-  rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-50p", est = cp, real = p$changepoint, N = p$N))


  return(output)
}



output_file = "./simulations/pre-change-known/results/dr_comp4.RData"

grid <- find_grid(0, 50, .01, 1.3)

sim_grid <- expand.grid(
  N = 2e6,
  changepoint = 1e5,
  delta = c(seq(.1, .55, by = 0.02), grid[grid > .1 & grid < .6])
)

load("simulations/pre-change-known/thresholds.RData")

tlist <- thresholds %>%
  filter(run_len == 1e6) %>%
  group_by(algo) %>%
  summarise(tres = mean(threshold)) %>%
  column_to_rownames(var = "algo")

#run_simulation(sim_grid[10, ], NREP, tlist = tlist)

if (T) {
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
    det_delay = ifelse(est - real >= 0, est - real, NA),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & is.na(det_delay), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positive
  )





## detection delay ####

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 1, 3, 4, 5, 6)]
detection_delay <- ggplot(summary_df %>% filter(true_positive == 1),
                           aes(x = magnitude, y = det_delay, group = algo, col = algo)) +
  geom_vline(xintercept = grid, col = "grey") +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  xlim(.3, .7) +
  scale_y_log10() +
  ylab("Detection Delay") +
  theme_idris()

detection_delay
