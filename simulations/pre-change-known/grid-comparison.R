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



output_file <- "./simulations/pre-change-known/results/dr_comp4.RData"

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

#tlist[,1] <- 16
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



summary_df %>% filter(false_alarm == 1) %>% select(algo) %>% table()

to_exclude <- summary_df %>% filter(false_alarm == 1) %>% select(sim, magnitude) %>% unique()


for(i in 1:nrow(to_exclude)) {
  summary_df[summary_df$sim == to_exclude[i, 1] & summary_df$magnitude == to_exclude[i, 2], ]$det_delay <- NA
}

## detection delay ####

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 1, 3, 4, 5, 6)]
detection_delay <- ggplot(summary_df,
                           aes(x = magnitude, y = det_delay, group = algo, col = algo)) +
  geom_vline(xintercept = grid, col = "grey") +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("magnitude") +
  xlim(.3, .55) +
  scale_y_log10() +
  ylab("Detection Delay") +
  theme_idris()

detection_delay

fa_rate <- ggplot(summary_df %>% filter(algo != "MOSUM"), aes(x = magnitude, y = false_alarm, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  stat_summary(fun.data = "mean_se", geom = "errorbar") +
  geom_vline(xintercept = grid, col = "grey") +
  scale_color_manual(values = cbPalette) +
  ylim(0, 1) +
  xlab("magnitude") +
  ylab("False Alarm Rate") +
  xlim(.3, .55) +
  theme_idris()
fa_rate








#################### experiment 2 ##############################


load("simulations/pre-change-known/thresholds.RData")
tlist <- thresholds %>%
  filter(run_len == 1e6) %>%
  group_by(algo) %>%
  summarise(tres = mean(threshold)) %>%
  column_to_rownames(var = "algo")

maxY <- 1e5 + 50

### plotting

plot_piecewise_quad <- function(x, quad) {
  for (q in quad)
    for (i in q$ints)
      if (i$l < x && i$u >= x)
        return(q$a * x^2 + q$b * x + q$c)
}
plot_piecewise_quad <- Vectorize(plot_piecewise_quad, vectorize.args = "x")


set.seed(32)
y <- c(rnorm(1e5,0 ,.1), rnorm(2e6, .50, .1))


resF <- FOCuS_offline(y[1:maxY], thres = tlist[1,], mu0 = 0)

gridP <- find_grid(0, 50, .01, 1.3)
gridP[41] <- .5

resP <- PageCUSUM_offline(y[1:maxY], thres = tlist[3, ], mu0 = 0, grid = gridP)

# ggplot(mapping = aes(x = t, y = y),
#        data = tibble(t = (1e5 - 200):maxY, y = y[(1e5 - 200):maxY])) +
#   geom_point(aes(x = t, y = y), tibble(t = (1e5 - 200):(maxY + 100), y = y[(1e5 - 200):(maxY + 100)]), col = "grey") +
#   geom_point() + geom_line()


plot1 <- ggplot(tibble(mu = -3:3)) +
  stat_function(aes(x = mu), fun = function(x) plot_piecewise_quad(x, quad = resF$Q1), col = 4) +
  theme_idris() +
  xlim(0, 2) +
  geom_vline(xintercept = gridP, col = "grey", alpha = .3) +
  xlab(expression(mu)) + ylab(expression(Q[t](mu))) +
  geom_text(aes(x = grid, y = y, label = Q), data = tibble(grid = gridP + .02, y = round(resP$Q, 2) + .1, Q = round(resP$Q, 2)), col = "grey", alpha = .8) +
  geom_text(aes(x = .47, y = round(tail(resF$maxs, 1), 2) + .1, label = round(tail(resF$maxs, 1), 2)), col = 4) +
  geom_vline(xintercept = .5, lty = 2)


# off grid

set.seed(32)
y2 <- c(rnorm(1e5,0 ,.1), rnorm(2e6, .58, .1))


resF2 <- FOCuS_offline(y2[1:maxY], thres = tlist[1,], mu0 = 0)

gridP <- find_grid(0, 50, .01, 1.3)
gridP[41] <- .5

resP2 <- PageCUSUM_offline(y2[1:maxY], thres = tlist[3, ], mu0 = 0, grid = gridP)

plot2 <- ggplot(tibble(mu = -3:3)) +
  stat_function(aes(x = mu), fun = function(x) plot_piecewise_quad(x, quad = resF2$Q1), col = 4) +
  theme_idris() +
  xlim(0, 2) +
  geom_vline(xintercept = gridP, col = "grey", alpha = .3) +
  xlab(expression(mu)) + ylab(expression(Q[t](mu))) +
  geom_text(aes(x = grid, y = y, label = Q), data = tibble(grid = gridP + .02, y = round(resP2$Q, 2) + .1, Q = round(resP2$Q, 2)), col = "grey", alpha = .8) +
  geom_text(aes(x = .56, y = round(tail(resF2$maxs, 1), 2) + .15, label = round(tail(resF2$maxs, 1), 2)), col = 4) +
  geom_vline(xintercept = .58, lty = 2)

plot2



ggsave(ggarrange(plot1, plot2, nrow = 2), filename = "grid-comp.pdf", width = 18, height = 12)
