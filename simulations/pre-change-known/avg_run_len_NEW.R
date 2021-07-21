library(parallel)
library(tidyverse)
source("simulations/helper_functions.R")

# calculates the run lenght, if it goes over the length of the sequence
run_len_calculator <- function (res, thres) {
  n <- length(res$maxs)
  cp <- which(res$maxs >= thres)[1]
  ifelse(is.na(cp), n, cp)
}

SEED <- 45
CORES <- 16
REP <- 100
N <- 2e6 # up to a million observations
set.seed(SEED)
#data <- lapply(1:REP, function (i) rnorm(N))
data <- lapply(1:REP, function (i) rnorm(N, sd = .5))

if (T) {
  grid <- find_grid(0, 26, .01, 1.74)
  FOCuSRUN <- mclapply(data, FOCuS_offline, thres = Inf, mu0 = 0, grid = NA, K = Inf, mc.cores = CORES)
  FOCuS10RUN <- mclapply(data, FOCuS_offline, thres = Inf, mu0 = 0, grid = grid[c(3, 6, 8, 11, 13, 14, 16, 19, 21, 24)], K = Inf, mc.cores = CORES)

  page25RUN <- mclapply(data, PageCUSUM_offline, thres = Inf, mu0 = 0, grid = grid, mc.cores = CORES)
}

totalRUN <- list(FOCuSRUN, FOCuS10RUN, page25RUN)
#save.image(file = "simulations/pre-change-known/results/avg_run_len_NEW.RData")

#thre_seq <- seq(1, 20, by = .05)
thre_seq <- seq(1, 5, by = .03)
avg_run_len <- matrix(nr = length(thre_seq), nc = length(totalRUN))

row.names(avg_run_len) <- thre_seq

for (i in seq_along(thre_seq)) {
  for (j in seq_along(totalRUN)) {
    cat(thre_seq[i], j, "\n")
    avg_run_len[i, j] <- min(sapply(totalRUN[[j]], run_len_calculator, thres = thre_seq[i]))
  }
}



colnames(avg_run_len) <- c("FOCuS", 'FOCuS 10', 'Page-CUSUM 25')

avg_run_len

save(avg_run_len, file = "simulations/pre-change-known/results/avg_run_len_NEW2.RData")

#load("simulations/pre-change-known/results/avg_run_len_NEW2.RData")
tlist <- apply(avg_run_len, 2, function (len) thre_seq[which(len == 1e6)][1])
tlist

plotDF <- as.data.frame(avg_run_len) %>%
  add_column(threshold = thre_seq) %>%
  pivot_longer(names_to = "algo", values_to = "avg_run_len", - threshold)

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(3, 4, 2, 5, 6)]
ggplot(plotDF %>% filter(algo != "FOCuSmelk")) +
  geom_line(aes(x = threshold, y = avg_run_len, group = algo, col = algo)) +
  scale_y_log10() +
  xlim(1, 19) +
  scale_color_manual(values = cbPalette) +
  ylab("Run Length") +
  theme_idris() +
  theme(legend.position = "none")





###### plot of the detailed grid ########

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
y <- c(rnorm(1e5,0 ,.1), rnorm(2e6, 0.48288614, .1))


resF <- FOCuS_offline(y[1:maxY], thres = tlist["FOCuS"], mu0 = 0)

gridP <-  find_grid(0, 26, .01, 1.74)
#gridP[41] <- .5

resP <- PageCUSUM_offline(y[1:maxY], thres = tlist["Page-CUSUM 25"], mu0 = 0, grid = gridP)

# ggplot(mapping = aes(x = t, y = y),
#        data = tibble(t = (1e5 - 200):maxY, y = y[(1e5 - 200):maxY])) +
#   geom_point(aes(x = t, y = y), tibble(t = (1e5 - 200):(maxY + 100), y = y[(1e5 - 200):(maxY + 100)]), col = "grey") +
#   geom_point() + geom_line()


plot1 <- ggplot(tibble(mu = -3:3)) +
  stat_function(aes(x = mu), fun = function(x) plot_piecewise_quad(x, quad = resF$Q1), col = 4) +
  theme_idris() +
  xlim(0, 1.5) +
  geom_vline(xintercept = gridP, col = "grey", alpha = .3) +
  xlab(expression(mu)) + ylab(expression(Q[t](mu))) +
  geom_text(aes(x = grid, y = y, label = Q), data = tibble(grid = gridP + .02, y = round(resP$Q, 2) + .1, Q = round(resP$Q, 2)), col = "grey", alpha = .8) +
  geom_text(aes(x = .45, y = round(tail(resF$maxs, 1), 2) + .1, label = round(tail(resF$maxs, 1), 2)), col = 4) +
  geom_vline(xintercept = .482, lty = 2)


# off grid

set.seed(32)
y2 <- c(rnorm(1e5,0 ,.1), rnorm(2e6, sum(grid[21:22])/2, .1))


resF2 <- FOCuS_offline(y2[1:maxY], thres = tlist["FOCuS"], mu0 = 0)
resP2 <- PageCUSUM_offline(y2[1:maxY], thres = tlist["Page-CUSUM 25"], mu0 = 0, grid = gridP)

plot2 <- ggplot(tibble(mu = -3:3)) +
  stat_function(aes(x = mu), fun = function(x) plot_piecewise_quad(x, quad = resF2$Q1), col = 4) +
  theme_idris() +
  xlim(0, 1.5) +
  geom_vline(xintercept = gridP, col = "grey", alpha = .3) +
  xlab(expression(mu)) + ylab(expression(Q[t](mu))) +
  geom_text(aes(x = grid, y = y, label = Q), data = tibble(grid = gridP + .02, y = round(resP2$Q, 2) + .1, Q = round(resP2$Q, 2)), col = "grey", alpha = .8) +
  geom_text(aes(x = .61, y = round(tail(resF2$maxs, 1), 2) + .15, label = round(tail(resF2$maxs, 1), 2)), col = 4) +
  geom_vline(xintercept = sum(grid[21:22])/2, lty = 2)

plot2


ggsave(ggarrange(plot1, plot2, nrow = 2), filename = "grid-comp.pdf", width = 18, height = 12)



#### the grid plot ########

par(mfrow = c(1, 2))

ggplot(data.frame()) + geom_point() + xlim(-5, 5) + ylim(0, 1) +
  geom_vline(xintercept = grid) + geom_vline(xintercept = grid) +
  xlab(expression(mu)) +
  theme_idris() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

ggplot(data.frame()) + geom_point() + xlim(-5, 5) + ylim(0, 1) +
  geom_vline(xintercept = grid[c(3, 6, 8, 11, 13, 14, 16, 19, 21, 24)]) +
  xlab(expression(mu)) +
  theme_idris() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())
