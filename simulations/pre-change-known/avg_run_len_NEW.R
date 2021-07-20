library(parallel)
library(tidyverse)
source("simulations/helper_functions.R")

# calculates the run lenght, if it goes over the length of the sequence
run_len_calculator <- function (res, thres) {
  n <- length(res$maxs)
  cp <- which(res$maxs >= thres)[1]
  ifelse(is.na(cp), n, cp)
}


CORES <- 16
REP <- 100
N <- 1e6 # up to a million observations
set.seed(42)
data <- mclapply(1:REP, function (i) rnorm(N), mc.cores = CORES)


if (T) {
  grid <- find_grid(0, 50, .01, 1.3)
  FOCuSRUN <- mclapply(data, FOCuS_melk, thres = Inf, mu0 = 0, grid = NA, K = Inf, mc.cores = CORES)
  FOCuS10RUN <- mclapply(data, FOCuS_melk, thres = Inf, mu0 = 0, grid = grid[round(seq(1, 50, length.out = 25))], K = Inf, mc.cores = CORES)

  page50RUN <- mclapply(data, PageCUSUM_offline, thres = Inf, mu0 = 0, grid = grid, mc.cores = CORES)
  page25RUN <- mclapply(data, PageCUSUM_offline, thres = Inf, mu0 = 0, grid = grid[round(seq(1, 50, length.out = 25))], mc.cores = CORES)
}

save.image(file = "simulations/pre-change-known/results/avg_run_len_NEW.RData")

#cusumTR <- mclapply(data, CUSUM_offline_kirch, threshold = Inf, mc.cores = CORES)
#mosumTR <- mclapply(data, MOSUM_offline_kirch, threshold = Inf, w = 50, mc.cores = CORES)
#pageTR <- mclapply(data, PageCUSUM_offline_kirch, threshold = Inf, mc.cores = CORES)
#
# tre_seq <- seq(0.5, 5, length.out = 100)
# avg_run_len <- matrix(nr = length(tre_seq), nc = 3)
#
# row.names(avg_run_len) <- tre_seq
#
# for (i in seq_along(tre_seq)) {
#   avg_run_len[i, 1] <- mean(sapply(cusumTR, run_len_calculator, thres = tre_seq[i]))
#   avg_run_len[i, 2] <- mean(sapply(mosumTR, run_len_calculator, thres = tre_seq[i]))
#   avg_run_len[i, 3] <- mean(sapply(pageTR, run_len_calculator, thres = tre_seq[i]))
# }
#
#
