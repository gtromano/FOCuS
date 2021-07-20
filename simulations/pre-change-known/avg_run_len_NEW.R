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
  grid <- find_grid(0, 26, .01, 1.74)
  FOCuSRUN <- mclapply(data, FOCuS_offline, thres = Inf, mu0 = 0, grid = NA, K = Inf, mc.cores = CORES)
  FOCuS10RUN <- mclapply(data, FOCuS_offline, thres = Inf, mu0 = 0, grid = grid[floor(seq(1, 25, length.out = 10))], K = Inf, mc.cores = CORES)

  page25RUN <- mclapply(data, PageCUSUM_offline, thres = Inf, mu0 = 0, grid = grid, mc.cores = CORES)
}

totalRUN <- list(FOCuSRUN, FOCuS10RUN, page25RUN)
#save.image(file = "simulations/pre-change-known/results/avg_run_len_NEW.RData")

thre_seq <- seq(14, 20, by = .05)
avg_run_len <- matrix(nr = length(thre_seq), nc = length(totalRUN))

row.names(avg_run_len) <- thre_seq

for (i in seq_along(thre_seq)) {
  for (j in seq_along(totalRUN)) {
    cat(thre_seq[i], j, "\n")
    avg_run_len[i, j] <- mean(sapply(totalRUN[[j]], run_len_calculator, thres = thre_seq[i]))
  }
}



colnames(avg_run_len) <- c("FOCuS", 'FOCuS 10', 'Page-CUSUM 25')

avg_run_len

save(avg_run_len, file = "simulations/pre-change-known/results/avg_run_len_NEW.RData")
tlist <- apply(avg_run_len, 2, function (len) thre_seq[which(len == 1e6)][1])

