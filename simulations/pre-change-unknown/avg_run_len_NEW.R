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
CORES <- 6
REP <- 100
N <- 2e6

# training
train <- lapply(1:REP, function (i) rnorm(1e5))

set.seed(SEED)
data <- lapply(1:REP, function (i) rnorm(N))

if (T) {
  grid <- find_grid(0, 26, .01, 1.74)
  FOCuSRUN <- mclapply(data, FOCuS_offline, thres = Inf, grid = NA, K = Inf, mc.cores = CORES)
  FOCuSRUNtrain <- mclapply(data, FOCuS_offline, thres = Inf, training_data = train[[k]][1:1e5], grid = NA, K = Inf, mc.cores = CORES)
  FOCuSRUN1e3 <- mclapply(1:REP, function (k) FOCuS_offline(data[[k]], thres = Inf, mu0 = mean(train[[k]][1:1e3]), grid = NA, K = Inf), mc.cores = CORES)
  FOCuSRUN1e4 <- mclapply(1:REP, function (k) FOCuS_offline(data[[k]], thres = Inf, mu0 = mean(train[[k]][1:1e4]), grid = NA, K = Inf), mc.cores = CORES)
  FOCuSRUN1e5 <- mclapply(1:REP, function (k) FOCuS_offline(data[[k]], thres = Inf, mu0 = mean(train[[k]][1:1e5]), grid = NA, K = Inf), mc.cores = CORES)
  FOCuSRUNInf <- mclapply(1:REP, function (k) FOCuS_offline(data[[k]], thres = Inf, mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)

}

totalRUN <- list(FOCuSRUN, FOCuSRUNtrain, FOCuSRUN1e3, FOCuSRUN1e4, FOCuSRUN1e5, FOCuSRUNInf)

thre_seq <- seq(1, 20, by = .05)
avg_run_len <- matrix(nr = length(thre_seq), nc = length(totalRUN))

row.names(avg_run_len) <- thre_seq

for (i in seq_along(thre_seq)) {
  for (j in seq_along(totalRUN)) {
    cat(thre_seq[i], j, "\n")
    avg_run_len[i, j] <- mean(sapply(totalRUN[[j]], run_len_calculator, thres = thre_seq[i]))
  }
}




colnames(avg_run_len) <- c("FOCuS", "FOCuS-t", paste("FOCuS0", 1e3), paste("FOCuS0", 1e4), paste("FOCuS0", 1e5), paste("FOCuS0", Inf))
avg_run_len

save(avg_run_len, file = "simulations/pre-change-unknown/results/avg_run_len_NEW2.RData")

#load("simulations/pre-change-unknown/results/avg_run_len_NEW.RData")
tlist <- apply(avg_run_len, 2, function (len) thre_seq[which(len >= 1e6-1)][1])
tlist
save(tlist, "simulations/pre-change-unknown/tlist.RData")


plotDF <- as.data.frame(avg_run_len) %>%
  add_column(threshold = thre_seq) %>%
  pivot_longer(names_to = "algo", values_to = "avg_run_len", - threshold)

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 3, 4, 5, 6)]
ggplot(plotDF %>% filter(algo != "FOCuSmelk")) +
  geom_line(aes(x = threshold, y = avg_run_len, group = algo, col = algo)) +
  scale_y_log10() +
  xlim(1, 15) +
  scale_color_manual(values = cbPalette) +
  ylab("Run Length") +
  geom_hline(yintercept = 1e6, col = "grey", lty = 2) +
  theme_idris() +
  theme(legend.position = "none")
