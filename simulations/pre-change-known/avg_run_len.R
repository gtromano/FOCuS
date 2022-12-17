# this scripts generates the average run length figure and thresholds

source("simulations/helper_functions.R")

# calculates the run length, if it goes over the length of the sequence
run_len_calculator <- function (res, thres) {
  n <- length(res$maxs)
  cp <- which(res$maxs >= thres)[1]
  ifelse(is.na(cp), n, cp)
}

SEED <- 45
CORES <- 16
REP <- 100
N <- 2e6
set.seed(SEED)
data <- lapply(1:REP, function (i) rnorm(N))


if (T) {
  grid <- find_grid(0, 21, .01, 1.74)
  FOCuSRUN <- mclapply(data, FOCuS, thres = Inf, mu0 = 0, grid = NA, K = Inf, mc.cores = CORES)
  FOCuS10RUN <- mclapply(data, FOCuS, thres = Inf, mu0 = 0, grid = grid[c(1, 3, 6, 8, 11, 10, 13, 15, 18, 20)], K = Inf, mc.cores = CORES)
  page25RUN <- mclapply(data, PageCUSUM_offline, thres = Inf, mu0 = 0, grid = grid, mc.cores = CORES)
  page10RUN <- mclapply(data, PageCUSUM_offline, thres = Inf, mu0 = 0, grid = grid[c(1, 3, 6, 8, 11, 10, 13, 15, 18, 20)], mc.cores = CORES)
  wins <- unique(abs(18 / grid)) %>% round()
  MOSUMRUN <- mclapply(data, MOSUM_offline_kirch, thres = Inf, W = wins, mc.cores = CORES)

  # merging the total runs
  totalRUN <- list(FOCuSRUN, FOCuS10RUN, page25RUN, page10RUN, MOSUMRUN)
  rm(FOCuSRUN, FOCuS10RUN, page25RUN, page10RUN, MOSUMRUN)

}

# summarising
thre_seq <- seq(1, 20, by = .1)
avg_run_len <- matrix(nr = length(thre_seq), nc = length(totalRUN))

row.names(avg_run_len) <- thre_seq
colnames(avg_run_len) <- c("FOCuS", 'FOCuS 10', 'Page-CUSUM 20', 'Page-CUSUM 10', 'MOSUM')

# a bit inefficient
if (T) {
  for (i in seq_along(thre_seq)) {
    for (j in seq_along(totalRUN)) {
      cat(thre_seq[i], j, "\n")

      avg_run_len[i, j] <- mean(mclapply(totalRUN[[j]], run_len_calculator, thres = thre_seq[i], mc.cores = CORES) %>% unlist)
    }
  }
  avg_run_len

  save(avg_run_len, file = "simulations/pre-change-known/results/avg_run_len_NEW6.RData")

  # getting the threshold list
  tlist <- apply(avg_run_len, 2, function (len) thre_seq[which(len >= 1e6)][1])
  tlist <- lapply(tlist, function (x) x)
  names(tlist) <- colnames(avg_run_len)
  save(tlist, file = "simulations/pre-change-known/results/tlist.RData")
  tlist

}

### ggplot ###
load("simulations/pre-change-known/results/avg_run_len_NEW7.RData")
plotDF <- as.data.frame(avg_run_len) %>%
   add_column(threshold = thre_seq) %>%
   pivot_longer(names_to = "algo", values_to = "avg_run_len", - threshold)

plotDF[plotDF$algo == "MOSUM", ]$avg_run_len <- sqrt(plotDF[plotDF$algo == "MOSUM", ]$avg_run_len)

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(3, 4, 2, 6, 5)]
ggplot(plotDF %>% filter(algo != "FOCuSmelk", avg_run_len < 1.5e6)) +
   geom_line(aes(x = threshold, y = avg_run_len, group = algo, col = algo)) +
   scale_y_log10() +
   xlim(1, 15) +
   scale_color_manual(values = cbPalette) +
   ylab("Run Length") +
   geom_hline(yintercept = 1e6, col = "grey", lty = 2) +
   theme_idris() +
   theme(legend.position = "none") + hugefonts()