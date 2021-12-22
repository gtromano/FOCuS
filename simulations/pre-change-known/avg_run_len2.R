
source("simulations/helper_functions.R")

# calculates the run length, if it goes over the length of the sequence
run_len_calculator <- function (res, thres) {
  n <- length(res)
  cp <- which(res >= thres)[1]
  ifelse(is.na(cp), n, cp)
}

SEED <- 45
CORES <- 16
REP <- 80   # replicates per experiment
N <- 2e6
change_factor <- 10
set.seed(SEED)
data <- lapply(1:REP, function (i) rnorm(N))



if (T) {
  W <- 100
  K <- 15

  simRUN <- mclapply(data, function (y) {
    grid <- find_grid(0, 21, .01, 1.74)

    FOCuSRUN <- FOCuS(y, thres = Inf, mu0 = 0, grid = NA, K = Inf)$maxs
    FOCuS10RUN <- FOCuS(y, thres = Inf, mu0 = 0, grid = grid[c(1, 3, 6, 8, 11, 10, 13, 15, 18, 20)], K = Inf)$maxs
    cat("FOCuS done\n")

    page20RUN <- PageCUSUM_offline(y, thres = Inf, mu0 = 0, grid = grid)$maxs
    page10RUN <- PageCUSUM_offline(y, thres = Inf, mu0 = 0, grid = grid[c(1, 3, 6, 8, 11, 10, 13, 15, 18, 20)])$maxs
    cat("Page done\n")

    wins <- unique(abs(18 / grid)) %>% round()
    MOSUMRUN <- MOSUM_offline_kirch(y, thres = Inf, W = wins)$maxs
    cat("MOSUM done\n")


    return(data.frame(FOCuS = FOCuSRUN, "FOCuS 10" = FOCuS10RUN, 'Page-CUSUM 20' = page20RUN, 'Page-CUSUM 10' = page10RUN, 'MOSUM' = MOSUMRUN))

  }, mc.cores = CORES)


  #save(simRUN, file = "simulations/pre-change-known/results/avg_run_len_NEW8.RData")

}

cat("All done\n")



#load("simulations/pre-change-known/results/avg_run_len_NEW8.RData")


# checking broken sims
# sapply(NPFOCuSRUN, is.character)


# summarising
thre_seq <- seq(0.1, 20, by = .1)


avg_run_len <- mclapply(thre_seq, function (thres) {
  focusdf <- lapply(simRUN, function (run) {
    runlen <- apply(run, 2, run_len_calculator, thres)
    data.frame(algo = names(run), threshold = thres, run_len = runlen)
  })

  Reduce(rbind, focusdf)

}, mc.cores = 16)
avg_run_len <- Reduce(rbind, avg_run_len)

save(avg_run_len, file = "simulations/pre-change-known/results/avg_run_len_summ.RData")

avg_run_len <- avg_run_len %>% mutate(fp = run_len < 500)

summary_1 <- avg_run_len %>% group_by(threshold, algo) %>% summarise(avgrun = mean(run_len), fpr = mean(fp))


cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(3, 4, 2, 6, 5)]
ggplot(summary_1, aes(x = threshold, y = avgrun, col = algo)) +
  geom_line() +
  scale_y_log10() +
   xlim(1, 15) +
   scale_color_manual(values = cbPalette) +
   ylab("Run Length") +
   geom_hline(yintercept = 1e6, col = "grey", lty = 2) +
   theme_idris() +
   theme(legend.position = "none")

ggplot(summary_1, aes(x = threshold, y = fpr, col = algo)) +
  geom_line() +
  scale_x_log10() +
  facet_grid(~scenario)

