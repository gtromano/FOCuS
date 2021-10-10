library(parallel)
library(tidyverse)
# this is the example in the introduction, find here the methods
# as described in the kirch paper


.H <- function (x, mu0 = 0) x - mu0

MOSUM_offline_kirch <- function (Y, threshold, w,  FUN = .H) {
  const <- 1 / sqrt(w)
  stat <- sapply((w+1):length(Y), function (n){
    abs(sum(FUN(Y[(n-w):n]))) * const
  })

  stat <- c(rep(0, w), stat) # adding w zeroes at the start of the statistic

  cp <- which(stat >= threshold)[1]
  cp <- ifelse(is.na(cp), -1, cp)
  return(list(cp = cp, maxs = stat))
}

CUSUM_offline_kirch <- function (Y, threshold,  FUN = .H) {
  #stat <- sapply(1:length(Y), function (n){
  #  abs(sum(FUN(Y[1:n]))) * (1 / sqrt(n))
  #})
  stat <- abs(cumsum(Y)) * (1/sqrt(1:length(Y)))
  cp <- which(stat >= threshold)[1]
  cp <- ifelse(is.na(cp), -1, cp)
  return(list(cp = cp, maxs = stat))
}

CUSUM_offline_kirch(rnorm(100), Inf)

PageCUSUM_offline_kirch <- function (Y, threshold,  FUN = .H) {
  stat <- sapply(1:length(Y), function (n){
    #Q <- sapply(1:n-1, function (w) abs(sum(FUN(Y[(n-w):n]))) * (1 / sqrt(w + 1)) )
    Q <- abs(cumsum(Y[n:1])) * (1/sqrt(1:n))
    return(max(Q))
  })

  cp <- which(stat >= threshold)[1]
  cp <- ifelse(is.na(cp), -1, cp)
  return(list(cp = cp, maxs = stat))
}



#### penalty tuning, for a avg. run length up to 2000 observations

run_len_calculator <- function (res, thres) {
  n <- length(res$maxs)
  cp <- which(res$maxs >= thres)[1]
  ifelse(is.na(cp), n, cp)
}

REP <- 100
set.seed(42)
data <- lapply(1:REP, function (i) rnorm(2000))


cusumTR <- mclapply(data, CUSUM_offline_kirch, threshold = Inf, mc.cores = 6)
mosumTR <- mclapply(data, MOSUM_offline_kirch, threshold = Inf, w = 50, mc.cores = 6)
pageTR <- mclapply(data, PageCUSUM_offline_kirch, threshold = Inf, mc.cores = 6)


tre_seq <- seq(0.5, 5, length.out = 100)
avg_run_len <- matrix(nr = length(tre_seq), nc = 3)

row.names(avg_run_len) <- tre_seq

for (i in seq_along(tre_seq)) {
  avg_run_len[i, 1] <- mean(sapply(cusumTR, run_len_calculator, thres = tre_seq[i]))
  avg_run_len[i, 2] <- mean(sapply(mosumTR, run_len_calculator, thres = tre_seq[i]))
  avg_run_len[i, 3] <- mean(sapply(pageTR, run_len_calculator, thres = tre_seq[i]))
}



cusumT <- tre_seq[which(avg_run_len[, 1] == 2000)][1]
mosumT <- tre_seq[which(avg_run_len[, 2] == 2000)][1]
pageT <- tre_seq[which(avg_run_len[, 3] == 2000)][1]




### this is the common noise
set.seed(45)
base <- rnorm(2e3)

#### case 1, length 2000, mag of 1 ###

y1 <- c(rep(0, 1e3), rep(1, 1e3)) + base


res1c <- CUSUM_offline_kirch(y1, cusumT)$cp
res1m <- MOSUM_offline_kirch(y1, mosumT, 50)$cp
res1p <- PageCUSUM_offline_kirch(y1, pageT)$cp



#### case 2, length 2000, mag of 2 ###

y2 <- c(rep(0, 1e3), rep(.2, 1e3)) + base

res2c <- CUSUM_offline_kirch(y2, cusumT)$cp
res2m <- MOSUM_offline_kirch(y2, mosumT, 50)$cp
res2p <- PageCUSUM_offline_kirch(y2, pageT)$cp



####### NOW WE MAKE THE SEQUENCE LARGER

REP <- 100
set.seed(42)
data <- lapply(1:REP, function (i) rnorm(1e4))


cusumTR <- mclapply(data, CUSUM_offline_kirch, threshold = Inf, mc.cores = 6)
mosumTR <- mclapply(data, MOSUM_offline_kirch, threshold = Inf, w = 50, mc.cores = 6)
pageTR <- mclapply(data, PageCUSUM_offline_kirch, threshold = Inf, mc.cores = 6)


tre_seq <- seq(2, 6, length.out = 100)
avg_run_len <- matrix(nr = length(tre_seq), nc = 3)

row.names(avg_run_len) <- tre_seq

for (i in seq_along(tre_seq)) {
  avg_run_len[i, 1] <- mean(sapply(cusumTR, run_len_calculator, thres = tre_seq[i]))
  avg_run_len[i, 2] <- mean(sapply(mosumTR, run_len_calculator, thres = tre_seq[i]))
  avg_run_len[i, 3] <- mean(sapply(pageTR, run_len_calculator, thres = tre_seq[i]))
}



cusumT2 <- tre_seq[which(avg_run_len[, 1] == 10000)][1]
mosumT2 <- tre_seq[which(avg_run_len[, 2] == 10000)][1]
pageT2 <- tre_seq[which(avg_run_len[, 3] == 10000)][1]


set.seed(42)
y3 <- c(rnorm(8e3), y1)

res3c <- CUSUM_offline_kirch(y3, cusumT2)$cp
res3m <- MOSUM_offline_kirch(y3, mosumT2, 50)$cp
res3p <- PageCUSUM_offline_kirch(y3, pageT2)$cp




##### PLOTTING #####
library(ggplot2)


theme_idris <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey20"),
    panel.border =  element_rect(fill = NA,
                                 colour = "grey20")
  )
}


anomalies_plot <- function (y, t, real, estimatedFOCuS, estimatedHTM, truechange = 1000) {

  sdy <-  sd(diff(y))
  ymin <- min(y) - 0.3 * sdy

  len <- 3 * sdy

  # real cps
  realdf <- data.frame(x1 = real, y1 = ymin, y2 = ymin - len)
  realdf <- realdf %>% mutate(yc = (y1 + y2)/2, xst = truechange)
  realcps <- geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = realdf, col = 3)
  s1 <- geom_segment(aes(x = xst, xend = x1, y = yc, yend = yc), data = realdf, col = "grey", lty = 2)

  # est cps
  estdf <- data.frame(x1 = estimatedFOCuS, y1 = (ymin - len), y2 = (ymin - 2 * len))
  estdf <- estdf %>% mutate(yc = (y1 + y2)/2, xst = truechange)
  FOCuS_est = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = estdf, col = 2)
  s2 <- geom_segment(aes(x = xst, xend = x1, y = yc, yend = yc), data = estdf, col = "grey", lty = 2)


  # est numenta
  estHTM <- data.frame(x1 = estimatedHTM, y1 =  (ymin - 2 * len), y2 = (ymin - 3 * len))
  estHTM <- estHTM %>% mutate(yc = (y1 + y2)/2, xst = truechange)
  HTM_est <- geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = estHTM, col = 4)
  s3 <- geom_segment(aes(x = xst, xend = x1, y = yc, yend = yc), data = estHTM, col = "grey", lty = 2)

  exe <- ggplot(data.frame(t = t, y), aes(x = t, y = y)) +
    geom_line() +
    ylab("y") +
    theme_idris() +
    theme(legend.position = "none")

  complete_plot <- exe + realcps + FOCuS_est + HTM_est + s1 + s2 + s3
  complete_plot
}


# first scenario
p1 <- anomalies_plot(y1, 1:2000, res1c, res1m, res1p) +
  geom_vline(xintercept = 1000, col = "grey")

ggsave(p1, filename = "simulations/intro_plot_1.png", width = 8, height = 2)

# second scenario
p2 <- anomalies_plot(y2, 1:2000, res2c, res2m, res2p) +
  geom_vline(xintercept = 1000, col = "grey")
ggsave(p2, filename = "simulations/intro_plot_2.png", width = 8, height = 2)

# thirds scenario
p3 <- anomalies_plot(y3, 1:1e4, res3c, res3m, res3p, truechange = 9e3) +
  geom_vline(xintercept = 9e3, col = "grey") +
  xlim(8000, 10000)
ggsave(p3, filename = "simulations/intro_plot_3.png", width = 8, height = 2)
