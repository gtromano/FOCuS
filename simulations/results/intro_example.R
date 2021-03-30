source("../helper_functions.R")
library(ggpubr)


anomalies_plot <- function (y, t, cp, estimatedMOSUM, estimatedCUSUM, estimatedPAGE) {

  sdy <-  sd(diff(y))
  ymin <- min(y) - 0.3 * sdy

  len <- 2 * sdy

  # real mosum
  mosumdf <- data.frame(x1 = estimatedMOSUM, y1 = ymin, y2 = ymin - len)
  mosumcps = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = mosumdf, col = 3)

  # est mosum
  cusumdf <- data.frame(x1 = estimatedCUSUM, y1 = (ymin - len), y2 = (ymin - 2 * len))
  cusumest = geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = cusumdf, col = 2)

  # est page
  pagedf <- data.frame(x1 = estimatedPAGE, y1 =  (ymin - 2 * len), y2 = (ymin - 3 * len))
  pageest <- geom_segment(aes(x = x1, xend = x1, y = y1, yend = y2), data = pagedf, col = 4)

  exe <- ggplot(data.frame(t = t, y), aes(x = t, y = y)) +
    geom_line() +
    ylab("y") +
    geom_vline(xintercept = cp, col = "grey", lty = 2) +
    theme_idris() +
    theme(legend.position = "none")

  complete_plot <- exe + mosumcps + cusumest + pageest
  complete_plot
}


THESEED <- 47



grid <- find_grid(0, 50, .01, 1.3)

load("../thresholds.RData")
tlist1e3 <- thresholds %>%
  filter(run_len == 1e3) %>%
  group_by(algo) %>%
  summarise(tres = quantile(threshold, .9)) %>%
  column_to_rownames(var = "algo")


# tweaking the penalties
REPS <- 100
N <- 1e3

data <- lapply(1:REPS, function (k) rnorm(N))
res <- mclapply(data, function (y) pageCUSUM_Kirch(y, Inf), mc.cores = 6)
max1e3 <- sapply(res, function (r) max(r$maxs))
tlist1e3["Page-CUSUM", 1] <- quantile(max1e3, .9)

# example 1 - normal on 0 with change at 1e3 at 1
set.seed(THESEED)
y1 <- c(rnorm(1e3, 0), rnorm(1e3, 1))

cpMOS <- MOSUMwrapper(y1, thres = tlist1e3["MOSUM", 1], bandw = 50)$t

cpCUSUM <- CUSUM_offline(y1, tlist1e3["CUSUM", 1], 0)$t

cpPage <- pageCUSUM_Kirch(y1, tlist1e3["Page-CUSUM", 1])$cp

plt1 <- anomalies_plot(y1, 1:2e3, 1e3, cpMOS, cpCUSUM, cpPage)
plt1

# example 2 - small magnitude (0.1)

set.seed(THESEED)
y2 <- c(rnorm(1e3, 0), rnorm(1e3, .5))

cpMOS <- MOSUMwrapper(y2, thres = tlist1e3["MOSUM", 1], bandw = 50)$t

cpCUSUM <- CUSUM_offline(y2, tlist1e3["CUSUM", 1], 0)$t

cpPage <- pageCUSUM_Kirch(y2, tlist1e3["Page-CUSUM", 1])$cp


plt2 <- anomalies_plot(y2, 1:2e3, 1e3, cpMOS, cpCUSUM, cpPage)
plt2


# example 3 - long term change
tlist1e5 <- thresholds %>%
  filter(run_len == 1e5) %>%
  group_by(algo) %>%
  summarise(tres = quantile(threshold, .9)) %>%
  column_to_rownames(var = "algo")


# takes roughly 30 mins to run
data <- lapply(1:REPS, function (k) rnorm(1e5))
res <- mclapply(data, function (y) pageCUSUM_Kirch(y, Inf), mc.cores = 6)
max1e3 <- sapply(res, function (r) max(r$maxs))
tlist1e5["Page-CUSUM", 1] <- quantile(max1e3, .9)


set.seed(THESEED)
y3 <- c(rnorm(99000, 0), y1)

cpMOS <- MOSUMwrapper(y3, thres = tlist1e5["MOSUM", 1], bandw = 50)$t

cpCUSUM <- CUSUM_offline(y3, tlist1e5["CUSUM", 1], 0)$t

cpPage <- pageCUSUM_Kirch(y3, tlist1e5["Page-CUSUM", 1])$cp


plt3 <- anomalies_plot(y3, 1:length(y3), 1e5, cpMOS, cpCUSUM, cpPage)
plt3 + xlim(9.1e4, 1e5 + 1e3)



# total plot with legend


legend <- ggplot(data.frame(y = rnorm(12), x = 1:12, algo = rep(c("MOSUM", "CUSUM", "Page-CUSUM"))), aes(x = x, y = y, col = algo)) +
  geom_line()
legend <- get_legend(legend)

final_plot <- ggarrange(ggarrange(plt1, plt2, plt3 + xlim(9.1e4, 1e5 + 1e3), labels = "AUTO", ncol = 1),
          legend, widths = c(3,1))
ggsave(final_plot,filename = "simulations/results/intro_example.pdf", device = "pdf", height = 4, width = 7)
