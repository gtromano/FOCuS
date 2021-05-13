library(microbenchmark)
library(FOCuS)
library(tidyverse)

N <- c(1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5)
  
res <- list()

if(F) {
  for (n in N) {
    print(n)
    y <- rnorm(n)
    res[[length(res) + 1]] <- microbenchmark(FOCuS_offline(y, Inf), YuCUSUM_offline(y, Inf), times=30, unit = "ms")
    save.image("TIMING.RData")
  }
  names(res) <- N
  save.image("TIMING.RData")
  
}


load("TIMING.RData")



res[[as.character(n)]]


timings <- data.frame(res[[1]], len = 100)


for (n in N[2:length(N)]) {
  
  timings <- rbind(timings, 
                   data.frame(res[[as.character(n)]], len = n))

}


timings <- timings %>% mutate(milliseconds = time * 1e-6, algo = if_else(expr == "FOCuS_offline(y, Inf)", "FOCuS", "Yu-CUSUM"))


cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 6)]
timing_plot <- ggplot(timings, aes(x = len, y = milliseconds, group = algo, lty = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  #stat_summary(fun.data = "mean_se", geom = "errorbar") +
  scale_color_manual(values = cbPalette) +
  xlab("Run Length") +
  scale_x_log10() +
  theme_idris()





ggsave("simulations/pre-change-unknown/results/timings.pdf", ggarrange(timing_plot, legend = F), width = 6, height = 6)
