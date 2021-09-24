library(microbenchmark)
library(FOCuS)
library(tidyverse)
source("simulations/helper_functions.R")

N <- c(1e2, 5e2, 1e3, 2e3, 3e3, 5e3, (1:5) * 1e4)
  
res <- list()

if(F) {
  for (n in N) {
    print(n)
    y <- rnorm(n)
    res[[length(res) + 1]] <- microbenchmark(FOCuS_offline(y, Inf, mu0 = 0),
                                             FOCuS_offline(y, Inf),
                                             FOCuS_offline(y, Inf, K = 5),
                                             YuCUSUM_offline(y, Inf),
                                             times = 30, unit = "ms")
    save.image("TIMING.RData")
  }
  names(res) <- N
  save(res, N, file = "TIMING.RData")
  
}



res[[as.character(n)]]


timings <- tibble(res[[1]], len = 100)


for (n in N[2:length(N)]) {
  
  timings <- rbind(timings, 
                   data.frame(res[[as.character(n)]], len = n))

}


levels(timings$expr) <- c("FOCuS0", "FOCuS", "R-FOCuS", "Yu-CUSUM")

timings <- timings %>% mutate(milliseconds = time * 1e-6) %>% rename(algo = expr)


theme_idris_remastered <- function(legend = c("nuke_it", "leave_it_this_time")) {
  th <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey20"),
    panel.border =  element_rect(fill = NA,
                                 colour = "grey20")
  )
  legend <- match.arg(legend)
  switch(legend,
         nuke_it = th + theme(legend.position = "none"),
         leave_it_this_time = th)
}




cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(4,2, 3, 6)]
timing_plot <- ggplot(timings, aes(x = len, y = milliseconds, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  scale_color_manual(values = cbPalette) +
  stat_function(fun = function (len) len/1000, colour = "grey", lty = 2) +
  stat_function(fun = function (len) (len/400)^2, colour = "grey", lty = 3) +
  xlab("Run Length") +
  scale_x_log10() +
  scale_y_log10() +
  theme_idris_remastered()
timing_plot


library(ggrepel)

data_label <- timings                              # Modify data
data_label$label <- NA
data_label$label[which(data_label$len == max(data_label$len))] <- as.character(data_label$algo[which(data_label$len == max(data_label$len))])

to_select <- NULL
for(algo in levels(data_label$algo))
  to_select <- c(to_select, which(algo == data_label$label)[1])

data_label <- data_label[to_select, ]

timing_plot <- ggplot(timings, aes(x = len, y = milliseconds, group = algo, col = algo)) +
  stat_summary(fun.data = "mean_se", geom = "line") +
  scale_color_manual(values = cbPalette) +
  stat_function(fun = function (len) len/1000, colour = "grey", lty = 2) +
  stat_function(fun = function (len) (len/400)^2, colour = "grey", lty = 3) +
  xlab("Run Length") +
    geom_label_repel(aes(label = label), nudge_x = 1000000, force = 10, show.legend = F, data = data_label) +
  scale_x_log10() +
  scale_y_log10() +
  theme_idris_remastered()
timing_plot

ggsave("simulations/pre-change-unknown/results/timings.pdf", ggarrange(timing_plot, legend = F), width = 6, height = 6)
