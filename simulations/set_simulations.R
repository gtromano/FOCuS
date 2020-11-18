library(tidyverse)
library(FOCuS)
library(parallel)
library(ggpubr)

library(compiler)
compiler::enableJIT(3)
pageCUSUM_offline <- compiler::cmpfun(pageCUSUM_offline)


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

run_simulation <- function(p, REPS, seed = 42) {
  print(p)

  set.seed(seed)
  data <- lapply(1:REPS, function (k) c(rnorm(p$changepoint,0), rnorm(p$N - p$changepoint, p$delta)))

  # FOCuS with no pruning costraint
  res <- lapply(data, function (y) FOCuS_offline_sim(y, p$threshold, grid = NA))
  cp <- sapply(res, function (r) r$cp)
  res_FOCuS <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("FOCus done")

  # FoCUS 5
  grid <- find_grid(0, 5, .1)
  res <- lapply(data, function (y) FOCuS_offline_sim(y, p$threshold, grid = grid))
  cp <- sapply(res, function (r) r$cp)
  res_FOCuS5 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS 5", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("page-CUSUM done")


  # Page CUSUM 5
  grid <- find_grid(0, 5, .1)
  cp <- unlist(mclapply(data, function (y) pageCUSUM_offline(y, p$threshold, grid = grid), mc.cores = 6))
  res_page25 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-CUSUM 5", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)
  #print("page-CUSUM done")

  # Page CUSUM 50
  grid <- find_grid(0, 50, .1)
  cp <- unlist(mclapply(data, function (y) pageCUSUM_offline(y, p$threshold, grid = grid), mc.cores = 6))
  res_page50 <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-CUSUM 50", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold)


  return(rbind(res_FOCuS, res_FOCuS5, res_page25, res_page50))
}
