library(FOCuS)
library(parallel)
library(ggpubr)

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

run_simulation <- function(p, REPS) {
  print(p)

  data <- lapply(1:REPS, function (k) c(rnorm(p$changepoint,0), rnorm(p$N - p$changepoint, p$delta)))

  # FOCuS with no pruning costraint
  res <- lapply(data, function (y) FOCuS_offline(y, p$threshold))
  cp <- sapply(res, function (r) r$cp)
  res_FOCuS <- data.frame(magnitude = p$delta, algo = "FOCuS", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold, grid_points = NA)
  #print("FOCus done")

  # Page CUSUM
  grid <- find_grid(0, p$grid_points, -1/1e10)
  cp <- sapply(data, function (y) pageCUSUM_offline(y, p$threshold, grid = grid))
  res_page <- data.frame(magnitude = p$delta, algo = "Page-CUSUM", est = cp, real = p$changepoint, N = p$N, threshold = p$threshold, grid_points = p$grid_points)
  #print("page-CUSUM done")

  return(rbind(res_FOCuS, res_page))
}
