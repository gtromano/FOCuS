# a simple vectorization of
.v_max <- function (v, ...) sapply(v, function (a) max(a, ...))

.get_Qn <- function(Qold, mu, x) {
  Qold <- Qold + mu * (x - (mu / 2))
  return(.v_max(Qold, 0))
}

.pageCUSUM_step <- function(Info, new_point) {
  .get_Qn(Info$Q, Info$grid, new_point)
}

pageCUSUM_offlineR <- function (Y, threshold, grid) {

  Info <- list(Q = 0, grid = grid)
  cp <- -1
  t <- 0
  max_storage = NULL

  for (y in Y) {
    t <- t + 1
    Info$Q <- .pageCUSUM_step(Info, y)
    maximum <- max(Info$Q)
    max_storage <- c(max_storage, maximum)
    if (maximum >= threshold) {
      cp <- t
      break
    }
  }
  return(list(cp = cp, maxs = max_storage))
}
