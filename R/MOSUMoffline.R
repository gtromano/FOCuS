.MOSUMstep <- function (Info, new_point, m) {

  Info$window <- c(Info$window[2:length(Info$window)],  new_point)
  Info$S <- abs(sum(.H(Info$window, m)))
  return(Info)
}

.H <- function (x, m) (x - m)


MOSUM_offline <- function (Y, threshold, w, m = 0) {

  Info <- list(window = Y[1:w] , S = grid)
  cp <- -1
  t <- w

  for (y in Y[(w+1):length(Y)]) {
    t <- t + 1
    Info <- .MOSUMstep(Info, y, m)
    if (Info$S >= threshold) {
      cp <- t
      break
    }
  }
  return(cp)
}


MOSUM_offline(c(rnorm(100), rnorm(100, -2)), 15, 20)
