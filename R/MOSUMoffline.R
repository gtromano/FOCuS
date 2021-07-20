.MOSUMstep <- function (Info, new_point, FUN) {
  Info$window <- c(Info$window[2:length(Info$window)],  new_point) # we update the window
  Info$S <- abs(sum(FUN(Info$window)))
  return(Info)
}

.H <- function (x) x

MOSUM_offline <- function (Y, threshold, w, mu0 = 0, FUN = .H) {

  Info <- list(window = Y[1:w] , S = 0)
  cp <- -1
  t <- w

  max_storage <- NULL

  for (y in Y[(w+1):length(Y)]) {
    t <- t + 1
    Info <- .MOSUMstep(Info, y - mu0, FUN)
    max_storage <- c(max_storage, Info$S)
    if (Info$S >= threshold) {
      cp <- t
      break
    }
  }
  return(list(cp = cp, maxs = max_storage))
}



MOSUM_offline_kirch <- function (Y, threshold, w, mu0 = 0, FUN = .H) {
  const <- 1 / sqrt(w)
  stat <- sapply((w+1):length(Y), function (t){
    abs(sum(FUN(Y[(t-w):t]))) * const
  })

  cp <- which(stat >= threshold)[1] + w
  cp <- ifelse(is.na(cp), -1, cp)
  return(list(cp = cp, maxs = stat[1:(w-cp)]))
}
