.CUSUMstep <- function (S, new_point, m) {

  S <- S + .H(new_point, m)
  return(S)
}

.H <- function (x, m) (x - m)

CUSUM_offline <- function (Y, threshold, m = 0) {

  S <- 0
  cp <- -1
  t <- 0

  for (y in Y) {
    t <- t + 1
    S <- .CUSUMstep(S, y, m)
    if (abs(S) >= threshold) {
      cp <- t
      break
    }
  }
  return(cp)
}


CUSUM_offline(c(rnorm(1e5), rnorm(100, -2)), 200, 0)
