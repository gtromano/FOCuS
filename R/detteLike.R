# this is a simple change in mean likelihood ratio test based
# from dette


detteStat <- function (j, k, m, sd, x) {
  (m^-3) * (((m + j) ^ 2 * (k - j) ^ 2)) * (mean(x[1:(m + j)]) - mean(x[(m+j+1):k])) ^ 2 / sd
}

detteStat <- function (j, k, m, mu0, sd, x) {

  pre_change <- (mean(x[m:j]) + mu0) / 2

  (sd ^ -1) * (pre_change - mean(x[(j+1):k])) ^ 2
}


# not truly online but avoids R for cycle
detteCUSUM <- function (x, thres, m = 100) {

  cp <- -1
  sd <- var(x[1:m])
  mu0 <- mean(x[1:m])

  max_storage <- sapply (m:length(x), function(t) {
    maximum <- max(sapply(1:t, detteStat, k = t, m = m, mu0 = mu0, sd = sd, x = x))
    maximum
  })

  return(list(cp = which(max_storage > thres)[1], maxs = max_storage))
}

res <- detteCUSUM(c(rnorm(1e3), rnorm(1000, 1)), Inf)
plot(res$maxs)
