# D <- function (s, t, x) {
#   abs(sqrt((t-s)/(t * s)) * sum( x[1:s] ) - sqrt((s)/(t * (t- s))) * sum( x[(s+1):t]))
# }
#
# yuCUSUM <- function (x, thres) {
#
#   cp <- -1
#   max_storage <- NULL
#
#   for (t in 2:length(x)) {
#     maximum <- max(sapply(1:(t-1), D, t = t, x = x))
#
#     max_storage <- c(max_storage, maximum)
#
#     if (maximum >= thres) {
#       cp <- t
#       break
#     }
#
#   }
#   return(list(cp = cp, maxs = max_storage))
# }



# system.time(res <- yuCUSUM(c(rnorm(1000, 10)), Inf))


# if (FALSE) {
#   N <- 1e3
#   REPS <- 100
#   data <- lapply(1:REPS, function (k) rnorm(N, 5))
#   res <- mclapply(data, function (y) yuCUSUM(y, Inf), mc.cores = 6)
#   max1e3 <- sapply(res, function (r) max(r$maxs))
#   save(max1e3yu, file = "yu.RData")
# }
#



####### yu efficient

D2 <- function (s, t, sums) {
    sqrt((t)/(s * (t - s))) * abs((s/t) * sums[t] -  sums[s])
}

yuCUSUM_v2 <- function (x, thres) {

  cp <- -1
  max_storage <- NULL
  sums <- cumsum(x)

  for (t in 2:length(x)) {
    maximum <- max(sapply(1:(t-1), D2, t = t, sums = sums))

    max_storage <- c(max_storage, maximum)

    if (maximum >= thres) {
      cp <- t
      break
    }

  }
  return(list(cp = cp, maxs = max_storage))
}

# not truly online but avoids R for cycle
yuCUSUM_v3 <- function (x, thres) {

  cp <- -1
  sums <- cumsum(x)

  max_storage <- sapply (2:length(x), function(t) {
    maximum <- max(sapply(1:(t-1), D2, t = t, sums = sums))
    maximum
  })

  return(list(cp = which(max_storage > thres)[1], maxs = max_storage))
}


if (FALSE) {
  N <- 1e4
  REPS <- 100
  data <- lapply(1:REPS, function (k) rnorm(N, 5))
  res <- mclapply(data, function (y) yuCUSUM_v2(y, Inf), mc.cores = 6)
  max1e4yu <- sapply(res, function (r) max(r$maxs))
  save(max1e4yu, file = "yu1e4.RData")
}
