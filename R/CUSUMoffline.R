# .CUSUMstep <- function (S, new_point, m) {
#
#   S <- S + .H(new_point, m)
#   return(S)
# }
#
# .H <- function (x, m) (x - m)
#
# CUSUM_offline <- function (Y, threshold, m = 0) {
#
#   S <- 0
#   cp <- -1
#   t <- 0
#   max_storage <- NULL
#
#   for (y in Y) {
#     t <- t + 1
#     S <- .CUSUMstep(S, y, m)
#
#     max_storage <- c(max_storage, abs(S))
#
#     if (abs(S) >= threshold) {
#       cp <- t
#       break
#     }
#   }
#   return(list(cp = cp, maxs = max_storage))
# }
#
#
# #CUSUM_offline(c(rnorm(1e5), rnorm(100, -2)), 200, 0)
