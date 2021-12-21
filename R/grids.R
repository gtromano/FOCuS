#' Compute a geometric grid of values.
#'
#' @param mu0 The center of the grid, corresponding to the pre-change location.
#' @param max_ints The maximum amounts of intervals to generate.
#' @param smallest_mag The smallest absolute difference from the center of the grid.
#' @param g Parameter to regulate the geometric increase of the grid.
#'
#' @return A vector of grid values.
#' @export
#'
#' @examples
#'
#' set.seed(42)
#' y <- c(rnorm(3e5, 0), rnorm(1e4, 1))
#' g <- find_grid()
#' FOCuS(y, 18, mu0 = 0, grid = g)

find_grid <- function(mu0 = 0, max_ints = 21, smallest_mag = .01, g = 1.74) {
  
  interval_lengths <- vector(length = max_ints/2)
  interval_lengths[1] <- smallest_mag
  for (i in 2:length(interval_lengths))
    interval_lengths[i] <- interval_lengths[i - 1] * g
  
  return(c(mu0 - interval_lengths[length(interval_lengths):1], mu0 + interval_lengths))
}
