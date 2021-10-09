#' Compute a geometric grid of values.
#'
#' @param mu0 
#' @param max_ints 
#' @param smallest_mag 
#' @param g 
#'
#' @return
#' @export
#'
#' @examples
find_grid <- function(mu0 = 0, max_ints = 20, smallest_mag = .001, g = .5) {
  
  interval_lengths <- vector(length = max_ints/2)
  interval_lengths[1] <- smallest_mag
  for (i in 2:length(interval_lengths))
    interval_lengths[i] <- interval_lengths[i - 1] * g
  
  return(c(h0_m - interval_lengths[length(interval_lengths):1], h0_m + interval_lengths))
}
