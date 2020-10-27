find_grid <- function(h0_m, max_ints, smallest_mag) {
  
  interval_lengths <- vector(length = max_ints/2)
  interval_lengths[1] <- smallest_mag
  for (i in 2:length(interval_lengths))
    interval_lengths[i] <- interval_lengths[i - 1] * 2
  
  return(c(h0_m - interval_lengths[length(interval_lengths):1], h0_m + interval_lengths))
}


find_grid_ocd <- function(h0_m, p, beta) {
  l <- 0:floor(log(p, 2))
  grid <- beta / sqrt((2 ^ l) * log(2 * p, 2))
  b0 <- beta / sqrt((2 ^ (floor(log(p, 2)) + 1)) * log(2 * p, 2))
  
  return(c(-grid[length(grid):1], -b0, h0_m, b0, grid))
}

