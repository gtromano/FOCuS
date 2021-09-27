find_grid <- function(h0_m, max_ints, smallest_mag, g) {
  
  interval_lengths <- vector(length = max_ints/2)
  interval_lengths[1] <- smallest_mag
  for (i in 2:length(interval_lengths))
    interval_lengths[i] <- interval_lengths[i - 1] * g
  
  return(c(h0_m - interval_lengths[length(interval_lengths):1], h0_m + interval_lengths))
}
