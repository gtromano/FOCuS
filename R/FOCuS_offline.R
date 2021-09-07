FOCuS_offline <- function(y, thres, mu0 = NA, training_data = NA, grid = NA, K = Inf) {
  
  # checks on the data generating function
  if( !is.numeric(y))
    stop("y must be a vector of observations")

  # checks on the threshold
  if( !is.numeric(thres) | thres < 0)
    stop("thres must be a positive numeric")
  
  # checks on the mu0
  if(!is.na(mu0))
    if(!is.numeric(mu0) | length(mu0) > 1)
      stop("mu0 must be a numeric value")
  
  # checks on the grid
  if(!is.na(grid[1]))
    if(!is.numeric(grid))
      stop("mu0 must be a numeric value")
  
  # checks on the K
  if(!is.na(K))
    if(!is.numeric(K) | K <= 0)
      stop("K must be a positive numeric")
  
  
  # running the function
  out <- .FoCUS_offline(y, thres, mu0, training_data, grid, K)
  out$changepoint <- out$t + out$Q1[[which.max(sapply(out$Q1, function(q) q$max))]]$a * 2
  return(out)
}