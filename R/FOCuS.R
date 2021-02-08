FOCuS <- function(dataFUN, thres, mu0 = NA, grid = NA, K = Inf) {
  
  # checks on the data generating function
  if(is.function(match.fun(dataFUN)))
    out_check <- dataFUN()
  else
    stop("Please provvide a data generating function.")
  if(!is.numeric(out_check)) stop("Data generating function provvided does not return a numeric output.")
  if(length(out_check) != 1) {
    warning("Length of the output from data generating function is greater than 1. Using only the first element.")
    f <- function() return(dataFUN()[1])
  } else {
    f <- match.fun(dataFUN)
  }
  
  # checks on the threshold
  if( !is.numeric(thres) | thres < 0)
    stop("thres must be a positive numeric")
  
  # checks on the mu0
  if(!is.na(mu0))
    if(!is.numeric(mu0) | length(mu0) > 1)
      stop("mu0 must be a numeric value")
  
  # checks on the grid
  if(!is.na(grid))
    if(!is.numeric(grid))
      stop("mu0 must be a numeric value")
  
  # checks on the K
  if(!is.na(K))
    if(!is.numeric(K) | K <= 0)
      stop("K must be a positive numeric")
  
  
  # running the function
  out <- .FoCUS(dataFUN, thres, mu0, grid, K)
  out$changepoint <- out$t + out$Q1[[which.max(sapply(out$Q1, function(q) q$max))]]$a * 2
  return(out)
}