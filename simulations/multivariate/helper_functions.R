library(ocd)
library(FOCuS)
library(parallel)

generate_sequence <- function(n = 1000, p = 100, cp = 200, sd = 1, magnitude = 1, dens = .1){

  noise <- matrix(rnorm(p * n, sd = sd), nr = p, nc = n)

  pre_change <- rep(0, p)#rnorm(p, mean = 0)
  #pre_change <- rnorm(p, mean = 0)
  
  post_change <- pre_change + c(sample(c(magnitude, -magnitude), floor(p * dens), replace = T), rep(0, floor(p * (1-dens))))
  noise + cbind(matrix(pre_change, nr = p, nc = cp), matrix(post_change, nr = p, nc = n - cp))

}


# this function takes treshold and known values for pre and post change mean and return an instance of
# an ocd detector ready to do some changepoint monitoring.
ocd_known <- function (thresh, mu0, sd0){
  detector <- ChangepointDetector(dim=nrow(Y_train), method='ocd', beta=1, thresh=thresh)
  detector <- setBaselineMean(detector, mu0)
  detector <- setBaselineSD(detector, sd0)
  
  setStatus(detector, 'monitoring')
}


# this function takes some training observations and return an instance of
# an ocd detector ready to do some changepoint monitoring.
ocd_training <- function (Y_train, thresh){
  detector <- ChangepointDetector(dim=nrow(Y_train), method='ocd', beta=1, thresh=thresh)
  detector <- setStatus(detector, 'estimating')
  for (t in seq_len(ncol(Y_train))) {
    detector <- getData(detector, Y_train[,t])
  }
  setStatus(detector, 'monitoring')
}

ocd_detecting <- function (Y, detector) {
  for (t in seq_len(ncol(Y))){
    detector <- getData(detector, Y[,t])
    if (is.numeric(attr(detector, "status")))
      break
  }
  list(t = t, det = detector)
}
