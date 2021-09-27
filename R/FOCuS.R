#'  Fast Online Changepoint Detection via Functional Pruning CUSUM statistics
#'
#' @param datasource Either a data generating function, for an online analysis, or a vector of observations, for offline testing. See \strong{Details}.
#' @param thres The threshold for a detection.
#' @param ... Other additional arguments passed to the method. 
#' @param mu0 The value of the pre-change mean, if known. Defaulting to \code{NA}. When \code{NA}, the pre-change-mean unknown recursion will be employed. Pre-change mean is therefore estimated iteratively.
#'
#' @return
#' @export
#' @details FOCuS employs S4 method dispatch based on the data source to perform either an online or offline analysis. For running the algorithm online, \code{datasource} requires an object of class \code{'function'}. For running the algorithm offline, useful for testing purposes, the method expects an object of class \code{'vector'}.
#'
#' @examples
setGeneric(
  "FOCuS",
  def = function(datasource, thres, ...) standardGeneric("FOCuS")
)
