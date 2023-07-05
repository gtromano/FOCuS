#' A Simple High-dimensional Extension of FOCuS
#' 
#' A simple implementation of the \link{FOCuS} to deal with multivariate, high-dimensional sequences. 
#' This is achieved by taking both the maximum and sum of the \link{FOCuS} statistics across data streams. 
#'
#' @param datasource A matrix of dimension \code{p x n}, containing \code{p} sequences of length \code{n}. 
#' @param thres The Thresholds for the max and the sum of the FOCuS sequences. This needs to be in the form of positive numeric vector or length 2.
#' @param ... Other parameters as described in \link{FOCuS}. 
#'
#' @return Returns an s3 object of class FOCuSout where:
#' \describe{
#' \item{\code{$t}}{is the stopping time at the detection.}
#' \item{\code{$stats}}{A matrix containing the FOCuS statistics from each data stream.}
#' \item{\code{$maxs}}{the max across the various data stream at each iteration.}
#' \item{\code{$sums}}{the sum across the various data stream at each iteration.}
#' }
#' @export
#'
#' @examples

setGeneric(
  "FOCuS_multi_JMLR",
  def = function(datasource, thres, ...) standardGeneric("FOCuS_multi_JMLR")
)

setMethod("FOCuS_multi_JMLR",
          signature(datasource = "matrix", thres = "numeric"),
          function (datasource, thres, mu0 = NA, training_data = NA, grid = NA, K = Inf)
          {
            # checks on the data generating function
            if( !is.numeric(datasource))
              stop("datasource must be a matrix of observations")
            
            # checks on the threshold
            if( length(thres) < 2)
              stop("thres must be a positive numeric vector or length 2, for the max and sum thresholds")
            
            # checks on the mu0
            if(!is.na(mu0[1])) {
              if(!is.numeric(mu0) | length(mu0) != nrow(datasource))
                stop(paste("mu0 must be a numeric vector of length", nrow(datasource)))
            }
            
            
            # checks on the grid
            if(!is.na(grid[1]))
              if(!is.numeric(grid))
                stop("mu0 must be a numeric value")
            
            # checks on the K
            if(!is.na(K))
              if(!is.numeric(K) | K <= 0)
                stop("K must be a positive numeric")
            
            p <- nrow(datasource)
            
            
            
            #warning("Going multivariate!")
            out <- .FoCUS_mult_offline(datasource, thres, mu0, training_data, grid, K)
            
            out$maxs <- apply(out$stats, 2, max)
            out$sums <- apply(out$stats, 2, sum)
            

            return(out)
            
          }
)