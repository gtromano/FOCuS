setMethod("FOCuS",
          signature(data = "function", thres = "numeric"),
          function (datasource, thres,  mu0 = NA, grid = NA, K = Inf) 
          {
            # checks on the data generating function
            if(is.function(match.fun(datasource)))
              out_check <- datasource()
            else
              stop("Please provvide a data generating function.")
            if(!is.numeric(out_check)) stop("Data generating function provvided does not return a numeric output.")
            if(length(out_check) != 1) {
              warning("Length of the output from data generating function is greater than 1. Using only the first element.")
              f <- function() return(datasource()[1])
            } else {
              f <- match.fun(datasource)
            }
            
            # checks on the thres
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
            
            
            out <- .FoCUS(datasource, thres, mu0, grid, K)
            out$changepoint <- out$t + out$Q1[[which.max(sapply(out$Q1, function(q) q$max))]]$a * 2
            class(out) <-  c("FOCuSout", class(out))
            class(out$Q1) <- c("PiecewiseQuadratic", class(out))
            
            if (!is.null(out$warning_message))
              warning(out$warning_message)
            
            return(out)
          }
)

  