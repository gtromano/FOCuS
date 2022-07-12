setMethod("FOCuS",
          signature(data = "numeric", thres = "numeric"),
          function (datasource, thres, mu0 = NA, training_data = NA, grid = NA, K = Inf) 
          {
            # checks on the data generating function
            if( !is.numeric(datasource))
              stop("datasource must be a vector of observations")
            
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
            if (is.matrix(datasource))
              out <- .FoCUS_mult_offline(datasource, thres, mu0, training_data, grid, K)
            else
              out <- .FoCUS_offline(datasource, thres, mu0, training_data, grid, K)
            out$changepoint <- out$t + out$Q1[[which.max(sapply(out$Q1, function(q) q$max))]]$a * 2
            class(out) <-  c("FOCuSout", class(out))
            class(out$Q1) <- c("PiecewiseQuadratic", class(out))
            
            if (!is.null(out$warning_message))
              warning(out$warning_message)
            
            return(out)
            
          }
)


setMethod("FOCuS",
          signature(data = "matrix", thres = "numeric"),
          function (datasource, thres, mu0 = NA, training_data = NA, grid = NA, K = Inf)
          {
            # checks on the data generating function
            if( !is.numeric(datasource))
              stop("datasource must be a vector or a matrix of observations")

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



            
            # if(is.na(MC_thres[1]))
            #   thresholds <- sapply(1:p, function (j) {
            #     min(
            #       2 * thres + (p * nu + 2 * sqrt(p * nu * thres)),
            #       2 * thres + a *  2 * j * log(p)
            #       )
            #   })
            # else 
            #   thresholds <- MC_thres
            
            

            #warning("Going multivariate!")
            out <- .FoCUS_mult_offline(datasource, thres, mu0, training_data, grid, K)

            out$maxs <- apply(out$stats, 2, max)
            out$sums <- apply(out$stats, 2, sum)
            

            # running the function
            # out <- .FoCUS_offline(datasource, thres, mu0, training_data, grid, K)
            # out$changepoint <- out$t + out$Q1[[which.max(sapply(out$Q1, function(q) q$max))]]$a * 2
            # class(out) <-  c("FOCuSout", class(out))
            # class(out$Q1) <- c("PiecewiseQuadratic", class(out))

            # if (!is.null(out$warning_message))
            #   warning(out$warning_message)

            return(out)

          }
)
