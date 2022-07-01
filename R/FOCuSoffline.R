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
          function (datasource, thres, a = .8, MC_thres = NA, mu0 = NA, training_data = NA, grid = NA, K = Inf)
          {
            # checks on the data generating function
            if( !is.numeric(datasource))
              stop("datasource must be a vector or a matrix of observations")

            # checks on the threshold
            if( !is.numeric(thres) | thres < 0)
              stop("thres must be a positive numeric")

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
            
            
            nu <- 2

            P1 <- function(j) 2 * thres + (p * nu + 2 * sqrt(p * nu * thres))
            P1 <- Vectorize(P1)
            
            P2 <- function(j) 2 * thres + a *  2 * j * log(p)
            
            
            P3 <- function(j) {
              cj <- qchisq(j/p, nu, lower.tail = F)
              fcj <-  dchisq(cj, nu)
              a * (2 * (thres + log(p))  + j * nu + 2*p*cj*fcj + 2 * sqrt((j * nu + 2*p*cj*fcj) * (thres + log(p))))
            }
            P3 <- Vectorize(P3)
            
            
            
            
            get_threshold <- function(j) min(P1(j),
                                             P2(j))#,
                                             #P3(j))
            get_threshold <- Vectorize(get_threshold)
            
            
            if(is.na(MC_thres[1]))
              thresholds <- get_threshold(1:nrow(datasource))
            else 
              thresholds <- MC_thres
            
            

            #warning("Going multivariate!")
            out <- .FoCUS_mult_offline(datasource, thresholds, mu0, training_data, grid, K)

            out$sortStat <- cumsum(sort(out$maxs[, ifelse(out$t == -1, ncol(datasource), out$t)], decreasing = T))
            out$threshold <- thresholds
            
            if (F) {
              plot(1:p, P1(1:p), type = "l", col = 2)
              lines(1:p, P3(1:p), col = 3)
              lines(1:p, P2(1:p), col = 1)
              lines(1:p, thresholds, col = 5, lty = 2, lwd = 2)
              abline(v = which(out$sortStat >= out$threshold))
            }
            
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
