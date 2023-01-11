#' Title
#'
#' @param out An object of class FOCuSout
#' @param verbose Prints the cost function at the last iteration as a data.frame.
#' @param ... further argument to be passed to the data.frame print call.
#'
#' @return
#' @export
#'
#' @examples
#' out <- FOCuS(rnorm(100), 18)
#' out

print.FOCuSout <- function(out, verbose = F, ...) {
  
  # if it detected a change
  if (out$changepoint > 0) {
    
    cat("A change was detected at iteration:", out$t, "\n")
    cat("The estimated changepoint is at time:", out$changepoint, "\n")
    
  } else {
    cat("No change detected. Terminating.", "\n")
    
  }
  
  cat("The max of the test statistics was: ", out$maxs[length(out$maxs)], "\n")
  
  
  if (verbose) {
    cat("The cost at the time of termination is: ", "\n")
    print(out$Q1, ...)
  }
  
  
}

print.PiecewiseQuadratic <- function(quad, ...) {
  o <- data.frame("x^2 coef" = sapply(quad, function(q) q$a), 
             "x coef"  = sapply(quad, function(q) q$b), 
             "const" = sapply(quad, function(q) q$c),
             "domain" = sapply(quad, function(q) {
               dom <- ""
               for (i in q$ints)
                 dom <- paste(dom, paste0("(", round(i$l, 2),", ", round(i$u, 2), ")"))
               dom
             }))
  print(o)
}


# set.seed(42)
# y <- c(rnorm(3e5, 1), rnorm(1e4, 0))
# o <- FOCuS(y, 18)
# print(o$Q1)
