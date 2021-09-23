###   This script produces the illustration of the cost function in Figure 2  ####

source("simulations/helper_functions.R")

########   some helper functions     #########

# for ggplot style plots
piecewise_quad <- function(x, quad) {
  for (q in quad)
    for (i in q$ints)
      if (i$l < x && i$u >= x)
        return(q$a * x^2 + q$b * x + q$c)
}

piecewise_quad <- Vectorize(piecewise_quad, vectorize.args = "x")


plot.Quadratic2 <- function(q, MIN, MAX, ...) {
  for (i in q$ints) {
    curve(.eval.Quadratic(q, x), from = max(i$l, MIN), to = min(i$u, MAX), add = T, col = (-q$a * 2), ...)
  }
}

.eval.Quadratic <- function(q, x) {
  q$a * (x^2) + q$b * x + q$c
}

# get the minimum of a quadratic (sets range of a plot)
getMinimum <- function(q, int) {
  if (q$a == 0) return(list(minim = q$c, at = int$l))
  else {
    at <- -q$b / (2 * q$a)

    if (at <= int$l)
      at <- int$l
    else if (at >= int$u)
      at <- int$u

    minim <- q$a * (at^2) + q$b * at + q$c
    return(list(minim = minim, at = at))
  }
}


########    plotting      #########

NPLOTS <- 4

set.seed(42)
y <- rnorm(1e3)
y[c(500, 501, 502)] <- c(.25, .37, 0.001)
res <- NULL
start <- 498

# getting 3 iterations (running FOCuS 3 times)
for (i in 1:NPLOTS) {

  res[[i]] <- FOCuS_offline(y = y[1:(start + i)], mu0 = 0, thres = Inf)$Q1

}


bigplotlist <- NULL
for (K in 1:NPLOTS) {
  bigplotlist[[K]] <- ggarrange(ggplot(tibble(mu = seq(0, 3, length.out = 1e4))) +
                                  stat_function(aes(x = mu), fun = function(x) piecewise_quad(x, quad = res[[K]]), col = 4) +
                                  theme_idris() +
                                  xlim(0, 1) +
                                  xlab(expression(mu)) +
                                  ylab(expression(Q(mu))))

  plot(NULL, type = "n", xlim = c(0, 2), ylim = c(0, max(piecewise_quad(seq(0, 2, length.out = 100), quad = res[[K]])) + .01),
       xlab = expression(mu), ylab = expression(Q(mu)))
  for (q in res[[K]]) {
    for (i in q$ints) {
      if (i$u < 0) next
      l <- max(0, i$l)
      u <- min(i$u, 5)
      xgr <- seq(l, u, length.out = 100)
      lines(xgr, piecewise_quad(xgr, quad = res[[K]]), col = (500 + K) + 2 * q$a)


      lab <- getMinimum(q, i)
      if (lab$at > 0) text(x = lab$at + 0.1, y = lab$minim + 0.003, labels = (500 + K) + 2 * q$a, col = (500 + K) + 2 * q$a)


    }
  }

}

ggarrange(bigplotlist[[1]], bigplotlist[[2]], bigplotlist[[3]], bigplotlist[[4]], nrow = 1)


#################################################
# ADDITIONAL STUFF


NPLOTS <- 4

set.seed(42)
y <- rnorm(1e3)
y[c(500, 501, 502)] <- c(.25, .37, 0.001)
res <- NULL
start <- 498
for (i in 1:NPLOTS) {

  res[[i]] <- FOCuS_offline(y = y[1:(start + i)], training_data = rnorm(1e3), thres = Inf)$Q1

}


for (K in 1:NPLOTS) {


  plot(NULL, type = "n", xlim = c(-2, 2), ylim = c(0, max(piecewise_quad(seq(-2, 2, length.out = 100), quad = res[[K]])) + .01),
       xlab = expression(mu), ylab = expression(Q(mu)))
  for (q in res[[K]]) {
    for (i in q$ints) {
      #if (i$u < 0) next
      l <- max(-5, i$l)
      u <- min(i$u, 5)
      xgr <- seq(l, u, length.out = 100)
      lines(xgr, piecewise_quad(xgr, quad = res[[K]]), col = (500 + K) + 2 * q$a)


      lab <- getMinimum(q, i)
      if (lab$at > 0) text(x = lab$at + 0.1, y = lab$minim + 0.003, labels = (500 + K) + 2 * q$a, col = (500 + K) + 2 * q$a)


    }
  }

}

