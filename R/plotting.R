# this function plots a quadratic between MIN and MAX
plot.Quadratic <- function(q, MIN, MAX, ...) {
  for (i in q$ints) {
    curve(.eval.Quadratic(q, x), from = max(i$l, MIN), to = min(i$u, MAX), add = T, ...)
  }
}

.eval.Quadratic <- function(q, x){
  q$a * (x ^ 2) + q$b * x + q$c
}
