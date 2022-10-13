source("simulations/helper_functions.R")

CORES <- 16
SEED <- 45

run_simulation <- function(p, REPS, noise, tlist) {
  print(p)
  grid <- find_grid(0, 21, .01, 1.74)
  data <- mclapply(noise, function (epsilon) c(rep(0, p$changepoint), rep(p$delta, p$N - p$changepoint)) + epsilon, mc.cores = CORES)

  # FOCuS with no pruning costraint
  print("FOCus0")
  res <- mclapply(data, function (y) FOCuS(y, tlist$"FOCuS", mu0 = 0, grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0", est = cp, real = p$changepoint, N = p$N)

  # FoCUS 10
  print("FOCus0 p10")
  res <- mclapply(data, function (y) FOCuS(y,  tlist$"FOCuS 10", mu0 = 0, grid = grid[c(1, 3, 6, 8, 10, 11, 13, 15, 17, 19)], K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0-10p", est = cp, real = p$changepoint, N = p$N))
  #print("page-CUSUM done")

  # Page CUSUM 25
  print("Page 20p")
  res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist$"Page-CUSUM 20", mu0 = 0, grid = grid), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <-  rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-20p", est = cp, real = p$changepoint, N = p$N))

  # Page CUSUM 10
  print("Page 10p")
  res <- mclapply(data, function (y) PageCUSUM_offline(y, tlist$"Page-CUSUM 10", mu0 = 0, grid = grid[c(1, 3, 6, 8, 10, 11, 13, 15, 17, 19)]), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <-  rbind(output, data.frame(sim = 1:REPS, magnitude = p$delta, algo = "Page-10p", est = cp, real = p$changepoint, N = p$N))

  # MOSUM
  print("MOSUM")
  #wins <- unique(10^2 / grid ^ 2) %>% round() # maybe remove the square here?
  wins <- unique(abs(14.4 / grid^2)) %>% round()
  res <- mclapply(data, function (y) MOSUM_offline_kirch(y, tlist$"MOSUM", wins), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = "MOSUM", est = cp, real = p$changepoint, N = p$N))


  return(output)
}


gg <- find_grid(0, 21, .01, 1.74)[11:20] # adding gripoints magnitudes to the simulation
N <- 2e6
sim_grid <- expand.grid(
  N = N,
  changepoint = 1e5,
  # delta = c(.05, .07, seq(.1, 1, by = 0.1), .25, gg),
  delta = unique(c(seq(from = .01, to = .09, by = .005), seq(.1, 1, by = 0.05), seq(0.006, 0.009, by = 0.001), 1.5, 2, .25, gg)) # fine grid
)


load(file = "simulations/pre-change-known/results/tlist.RData")

output_file <- "./simulations/pre-change-known/results/dr_new15.RData"

if (F) {
  NREP <- 100
  set.seed(SEED)
  noise <- lapply(1:NREP, function (i) rnorm(N))
  #run_simulation(sim_grid[10, ], NREP, noise, tlist = tlist)
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP, noise = noise, tlist = tlist))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}


load(output_file)

25






###### plot of the detailed cost comparison (Appendix) ########

maxY <- 1e5 + 50
### plotting

plot_piecewise_quad <- function(x, quad) {
  for (q in quad)
    for (i in q$ints)
      if (i$l < x && i$u >= x)
        return(q$a * x^2 + q$b * x + q$c)
}
plot_piecewise_quad <- Vectorize(plot_piecewise_quad, vectorize.args = "x")


set.seed(32)
y <- c(rnorm(1e5,0 ,.1), rnorm(2e6, 0.48288614, .1))


resF <- FOCuS(y[1:maxY], thres = 18.6, mu0 = 0)


gridP <- find_grid(0, 21, .01, 1.74)
#gridP[41] <- .5

resP <- PageCUSUM_offline(y[1:maxY], thres = 18.25, mu0 = 0, grid = gridP)


plot1 <- ggplot(tibble(mu = -3:3)) +
  stat_function(aes(x = mu), fun = function(x) plot_piecewise_quad(x, quad = resF$Q1), col = 4) +
  theme_idris() +
  xlim(0, 1.5) +
  geom_vline(xintercept = gridP, col = "grey", alpha = .3) +
  xlab(expression(mu)) + ylab(expression(Q[t](mu))) +
  geom_text(aes(x = grid, y = y, label = Q), data = tibble(grid = gridP + .02, y = round(resP$Q, 2) + .1, Q = round(resP$Q, 2)), col = "grey", alpha = .8) +
  geom_text(aes(x = .45, y = round(tail(resF$maxs, 1), 2) + .1, label = round(tail(resF$maxs, 1), 2)), col = 4) +
  geom_vline(xintercept = .482, lty = 2)

plot1
# off grid

set.seed(32)
y2 <- c(rnorm(1e5,0 ,.1), rnorm(2e6, sum(gridP[18:19])/2, .1))


resF2 <- FOCuS(y2[1:maxY], thres = 18.6, mu0 = 0)
resP2 <- PageCUSUM_offline(y2[1:maxY], thres = 18.25, mu0 = 0, grid = gridP)

plot2 <- ggplot(tibble(mu = -3:3)) +
  stat_function(aes(x = mu), fun = function(x) plot_piecewise_quad(x, quad = resF2$Q1), col = 4) +
  theme_idris() +
  xlim(0, 1.5) +
  geom_vline(xintercept = gridP, col = "grey", alpha = .3) +
  xlab(expression(mu)) + ylab(expression(Q[t](mu))) +
  geom_text(aes(x = grid, y = y, label = Q), data = tibble(grid = gridP + .02, y = round(resP2$Q, 2) + .1, Q = round(resP2$Q, 2)), col = "grey", alpha = .8) +
  geom_text(aes(x = .61, y = round(tail(resF2$maxs, 1), 2) + .15, label = round(tail(resF2$maxs, 1), 2)), col = 4) +
  geom_vline(xintercept = sum(gridP[18:19])/2, lty = 2)

plot2


ggsave(ggarrange(plot1, plot2, nrow = 2), filename = "grid-comp.pdf", width = 18, height = 12)