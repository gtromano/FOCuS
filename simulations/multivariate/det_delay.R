source("simulations/multivariate/helper_functions.R")

CORES <- 16



run_simulation <- function(simu, REPS) {
  print(simu)

  Y <- lapply(1:REPS, function(i) generate_sequence(n = simu$N, cp =  simu$changepoint, magnitude = simu$delta, dens = simu$prop, seed = i))

  
  # FOCuS0 - oracle mean
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    r <- FOCuS(y, foc0_thres, a = .7, mu0 = rep(0, 100))
    ifelse(r$t == -1, simu$N, r$t)
  }, mc.cores = CORES)
  res <- unlist(res)
  output <- data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "FOCuS0 Inf", est = res, real = simu$changepoint, N = simu$N)
  #print("FOCus done")

  # FOCuS0 - estimating mu0
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    mu0hat <- apply(Y_train[[i]], 1, mean)
    r <- FOCuS(y, foc0_est_thres, mu0 = mu0hat)
    ifelse(r$t == -1, simu$N, r$t)
  }, mc.cores = CORES)
  res <- unlist(res)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "FOCuS0 est", est = res, real = simu$changepoint, N = simu$N))
  
  # FOCuS - pre-change unkown
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    r <- FOCuS(y, foc_thres)
    ifelse(r$t == -1, simu$N, r$t)
  }, mc.cores = CORES)
  res <- unlist(res)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "FOCuS", est = res, real = simu$changepoint, N = simu$N))
  
  # ocd oracle
  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    ocd_det <- ocd_known(ocd_thres, rep(0, 100), rep(1, 100))
    r <- ocd_detecting(y, ocd_det)
    r$t
  }, mc.cores = CORES)
  res <- unlist(res)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "ocd Inf", est = res, real = simu$changepoint, N = simu$N))

  res <- mclapply(1:REPS, function(i) {
    y <- Y[[i]]
    y_tr <- Y_train[[i]]
    
    ocd_det <- ocd_training(y_tr, ocd_est_thres)
    r <- ocd_detecting(y, ocd_det)
    r$t
  }, mc.cores = CORES)
  res <- unlist(res)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, magnitude = simu$delta, density = simu$prop, algo = "ocd est", est = res, real = simu$changepoint, N = simu$N))
  
  
  return(output)
}



sim_grid <- expand.grid(
  delta = c(1, .5, .25, .1),  # magnitude of a change
  prop = c(.05, .1, .15, .5),   # proportion of sequences with a change
  changepoint = 200,
  N = 4000
)


# training data for reconstructing the value of mu0
Y_train <- lapply(1:100, function(i) generate_sequence(n = 500, cp = 199, magnitude = 0, dens = 0, seed = 600 + i))


load("simulations/multivariate/thres.RData")

output_file <- "./simulations/multivariate/results/r6.RData"


if (T) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    simu <- sim_grid[i, ]
    return(run_simulation(simu, NREP))
  })
  
  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}

load(output_file)


library(tidyverse)
summ <- outDF %>%
  mutate(det_delay = ifelse(est - real < 0, NA, est - real), false_positive = ifelse(is.na(det_delay), T, F)) %>%
  group_by(density, magnitude, algo) %>%
  summarise(avg_det_delay = mean(det_delay, na.rm = T), fps = sum(false_positive))

summ %>% pivot_wider(names_from = algo, values_from = fps, - avg_det_delay)
(out_table <- summ %>% pivot_wider(names_from = algo, values_from = avg_det_delay, - fps))

write_csv(out_table, file = "./simulations/multivariate/results/summary.csv")
