
source("simulations/multivariate/helper_functions.R")

set.seed(42)
Y <- generate_sequence(n = 2000, cp = 1500, magnitude = 2, dens = .05)

# the magic is with 5 magn, 0.05
Y_train <- Y[, 1:1000]
Y_test <- Y[, 1001:2000]

# mu0 estimation phase

mu0hat <- apply(Y_train, 1, mean)
res_focus0 <- FOCuS(Y_test, 12, a = .8, mu0 = mu0hat)
res_focus0$t

ocd_det <- ocd_training(Y_train, thresh = setNames(c(2 * log(500), 179.48, 54.87), c('diag', 'off_d', 'off_s')))
res_ocd <- ocd_detecting(Y_test, ocd_det)
res_ocd$t