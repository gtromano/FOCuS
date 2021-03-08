library(FOCuS)
library(parallel)


set.seed(42)
Y <- c(rnorm(1e6), rnorm(250, -0.4))
res1 <- FOCuS_melk(Y, thres = 15, mu0 = 0, grid = NA, K = Inf)


test_that("Melk algorithm - no pruning - no robust, error", {
  expect_equal(res1$t, 1000166)
  expect_equal(res1$maxs[res1$t], 16.233536363712385508506486075930297374725341796875)
})



res2 <- FOCuS_melk(Y, thres = 15, mu0 = 0, grid = FOCuS::find_grid(0, 10, .1, 2), K = Inf)

test_that("Melk algorithm - pruning - no robust, error", {
  expect_equal(res2$t, 1000166)
  expect_equal(res2$maxs[res1$t], 16.233536363712385508506486075930297374725341796875)
  expect_equal(length(res2$Qleft), 7)
})
