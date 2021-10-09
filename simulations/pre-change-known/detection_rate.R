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
  res <- mclapply(data, function (y) MOSUM_offline_kirch2(y, tlist$"MOSUM", wins), mc.cores = CORES)
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
  delta = unique(c(seq(from = .01, to = .09, by = .005), seq(.1, 1, by = 0.05), .25, gg)) # fine grid
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

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    det_delay = ifelse(est - real > 0, est - real, NA),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (is.na(det_delay)), 1, 0),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positive
  )


det_del_table <- summary_df %>% filter(magnitude > 0, magnitude < 2) %>% group_by(magnitude, algo) %>% summarise(dd = mean(det_delay, na.rm = T), no_det = mean(no_detection, na.rm = T), fa = mean(false_alarm, na.rm = T))
print(det_del_table, n = 100)

det_del_table[is.na(det_del_table)] <- N - 1e5

summ_table <- pivot_wider(det_del_table[1:3], names_from = algo, values_from = dd) #%>% mutate(FOCuSvPage = FOCuS0 - `Page-20p`, FOCuSvMOSUM = FOCuS0 - MOSUM) %>%  print(n = 100)
summ_table

#write_excel_csv(summ_table, file = "table_out.csv")

algs_to_compare <- names(summ_table)[-1]
possible_comparisons <- expand.grid(alg1 = algs_to_compare,  alg2 = algs_to_compare) # this is the ratios to calculate
plot_mat <- matrix(seq_len(nrow(possible_comparisons)), 5, 5) # this is to map the plots with the comparisons
to_plot <- upper.tri(plot_mat)

# now we make a dataframe with all the ratios
comp_table <- summ_table %>% select(magnitude)

for (i in plot_mat[to_plot]) {
  comp <- possible_comparisons[i, ]
  comp_name <- paste0(as.character(comp$alg1), "v", as.character(comp$alg2))
  # comp_table[, comp_name] <- summ_table[, comp$alg1] / summ_table[, comp$alg2] # regular ratio
  comp_table[, comp_name] <- log(summ_table[, comp$alg1] / summ_table[, comp$alg2]) # log ratio
}

comp_table <- comp_table %>% as.data.frame
to_plot2 <- lower.tri(plot_mat, diag = T)

# this for cycle construct the matrix of plots that will then be arranged by ggarrange
ggrid <- find_grid(0, 21, .01, 1.74)

##### version with the residuals #########
plot_list <- NULL
for (i in seq_along(plot_mat)) {

  if(to_plot2[i]) { # this is if we have to actually make a visualization

    comp <- possible_comparisons[t(plot_mat)[i], ]
    if (comp$alg1 == comp$alg2) { # case one, the title plot
      plot_list[[length(plot_list) + 1]] <- ggplot(comp_table) +
        annotate("text", x = 1, y = 1, size = 6, vjust = .5, hjust = 0.5,  label = comp$alg1) +
        theme_idris() +
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank())
    } else {
      comp_name <- paste0(as.character(comp$alg1), "v", as.character(comp$alg2)) # case two, the actual visualization of the ratio
      p <- ggplot(comp_table) +
        geom_hline(yintercept = 0, col = "grey", lty = 2) +
        geom_point(aes(x = g, y = 0), data = data.frame(g = ggrid), alpha = .8, pch = 1) +
        geom_point(aes(x = g, y = 0), data = data.frame(g = ggrid[c(1, 3, 6, 8, 10, 11, 13, 15, 17, 19)]), alpha = .8, pch = 16) +
        geom_line(aes(x = magnitude, y = comp_table[, comp_name])) +
        scale_x_log10() +
        ylim(-.5, .5) +
        theme_idris() +
        theme(axis.title=element_blank())
      plot_list[[length(plot_list) + 1]] <- ggarrange(p)
    }

  } else { # to plot some instograms
      comp <- possible_comparisons[plot_mat[i], ]

      comp_name <- paste0(as.character(comp$alg1), "v", as.character(comp$alg2)) # case two, the actual visualization of the ratio
      p2 <- ggplot(comp_table) +
        geom_histogram(aes(x = comp_table[, comp_name])) +
        geom_vline(xintercept = 0, lty = 2, col = "grey") +
        xlim(-max(abs(comp_table[, comp_name])), max(abs(comp_table[, comp_name]))) +
        #xlim(.5, 1.5) +
        theme_idris() +
        theme(axis.title=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
      plot_list[[length(plot_list) + 1]] <- ggarrange(p2)
  }

}

ggarrange(plotlist=plot_list, ncol = 5, nrow = 5)








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