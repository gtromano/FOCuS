source("simulations/helper_functions.R")

CORES <- 16

run_simulation <- function(p, REPS, seed = 42, tlist) {
  print(p)
  #grid <- find_grid(0, 50, .01, 1.3)

  set.seed(seed)
  #means <- runif(REPS, 1, 10)
  means <- rep(0, REPS)
  data <- lapply(1:REPS, function (k) c(rnorm(1e5 + p$changepoint, means[k]), rnorm(p$N - p$changepoint, means[k] + p$delta)))

  # FOCuS with no pruning constraint
  res <- mclapply(data, function (y) FOCuS(y[(1e5 + 1):length(y)], tlist["FOCuS"], grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS", est = cp, real = p$changepoint, N = p$N)
  #print("FOCus done")

  res <- mclapply(data, function (y) FOCuS(y[(1e5 + 1):length(y)], tlist["FOCuS-t"], training_data = y[1:1e5], grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS-t", est = cp, real = p$changepoint, N = p$N))

  # FOCuS0 1000 estimate
  m <- 1000
  res <- mclapply(data, function (y) FOCuS(y[(1e5 + 1):length(y)], tlist[paste("FOCuS0", m)], mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, 
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = paste("FOCuS0", m), est = cp, real = p$changepoint, N = p$N)) 
  
  m <- 10000
  res <- mclapply(data, function (y) FOCuS(y[(1e5 + 1):length(y)], tlist[paste("FOCuS0", m)], mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, 
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = paste("FOCuS0", m), est = cp, real = p$changepoint, N = p$N)) 
  
  m <- 100000
  res <- mclapply(data, function (y) FOCuS(y[(1e5 + 1):length(y)], tlist[paste("FOCuS0", m)], mu0 = mean(y[1:m]), grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output, 
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = paste("FOCuS0", m), est = cp, real = p$changepoint, N = p$N)) 

  res <- mclapply(1:REPS, function (k) FOCuS(data[[k]][(1e5 + 1):length(data[[k]])], tlist["FOCuS0 Inf"], mu0 = means[k], grid = NA, K = Inf), mc.cores = CORES)
  cp <- sapply(res, function (r) r$t)
  output <- rbind(output,
                  data.frame(sim = 1:REPS, magnitude = p$delta, algo = "FOCuS0 Inf", est = cp, real = p$changepoint, N = p$N))

  return(output)
}


output_file <- "./simulations/pre-change-unknown/results/dr_ukn9_bis.RData"

sim_grid <- expand.grid(
  N = 4e6,
  changepoint = 1e5, # saved on dr_ukn9.RData
  delta = c(- 0.5 ^ seq(5,0, length.out = 10), 0.5 ^ seq(5,0, length.out = 10)) %>% sort
)

load("simulations/pre-change-unknown/tlist.RData")


if (F) {
  NREP <- 100
  outDF <- lapply(seq_len(nrow(sim_grid)), function (i) {
    p <- sim_grid[i, ]
    return(run_simulation(p, NREP, tlist = tlist))
  })

  outDF <- Reduce(rbind, outDF)
  save(outDF, file = output_file)
}


load(output_file)

summary_df <- outDF %>% mutate(
    run_len = if_else(est == -1, N, est),
    no_detection = if_else(est == -1, 1, 0),
    false_alarm = if_else(!no_detection & (est - real < 0), 1, 0),
    det_delay = ifelse((est - real >= 0) & (!false_alarm), est - real, NA),
    true_positive = if_else(!no_detection & !false_alarm,1, 0), # if it's not a missed detection nor it's a false alarm, then it's a true positive
  )


# plug (N - tau) in case we get a missed detection!
summary_df[(summary_df$no_detection) == 1, "det_delay"] <- sim_grid$N[1] - sim_grid$changepoint[1]

grouped <- summary_df %>% group_by(magnitude, algo) %>%
  summarise(no_detection = mean(no_detection), false_alarm = mean(false_alarm), tp_rate = mean(true_positive), det_del = mean(det_delay, na.rm = T))
print(grouped, n = 200)

grouped <- summary_df %>% mutate(magnitude = abs(magnitude)) %>% group_by(magnitude, algo) %>%
  summarise(no_detection = mean(no_detection), false_alarm = mean(false_alarm), tp_rate = mean(true_positive), det_del = mean(det_delay, na.rm = T))
print(grouped, n = 200)


### detection delay ####

generate_labels <- function (dataset, var, XFUNC) {
  library(ggrepel)
  dataset <- dataset[which(dataset[, var] == XFUNC(dataset[, var])), ] %>% mutate(label = as.character(algo))
  dataset
}

data_label <- generate_labels(grouped, "magnitude", max)

cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2,1, 3, 4, 5, 6)]
cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(2, 3, 4, 5, 6)]

detection_delay <-
  ggplot(
    #grouped ,
    #grouped %>% filter(!(algo %in% c("FOCuS-t", "FOCuS0 1000"))),
    grouped %>% filter(!(algo %in% c("FOCuS-t"))),
    aes(
      x = magnitude,
      y = det_del,
      group = algo,
      col = algo
    )
  ) +
    geom_line() +
  scale_color_manual(values = cbPalette) +
  #geom_label_repel(aes(label = label), nudge_x = 100000, force = 10, show.legend = F, data = data_label %>% filter(algo != "FOCuS-t")) +
  xlab("magnitude") +
  ylab("Detection Delay") +
  scale_y_log10() +
  scale_x_log10() +
  theme_idris() #+ theme(legend.position = "none")
detection_delay



##### ratio plot ####

det_del_table <- summary_df %>% filter(magnitude < 2) %>% group_by(magnitude, algo) %>% summarise(dd = mean(det_delay, na.rm = T), no_det = mean(no_detection, na.rm = T), fa = mean(false_alarm, na.rm = T))

summ_table <- pivot_wider(det_del_table[1:3], names_from = algo, values_from = dd) #%>% mutate(FOCuSvPage = FOCuS0 - `Page-20p`, FOCuSvMOSUM = FOCuS0 - MOSUM) %>%  print(n = 100)
summ_table

comp_table <- summ_table %>% mutate(`FOCuS/FOCuS0 1000` = FOCuS / `FOCuS0 1000`,
                                    `FOCuS/FOCuS-t` = `FOCuS`/ `FOCuS-t`,
                                    `FOCuS/FOCuS0 10000` = FOCuS / `FOCuS0 10000`,
                                    `FOCuS/FOCuS0 1e+05` = FOCuS / `FOCuS0 1e+05`,
                                    `FOCuS/FOCuS0 Inf` = FOCuS / `FOCuS0 Inf`)

comp_table <- comp_table[, c(1, 8:12)] %>%
  pivot_longer(names_to = "ratio", values_to = "rval", -magnitude) %>%
  mutate(rval = log(rval)) %>%
  filter(ratio != "FOCuS/FOCuS0 1000")



data_label <- comp_table %>%
  filter(magnitude == max(comp_table$magnitude)) %>%
  mutate(label = as.character(ratio))


cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(1, 4, 5, 6)]
ggplot(comp_table %>% filter(ratio != "FOCuS/FOCuS0 1000")) +
        geom_hline(yintercept = 0, col = "grey", lty = 2) +
        geom_line(aes(x = magnitude, y = rval, col = ratio)) +
        scale_x_log10() +
        geom_label_repel(aes(label = label, x = magnitude, y = rval, col = ratio),min.segment.length = .1, nudge_y = .1, force = 200, nudge_x = 1e3, data = data_label) +
        ylim(-.5, .5) +
        ylab("log ratio") +
        scale_color_manual(values = cbPalette) +
        theme_idris() + theme(legend.position = "none")


# appendix plot
cbPalette <- RColorBrewer::brewer.pal(6, "Paired")[c(1, 5)]
ggplot(comp_table %>% filter(!(ratio %in% c("FOCuS/FOCuS0 1000", "FOCuS/FOCuS0 10000", "FOCuS/FOCuS0 Inf")) )) +
        geom_hline(yintercept = 0, col = "grey", lty = 2) +
        geom_line(aes(x = magnitude, y = rval, col = ratio)) +
        scale_x_log10() +
        geom_label_repel(aes(label = label, x = magnitude, y = rval, col = ratio), min.segment.length = .01, force = 1000, nudge_x = 1, data = data_label %>% filter(!(ratio %in% c("FOCuS/FOCuS0 1000", "FOCuS/FOCuS0 10000", "FOCuS/FOCuS0 Inf")))) +
        ylim(-.05, .4) +
        ylab("log ratio") +
        scale_color_manual(values = cbPalette) +
        theme_idris() + theme(legend.position = "none")
