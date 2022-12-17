library(tidyverse)
library(FOCuS)
library(parallel)
library(ggpubr)
library(mosum)

hugefonts <- function (size = 18) {
 theme(text=element_text(size=size), #change font size of all text
        axis.text=element_text(size=size), #change font size of axis text
        axis.title=element_text(size=size), #change font size of axis titles
        plot.title=element_text(size=size), #change font size of plot title
        legend.text=element_text(size=size), #change font size of legend text
        legend.title=element_text(size=size)) #change font size of legend title
}
theme_idris <- function() {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey20"),
    panel.border =  element_rect(fill = NA,
                                 colour = "grey20")
  )
}


### I could not find any true online implementation,
### so this is just a wrapper to pick the first point that passes the threshold

# MOSUMwrapper <- function (x, bandw, thres) {
#   m <- mosum(x, G = bandw)
#
#   cp <- which(m$stat > thres)[1]
#   t <- cp + bandw
#
#   if (is.na(cp)) {
#     cp <- -1
#     t <- length(x)
#   }
#
#   return(list(cp = cp, t = t, out = m))
# }
#

# .H <- function (x, mu0 = 0) x - mu0
#
# .mosum_stat <- function (Y, n, w, FUN) {
#     const <- 1 / sqrt(w)
#     abs(sum(FUN(Y[(n-w):n]))) * const
# }
#
# .v_min <- function (v, ...) sapply(v, function (a) min(a, ...))
#
# MOSUM_offline_kirch <- function (Y, threshold, W, FUN = .H) {
#   stat <- sapply((max(W)+1):length(Y), function (n){
#     Q <- sapply(W, function(w) .mosum_stat(Y, n, w, FUN))
#     max(Q)
#   })
#
#   stat <- c(rep(0, max(W)), stat) # adding w zeroes at the start of the statistic
#
#   cp <- which(stat >= threshold)[1]
#   cp <- ifelse(is.na(cp), -1, cp)
#   return(list(cp = cp, maxs = stat))
# }
#
# # more efficient version
# MOSUM_offline_kirch <- function (Y, threshold, W, FUN = .H) {
#
#   W <- .v_min(W, length(Y)-1)
#
#   tot <- lapply(W, function (w) {
#     stat <- sapply((w+1):length(Y), function (n){
#       .mosum_stat(Y, n, w, FUN)
#     })
#     c(rep(0, w), stat) # adding w zeroes at the start of the statistic
#   })
#
#   stat <- apply(Reduce(cbind, tot), 1, FUN = max)
#
#   cp <- which(stat >= threshold)[1]
#   cp <- ifelse(is.na(cp), -1, cp)
#   return(list(cp = cp, maxs = stat))
# }
#

#MOSUM_offline_kirch(rnorm(1e3), threshold = Inf, wins)
