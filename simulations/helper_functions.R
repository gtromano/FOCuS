library(tidyverse)
library(FOCuS)
library(parallel)
library(ggpubr)
library(mosum)


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

MOSUMwrapper <- function (x, bandw, thres) {
  m <- mosum(x, G = bandw)
  
  cp <- which(m$stat > thres)[1]
  t <- cp + bandw
  
  if (is.na(cp)) {
    cp <- -1
    t <- length(x)
  }
  
  return(list(cp = cp, t = t, out = m))
}

