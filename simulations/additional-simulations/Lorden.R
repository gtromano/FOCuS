####Computational complexity of Lorden

Lorden = function(mu,z){
n=length(z)
quad.st=rep(NA,n)
P=0
q=0
for(i in 1:n){
  P=P+mu*(z[i]-mu/2) ###Page update
  if(P<0){ #check whether we reset
    P=0
    q=1
  }else{
    q=q+1 ##if not we have an extra quadratic
  }
  quad.st[i]=q
}
return(quad.st)
}

###########################################################
###### comparison on the number of quadratics #############
###########################################################

library(FOCuS)

library(parallel)
N <- 1e6 + 4e4

mus <- c(0.01, 0.05, 0.1, 0.5, 1)

Y <- lapply(1:50, function(i) {
  set.seed(i)
  rnorm(N)
})

df <- data.frame(t = integer(), nq = integer(), algo = character())

for (i in 1:50) {

  t <- c(1, 100, which(1:N %% 1e4 == 0)) # subsampling

  y <- Y[[i]]

  res_foc <- FOCuS_Melk(y, Inf, 0, NaN, NaN)$nquads - 1
  df <- rbind(df, data.frame(t = t, nq = res_foc[t], algo = "FOCuS"))

  for (m in mus) {
    res_Lor <- Lorden(m, y) - 1
    df <- rbind(df, data.frame(t = t, nq = res_Lor[t], algo = paste0("Lor ", m)))
  }

  print(i)

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



library(isotone)
library(tidyverse)



# subsampling
df_summ <- df %>% group_by(algo) %>%
  summarise(t = t, iso = gpava(t, nq)$x)

df_summ <- df_summ %>% group_by(algo, t) %>% summarise(mean = mean(iso))

pava_est <- ggplot(df_summ, aes(x = t, y = mean, col = algo)) +
  geom_line() +
  scale_y_log10() +
  ylab("Ave. Number of Test Statistics/Quadratics") +
  xlim(0, 1e6) +
  scale_color_manual(values = cbPalette) +
  theme_idris()

pava_est
ggsave("simulations/additional-simulations/nquads_comp_iso.png", pava_est, width = 7, height = 4)

###############################################

new_df <- data.frame(t = integer(), fitted = integer(), algo = character())
aaaa <- df %>% pivot_wider(values_from = nq, names_from = algo)

for(algo in unique(df$algo)) {
  bbbb <- Reduce(rbind, aaaa[algo])
  res <- gpava(aaaa$t, bbbb)
  
  new_df <- rbind(new_df, data.frame(t = res$z, fitted = res$x, algo = algo))
}

pava_est <- ggplot(new_df, aes(x = t, y = fitted, col = algo)) +
  geom_line() +
  scale_y_log10() +
  ylab("Avg. Number of quadratics") +
  xlim(0, 1e6) +
  scale_color_manual(values = cbPalette) +
  theme_idris()

pava_est

#############################################################################

cbPalette <- RColorBrewer::brewer.pal(7, "Paired")[c(2, 3, 4, 5, 6, 7)]
final_plot <- ggplot(df, aes(x = t, y = nq, col = algo)) +
  geom_smooth(span = 0.8, se = FALSE) +
  scale_y_log10() +
  ylab("Avg. Number of quadratics") +
  scale_color_manual(values = cbPalette) +
  theme_idris()

final_plot

ggsave("simulations/additional-simulations/nquads_comp.png", final_plot, width = 7, height = 4)
ggsave("simulations/additional-simulations/nquads_comp_loglog.png", final_plot + scale_x_log10(), width = 7, height = 4)
