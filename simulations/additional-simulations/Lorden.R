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
#
# ##find expected number of quadratics as we vary mu.
# mus=c(1:10/100,0.1+1:10/25,0.5+1:5/10)
#
# K=50
# maxqs=matrix(NA,nrow=K,ncol=length(mus))
# meanqs=matrix(NA,nrow=K,ncol=length(mus))
#
# for(k in 1:K){
# set.seed(k)
# z=rnorm(1e6)
# for(j in 1:length(mus)){
#   out=Lorden(mus[j],z)
#   maxqs[k,j]=max(out)
#   meanqs[k,j]=mean(out)
# }
# cat(".")
# if(k%%10==0) cat("\n")
# }
#
#
# savepdf <- function(file, width=16, height=10)
# {
#   fname <- paste(file,".pdf",sep="")
#   pdf(fname, width=width/2.54, height=height/2.54,
#       pointsize=10)
#   par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
# }
#
# #savepdf("Quadratics",width=10,height=6)
# ##
# par(mfrow=c(1,2))
# plot(mus,apply(meanqs,2,mean),log="y",ylab="Ave. Number of Quadratics",xlab="mu",type="l")
# plot(mus,apply(maxqs,2,mean),log="y",ylab="Max Number of Quadratics",xlab="mu",type="l")
# #dev.off()

###########################################################
###### comparison on the number of quadratics #############
###########################################################

library(FOCuS)

library(parallel)
N <- 1e4

mus <- c(0.01, 0.05, 0.1, 0.5, 1)

Y <- lapply(1:20, function(i) {
  set.seed(i)
  rnorm(N)
})

totsim <- mclapply(Y, function (y) {
  res_foc <- FOCuS_Melk(y, Inf, 0, NaN, NaN)
  df <- data.frame(t = 1:N, nq = res_foc$nquads, algo = "FOCuS")

  for (m in mus) {
    res_Lor <- Lorden(m, y)
    df <- rbind(df,
                data.frame(t = 1:N, nq = res_Lor, algo = paste0("Lor ", m)))

  }
  df
}, mc.cores = 4)

df <- totsim[[1]]
for (sim in totsim[-1]) {
  df <- rbind(df, sim)
}

rm(totsim, Y)

library(ggplot2)
ggplot(df, aes(x = t, y = nq, col = algo)) +
  stat_summary(fun = mean, geom="line", alpha = 0.8) +
  scale_y_log10() +
  ylab("Number of quadratics") +
  theme_bw()
