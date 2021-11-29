library(dplyr)
library(ggplot2)
library(stringr)
library(entropy)

arr = array(0,dim=c(7,27,2))
evals = readRDS("flusight-natreg-spline-run/sswgtfc.forecast.evaluations.rds")
mode(evals) = "numeric"
evals = apply(evals[,,,,,,3],3:6,diag)
evals = apply(pmax(evals,-10),c(4,5),mean,na.rm=T)
arr[,,1] = evals
evals2 = readRDS("flusight-natreg-spline-run/swgtf.calibrated.evaluations.rds")
evals2 = apply(pmax(evals2,-10),c(4,5),mean,na.rm=T)
arr[,,2] = evals2
dimnames(arr) = list(Target=dimnames(evals)[[1]],Forecaster=dimnames(evals)[[2]],Calibration=c("Original","Calibrated"))
dn = dimnames(arr)
evals = plyr::adply(arr,.margins=1:2)
evals = evals %>%
  rename(OriMLS=Original,CalMLS=Calibrated)

entropy = function(x, nbins=100) {
  x = x[!is.na(x)]
  x = pmax(pmin(x,1-1/(2*nbins)),1/(2*nbins))
  y = hist(x, breaks=seq(0,1,length.out=nbins+1), plot=FALSE)
  p = nbins*y$counts/length(x)
  return(-sum(p*log(p),na.rm=T)/nbins)
}

quantiles = readRDS("flusight-natreg-spline-run/swgtfc.forecast.quantiles.rds")
quantiles1 = readRDS("flusight-natreg-spline-run/swgtf.calibrated.quantiles.rds")
arre = array(0,dim=c(7,27,2))
arre[,,1] = apply(quantiles[,,,,,3],4:5,entropy)
arre[,,2] = apply(quantiles1,4:5,entropy)
dimnames(arre) = dn
ents = plyr::adply(arre,.margins=1:2)
ents = ents %>%
  rename(OriEnt=Original,CalEnt=Calibrated)
evals = evals %>%
  inner_join(ents)

unifs = array(0,dim=c(100))
for (i in 1:100) {
  unifs[[i]] = entropy3(runif(9*29*11,0,1))
}
qs = quantile(unifs,c(0.05,0.95))

ggplot(evals %>% filter(str_sub(Target,start=3) == "wk ahead")) +
  geom_segment(aes(x=OriEnt,y=OriMLS,xend=CalEnt,yend=CalMLS,color=Target),
    arrow = arrow(length=unit(0.01,"npc"))) +
  geom_point(aes(x=OriEnt,y=OriMLS,color=Target), show.legend=FALSE, size=0.8) +
  geom_vline(xintercept=qs,linetype="dotted") +
  labs(color="Target",x="Entropy",y="Mean Log Score") +
  theme(aspect.ratio=0.75)

ggsave("fig7.png",width=8,height=6)
