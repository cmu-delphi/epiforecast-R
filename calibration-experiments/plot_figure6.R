library(dplyr)
library(ggplot2)
library(stringr)

arr = array(0,dim=c(7,27,4))
evals = readRDS("flusight-natreg-spline-run/sswgtfc.forecast.evaluations.rds")
mode(evals) = "numeric"
evals = apply(evals,3:7,diag)
evals = apply(pmax(evals,-10),c(4,5,6),mean,na.rm=T)
arr[,,1:3] = evals
evals2 = readRDS("flusight-natreg-spline-run/swgtf.calibrated.evaluations.rds")
evals2 = apply(pmax(evals2,-10),c(4,5),mean,na.rm=T)
arr[,,4] = evals2
dimnames(arr) = list(Target=dimnames(evals)[[1]],Forecaster=dimnames(evals)[[2]],Calibration=c(dimnames(evals)[[3]],"Ensemble"))
dn = dimnames(arr)
arr = sweep(arr[,,c(1,2,4)],c(1,2),arr[,,3],FUN="-")
print(arr)
evals = plyr::adply(arr,.margins=1:3)
evals[["Calibration"]] = as.character(evals[["Calibration"]])
evals[["Calibration"]][evals[["Calibration"]] == "np"] = "Nonparametric"
evals[["Calibration"]][evals[["Calibration"]] == "beta"] = "Parametric"
evals[["Calibration"]] = factor(evals[["Calibration"]],levels=c("Nonparametric","Parametric","Ensemble"))

entropy = function(x, nbins=100) {
  x = x[!is.na(x)]
  x = pmax(pmin(x,1-1/(2*nbins)),1/(2*nbins))
  y = hist(x, breaks=seq(0,1,length.out=nbins+1), plot=FALSE)
  p = nbins*y$counts/length(x)
  return(-sum(p*log(p),na.rm=T)/nbins)
}

quantiles = readRDS("flusight-natreg-spline-run/swgtfc.forecast.quantiles.rds")
quantiles1 = readRDS("flusight-natreg-spline-run/swgtf.calibrated.quantiles.rds")
arre = array(0,dim=c(7,27,4))
arre[,,1:3] = apply(quantiles,4:6,entropy)
arre[,,4] = apply(quantiles1,4:5,entropy)
print(arre[,,4])
dimnames(arre) = dn
arre = sweep(arre[,,c(1,2,4)],c(1,2),arre[,,3],FUN="-")
#print(arre)
ents = plyr::adply(arre,.margins=1:3)
ents[["Calibration"]] = as.character(ents[["Calibration"]])
ents[["Calibration"]][ents[["Calibration"]] == "np"] = "Nonparametric"
ents[["Calibration"]][ents[["Calibration"]] == "beta"] = "Parametric"
ents[["Calibration"]] = factor(ents[["Calibration"]],levels=c("Nonparametric","Parametric","Ensemble"))
ents = ents %>%
  rename(Entropy=V1)
evals = evals %>%
  rename(MLS=V1)
evals = evals %>%
  inner_join(ents)

ggplot(evals %>% filter(str_sub(Target,start=3) == "wk ahead", Calibration=="Ensemble")) +
  geom_point(aes(x=Entropy,y=MLS,color=Target)) +
  labs(color="Target",x="Entropy Improvement",y="Mean Log Score Improvement") +
  theme(aspect.ratio=0.75)

ggsave("fig6.png",width=8,height=6)
