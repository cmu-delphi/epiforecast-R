library(dplyr)
library(ggplot2)
library(stringr)
library(entropy)

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
arr = sweep(arr[,,c(1,2,4)],c(1,2),arr[,,3],FUN=">")
arr = apply(arr,c(1,3),mean)
print(arr)
evals = plyr::adply(arr,.margins=1:2)
evals[["Calibration"]] = as.character(evals[["Calibration"]])
evals[["Calibration"]][evals[["Calibration"]] == "np"] = "Nonparametric"
evals[["Calibration"]][evals[["Calibration"]] == "beta"] = "Parametric"
evals[["Calibration"]] = factor(evals[["Calibration"]],levels=c("Nonparametric","Parametric","Ensemble"))

entropy = function(x, nbins=100) {
  x = x[!is.na(x)]
  y = discretize(x, numBins=nbins, r=c(0,1))
  return(entropy.MillerMadow(y) - log(nbins))
}
quantiles = readRDS("flusight-natreg-spline-run/swgtfc.forecast.quantiles.rds")
quantiles1 = readRDS("flusight-natreg-spline-run/swgtf.calibrated.quantiles.rds")
arre = array(0,dim=c(7,27,4))
arre[,,1:3] = apply(quantiles,4:6,entropy)
arre[,,4] = apply(quantiles1,4:5,entropy)
dimnames(arre) = dn
arre = sweep(arre[,,c(1,2,4)],c(1,2),arre[,,3],FUN=">")
arre = apply(arre,c(1,3),mean)
print(arre)
ents = plyr::adply(arre,.margins=1:2)
ents[["Calibration"]] = as.character(ents[["Calibration"]])
ents[["Calibration"]][ents[["Calibration"]] == "np"] = "Nonparametric"
ents[["Calibration"]][ents[["Calibration"]] == "beta"] = "Parametric"
ents[["Calibration"]] = factor(ents[["Calibration"]],levels=c("Nonparametric","Parametric","Ensemble"))
ents[["Type"]] = "Entropy"
evals[["Type"]] = "Mean Log Score"
evals = bind_rows(evals,ents)
evals[["Type"]] = factor(evals[["Type"]],levels=c("Mean Log Score","Entropy"))

ggplot(evals) +
  geom_bar(aes(x=Calibration,y=V1,fill=Target),position="dodge",stat="identity") +
  facet_wrap(~Type) +
  labs(fill="Target",x="Recalibration Method",y="Proportion of Models Improved") +
  theme(aspect.ratio=0.75)

ggsave("fig5.png",width=8,height=5)
