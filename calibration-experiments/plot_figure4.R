library(dplyr)
library(ggplot2)
library(stringr)

arr = array(0,dim=c(7,4))
evals = readRDS("flusight-natreg-spline-run/sswgtfc.forecast.evaluations.rds")
mode(evals) = "numeric"
evals = apply(evals,3:7,diag)
evals = apply(pmax(evals,-10),c(4,6),mean,na.rm=T)
arr[,1:3] = evals
evals2 = readRDS("flusight-natreg-spline-run/swgtf.calibrated.evaluations.rds")
evals2 = apply(pmax(evals2,-10),4,mean,na.rm=T)
arr[,4] = evals2
dimnames(arr) = list(Target=dimnames(evals)[[1]],Calibration=c(dimnames(evals)[[2]],"Ensemble"))
evals = plyr::adply(arr,.margins=1:2)
evals[["Calibration"]] = as.character(evals[["Calibration"]])
evals[["Calibration"]][evals[["Calibration"]] == "np"] = "Nonparametric"
evals[["Calibration"]][evals[["Calibration"]] == "beta"] = "Parametric"
evals = evals %>%
  inner_join(evals %>%
               filter(Calibration == "none") %>%
               mutate(Baseline=V1) %>%
               select(Target,Baseline)) %>%
  mutate(Delta=V1-Baseline) %>%
  filter(Calibration != "none")

evals[["Calibration"]] = factor(evals[["Calibration"]],levels=c("Nonparametric","Parametric","Ensemble"))

ggplot(evals) +
  geom_bar(aes(x=Calibration,y=Delta,fill=Target),position="dodge",stat="identity") +
  labs(fill="Target",x="Recalibration Method",y="Improvement in Mean Log Score") +
  theme(aspect.ratio=0.75, legend.position="bottom")

ggsave("fig4.png",width=8,height=6)
