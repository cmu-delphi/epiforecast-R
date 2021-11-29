library(dplyr)
library(ggplot2)
library(stringr)

arr = array(0,dim=c(4,3))
e1 = readRDS("flusight-natreg-spline-run/swgt.ensemble.evaluations.rds")
e1 = apply(pmax(e1,-10),4,mean,na.rm=T)
arr[,1] = e1[4:7]
e2 = readRDS("flusight-natreg-spline-run/swgt.calibrated.ensemble.evaluations.rds")
e2 = apply(pmax(e2,-10),4,mean,na.rm=T)
arr[,2] = e2[4:7]
e3 = readRDS("flusight-natreg-spline-run/swgt.ensemble.calibrated.evaluations.rds")
e3 = apply(pmax(e3,-10),4,mean,na.rm=T)
arr[,3] = e3[4:7]
dimnames(arr) = list(Target=1:4,Calibration=c("None","E-C","C-E"))
evals = plyr::adply(arr,.margins=1:2) %>%
  rename(MLS = V1) %>%
  mutate(Target = as.numeric(Target))

xs = c(14,21,28)
ys = e2[5:7]
ms = (e1[5:7]-e1[4:6])/7
bs = e1[5:7] - xs*ms
xends = (ys-bs)/ms
seg_df = tibble(xs=xs,ys=ys,xends=xends)

ggplot(evals) +
  geom_line(aes(x=7*Target,y=MLS,color=Calibration,group=Calibration)) +
  geom_point(aes(x=7*Target,y=MLS,color=Calibration)) +
  geom_segment(aes(x=xs,y=ys,xend=xends,yend=ys,color="E-C"),data=seg_df,linetype="dotted") +
  labs(color="Calibration",x="Days Ahead",y="Mean Log Score") +
  theme(aspect.ratio=0.75)

ggsave("fig8.png",width=8,height=6)
