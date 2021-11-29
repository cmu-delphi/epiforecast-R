library(ggplot2)
library(stringr)

protea_window_evals = readRDS("flusight-natreg-windows-run/swgtfvc.forecast.evaluations.rds")[,,,4:7,"Protea_Cheetah",,]
protea_window_evals = apply(pmax(protea_window_evals,-10),4:6,mean,na.rm=T)
protea_window_evals = plyr::adply(protea_window_evals,.margins=1:3)
protea_window_evals[["Window"]] =
  as.integer(str_sub(protea_window_evals[["Window"]],start=3))
protea_window_evals[["Calibration"]] = as.character(protea_window_evals[["Calibration"]])
protea_window_evals[["Calibration"]][protea_window_evals[["Calibration"]] == "none"] = "None"
protea_window_evals[["Calibration"]][protea_window_evals[["Calibration"]] == "np"] = "Nonparametric"
protea_window_evals[["Calibration"]][protea_window_evals[["Calibration"]] == "beta"] = "Parametric"


ggplot(protea_window_evals) +
  geom_line(aes(x=Window,y=V1,color=Target,linetype=Calibration)) +
  labs(color="Target",linetype="Calibration",y="Mean Log Score",x="Window Size") +
  theme(aspect.ratio=0.75)

ggsave("img.png",width=8,height=6)
