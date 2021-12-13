library("pipeR")
devtools::load_all("/home/ec2-user/shared/Delphi/epiforecast-R/epiforecast")

options(mc.cores=parallel::detectCores()-1L)
## options(mc.cores=parallel::detectCores()-3L)

## different location naming schemes:
fluview.location.epidata.names = c("nat", paste0("hhs",1:10))
fluview.location.spreadsheet.names = c("US National", paste0("HHS Region ",1:10))

s.retro.seasons = seq.int(2010L,2018L) %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
w.retro.model.weeks = (0:28) %>>%
  stats::setNames(paste0("MW",.+40)) %>>%
  with_dimnamesnames("Model Week")
g.epigroups = fluview.location.spreadsheet.names %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Location")
target_trajectory_preprocessor = flusight2016ilinet_target_trajectory_preprocessor
t.target.specs = flusight2016.target.specs %>>%
  with_dimnamesnames("Target")
m.forecast.types = flusight2016.proxy.forecast.types %>>%
  with_dimnamesnames("Type")

experiment.cache.dir = "flusight-natreg-spline-run"
experiment.cache.dir = "test"

f.forecasters = c("CUBMA","CU_EAKFC_SEIRS","CU_EAKFC_SIRS","CU_EKF_SEIRS",
  "CU_EKF_SIRS","CU_RHF_SEIRS","CU_RHF_SIRS","Delphi_BasisRegression",
  "Delphi_EmpiricalFutures","Delphi_EmpiricalTrajectories",
  "Delphi_ExtendedDeltaDensity","Delphi_MarkovianDeltaDensity",
  "Delphi_Uniform","FluOutlook_Mech","FluOutlook_MechAug","FluX_ARLR",
  "FluX_LSTM","LANL_DBMplus","Protea_Cheetah","Protea_Kudu","Protea_Springbok",
  "ReichLab_kcde_backfill_none","ReichLab_kcde_backfill_post_hoc",
  "ReichLab_kde","ReichLab_sarima_seasonal_difference_FALSE",
  "ReichLab_sarima_seasonal_difference_TRUE","UA_EpiCos") %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Forecaster")

seasonModelWeekToEpiweek = function(season, model.week, first.week=40, owning.wday=3) {
  n.weeks = lastWeekNumber(season,owning.wday)
  year = ifelse(model.week > n.weeks-first.week, season+1L, season)
  week = ifelse(model.week > n.weeks-first.week, model.week-n.weeks, model.week) + first.week
  return(c(year, week))
}

load_file = function(season, model.week, forecaster) {
  ew = seasonModelWeekToEpiweek(season,model.week)
  filename = sprintf("cdc-flusight-ensemble/model-forecasts/component-models/%s/EW%02d-%d-%s.csv",forecaster,ew[[2]],ew[[1]],forecaster)
  result = tryCatch({
    read.csv(filename)
  }, error = function(e) { 
    NA
  })
  return(result)
}

swf.forecast.dfs = map_join(load_file,
  s.retro.seasons,w.retro.model.weeks,f.forecasters)

df_to_val = function(forecast.df, loc, target, type) {
  if(length(forecast.df) == 1 && is.na(forecast.df)) {
    return(NA)
  }
  names(forecast.df) = tolower(names(forecast.df))
  mask = which(forecast.df[["location"]]==loc & forecast.df[["target"]]==target[["Target"]]
                & forecast.df[["type"]] == type[["Type"]])
  return(forecast.df[mask,"value"])
}

swgtmf.forecast.values = map_join(df_to_val,
  swf.forecast.dfs,g.epigroups,t.target.specs,m.forecast.types) %>>%
  aperm(c(1:2,4:6,3))

saveRDS(swgtmf.forecast.values,file=file.path(experiment.cache.dir,"swgtmf.forecast.values.rds"))
