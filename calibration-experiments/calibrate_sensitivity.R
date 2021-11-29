library("pipeR")
devtools::load_all("../epiforecast")

options(mc.cores=parallel::detectCores()-1L)

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
v.windows = (0:10) %>>%
  stats::setNames(paste0("WW",.)) %>>%
  with_dimnamesnames("Window")

epiproject.cache.dir = "flusight-natreg-spline-run"
experiment.cache.dir = "flusight-natreg-windows-run"

swgt.retro.observed.multivals.file = file.path(epiproject.cache.dir,"swgt.retro.observed.multivals.rds")
swgt.retro.observed.multivals = readRDS(swgt.retro.observed.multivals.file)[8:16,6:34,,]

swgtm.retro.observed.values.file = file.path(epiproject.cache.dir,"swgtm.retro.observed.values.rds")
swgtm.retro.observed.values = readRDS(swgtm.retro.observed.values.file)[8:16,6:34,,,]

map_join_tc = function(f, ...) {
  map_join(function(...) { tryCatch(f(...), error = function(e) { NA })}, ...)
}

swgtmf.forecast.values = readRDS(file.path(experiment.cache.dir,"swgtmf.forecast.values.rds"))
swgtf.forecast.quantiles = map_join_tc(get_observed_quantile_cdc,
  swgtmf.forecast.values[,,,,"Bin",],swgtm.retro.observed.values[,,,,"Bin"])
mode(swgtf.forecast.quantiles) = "numeric"

named_idx_list = function(lst) {
  return(with_dimnames(1:length(lst),dimnames(lst)))
}
s.idx.seasons = named_idx_list(s.retro.seasons)
w.idx.weeks = named_idx_list(w.retro.model.weeks)
g.idx.groups = named_idx_list(g.epigroups)
t.idx.targets = named_idx_list(t.target.specs)
f.forecasters = dimnames(swgtmf.forecast.values)[[6]] %>>% stats::setNames(.) %>>% with_dimnamesnames("Forecaster")
f.idx.forecasters = named_idx_list(f.forecasters)

c.calibrations = list('np'=calibrate_forecast,
`beta`=calibrate_forecast_beta,
`none`=calibrate_forecast_null) %>>%
  with_dimnamesnames("Calibration")

swgtfvc.forecast.values.file = file.path(experiment.cache.dir,"swgtfvc.forecast.values.rds")
swgtfvc.forecast.values =
  if (file.exists(swgtfvc.forecast.values.file)) {
    readRDS(swgtfvc.forecast.values.file)
  } else {
    map_join_tc(function(s,w,g,t,f,v,cal) {
      week.idxs = w + -v:v
      week.idxs = week.idxs[1 <= week.idxs & week.idxs <= length(w.retro.model.weeks)]
      qs = swgtf.forecast.quantiles[-s,week.idxs,,t,f]
      qs = qs[!is.na(qs)]
      return(cal(swgtmf.forecast.values[[s,w,g,t,"Bin",f]],qs))},
      s.idx.seasons,w.idx.weeks,g.idx.groups,t.idx.targets,f.idx.forecasters,
      v.windows,c.calibrations)
  }

if (!file.exists(swgtfvc.forecast.values.file)) {
  saveRDS(swgtfvc.forecast.values,swgtfvc.forecast.values.file)
}

rm(swgtmf.forecast.values)
rm(swgtf.forecast.quantiles)
gc()

swgtfvc.forecast.evaluations.file = file.path(experiment.cache.dir,"swgtfvc.forecast.evaluations.rds")
swgtfvc.forecast.evaluations =
  if (file.exists(swgtfvc.forecast.evaluations.file)) {
    readRDS(swgtfvc.forecast.evaluations.file)
  } else {
    map_join_tc(get_evaluation,swgtfvc.forecast.values,swgtm.retro.observed.values[,,,,"Bin"],
                no_join(m.forecast.types[["Bin"]]),lapply_variant=lapply)
  }
mode(swgtfvc.forecast.evaluations) = "numeric"

if (!file.exists(swgtfvc.forecast.evaluations.file)) {
  saveRDS(swgtfvc.forecast.evaluations,swgtfvc.forecast.evaluations.file)
}
