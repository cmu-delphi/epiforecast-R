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

epiproject.cache.dir = "flusight-natreg-spline-run"
experiment.cache.dir = "flusight-natreg-spline-run"

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

calibrate_forecast_spline = function(forecast, qqs, bins=-1, alpha=0) {
  return(calibrate_forecast(forecast, qqs, bins=bins, alpha=alpha, fit_spline=TRUE))
}

c.calibrations = list('np'=calibrate_forecast_spline,
`beta`=calibrate_forecast_beta,
`none`=calibrate_forecast_null) %>>%
  with_dimnamesnames("Calibration")

s.outer.seasons = s.retro.seasons %>>% with_dimnamesnames("Outer Season")
s.outer.idx.seasons = named_idx_list(s.outer.seasons)

week.window = 3
sswgtfc.forecast.values.file = file.path(experiment.cache.dir,"sswgtfc.forecast.values.rds")
sswgtfc.forecast.values =
  if (file.exists(sswgtfc.forecast.values.file)) {
    readRDS(sswgtfc.forecast.values.file)
  } else {
    map_join_tc(function(so,s,w,g,t,f,cal) {
      week.idxs = w + -week.window:week.window
      week.idxs = week.idxs[1 <= week.idxs & week.idxs <= length(w.retro.model.weeks)]
      qs = swgtf.forecast.quantiles[c(-so,-s),week.idxs,,t,f]
      qs = qs[!is.na(qs)]
      return(cal(swgtmf.forecast.values[[s,w,g,t,"Bin",f]],qs))},
      s.outer.idx.seasons,s.idx.seasons,w.idx.weeks,g.idx.groups,t.idx.targets,f.idx.forecasters,
      c.calibrations)
  }

if (!file.exists(sswgtfc.forecast.values.file)) {
  saveRDS(sswgtfc.forecast.values,sswgtfc.forecast.values.file)
}

gc()

sswgtfc.forecast.evaluations.file = file.path(experiment.cache.dir,"sswgtfc.forecast.evaluations.rds")
sswgtfc.forecast.evaluations =
  if (file.exists(sswgtfc.forecast.evaluations.file)) {
    readRDS(sswgtfc.forecast.evaluations.file)
  } else {
    map_join_tc(get_evaluation,sswgtfc.forecast.values,swgtm.retro.observed.values[,,,,"Bin"],
                no_join(m.forecast.types[["Bin"]]))
  }
mode(sswgtfc.forecast.evaluations) = "numeric"

if (!file.exists(sswgtfc.forecast.evaluations.file)) {
  saveRDS(sswgtfc.forecast.evaluations,sswgtfc.forecast.evaluations.file)
}

swgtfc.indexer = list(all=NULL,smear=-3:3,all=NULL,each=NULL,each=NULL,all=NULL)

options(mc.cores=1L)

get_calibration_weights = function(sswgtfc.forecast.evaluations, indexer.list) {
  swtfc.result = array(0,dim=dim(sswgtfc.forecast.evaluations)[c(1,3,5,6,7)])
  for (s in 1:dim(sswgtfc.forecast.evaluations)[[1]]) {
    swtfc.result[s,,,,] =
      cv_apply(sswgtfc.forecast.evaluations[s,-s,,,,,],
        indexer.list,
        function(train, test) {
          train = R.utils::wrap(train,list(1:5,6))
          train = pmax(train[apply(train,1,function(lst) { all(!is.na(lst))}),],-8)
          degenerate.em.weights = degenerate_em_weights(exp(train))
          return(degenerate.em.weights)
        })[,,,,,,1] %>>% aperm(c(2:4,1))
  }
  dimnames(swtfc.result) <- dimnames(sswgtfc.forecast.evaluations)[c(2,3,5,6,7)]
  return(swtfc.result)
}

options(mc.cores=parallel::detectCores()-1L)

swtfc.weights.file = file.path(experiment.cache.dir,"swtfc.weights.rds")
swtfc.weights = get_calibration_weights(sswgtfc.forecast.evaluations,swgtfc.indexer)
if (!file.exists(swtfc.weights.file)) {
  saveRDS(swtfc.weights,swtfc.weights.file)
}
gc()
swgtfc.forecast.values = apply(sswgtfc.forecast.values,3:7,diag)
swgtfc.forecast.values = unlist(swgtfc.forecast.values,recursive=FALSE)
dim(swgtfc.forecast.values) = dim(sswgtfc.forecast.values)[2:7]
dimnames(swgtfc.forecast.values) = dimnames(sswgtfc.forecast.values)[2:7]
rm(sswgtfc.forecast.values)
rm(sswgtfc.forecast.evaluations)
gc()

swgtf.calibrated.forecasts.file = file.path(experiment.cache.dir,"swgtf.calibrated.forecasts.rds")
swgtf.calibrated.forecasts = 
  if (file.exists(swgtf.calibrated.forecasts.file)) {
    readRDS(swgtf.calibrated.forecasts.file)
  } else {
    map_join_tc("*",swgtfc.forecast.values,swtfc.weights,lapply_variant=lapply) %>>%
      apply(1:5,Reduce,f="+")
  }

if (!file.exists(swgtf.calibrated.forecasts.file)) {
  saveRDS(swgtf.calibrated.forecasts,swgtf.calibrated.forecasts.file)
}

swgtf.calibrated.evaluations.file = file.path(experiment.cache.dir,"swgtf.calibrated.evaluations.rds")
swgtf.calibrated.evaluations =
  if (file.exists(swgtf.calibrated.evaluations.file)) {
    readRDS(swgtf.calibrated.evaluations.file)
  } else {
    map_join_tc(get_evaluation,swgtf.calibrated.forecasts,swgtm.retro.observed.values[,,,,"Bin"],
                no_join(m.forecast.types[["Bin"]]),lapply_variant=lapply)
  }
mode(swgtf.calibrated.evaluations) = "numeric"

if (!file.exists(swgtf.calibrated.evaluations.file)) {
  saveRDS(swgtf.calibrated.evaluations,swgtf.calibrated.evaluations.file)
}

swgtf.original.evaluations = map_join_tc(get_evaluation,swgtmf.forecast.values[,,,,"Bin",],
                                 swgtm.retro.observed.values[,,,,"Bin"],no_join(m.forecast.types[["Bin"]]))
