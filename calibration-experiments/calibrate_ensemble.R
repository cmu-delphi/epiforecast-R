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
swgtf.forecast.evaluations = map_join_tc(get_evaluation,
  swgtmf.forecast.values[,,,,"Bin",],swgtm.retro.observed.values[,,,,"Bin"],no_join(m.forecast.types[["Bin"]]))

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

swgtf.indexer = list(loo=NULL,smear=-3:3,all=NULL,each=NULL,all=NULL)
mode(swgtf.forecast.evaluations) = "numeric"

swtf.weights.file = file.path(experiment.cache.dir,"swtf.ensemble.weights.rds")
swtf.weights =
  if (file.exists(swtf.weights.file)) {
    readRDS(swtf.weights.file)
  } else {
    cv_apply(swgtf.forecast.evaluations,swgtf.indexer,
      function(train, test) {
        train = R.utils::wrap(train,list(1:4,5))
        train = pmax(train[apply(train,1,function(lst) { all(!is.na(lst))}),],-8)
        degenerate.em.weights = degenerate_em_weights(exp(train))
        return(degenerate.em.weights)
      })[,,,,,1] %>>% aperm(c(2:4),1)
dimnames(swft.weights) = dimnames(swgtf.forecast.evaluations)[c(1,2,4,5)]

if (!file.exists(swtf.weights.file)) {
  saveRDS(swtf.weights,swtf.weights.file)
}

swgt.ensemble.forecasts.file = file.path(experiment.cache.dir,"swgt.ensemble.forecasts.rds")
swgt.ensemble.forecasts =
  if (file.exists(swgt.ensemble.forecasts.file)) {
    readRDS(swgt.ensemble.forecasts.file)
  } else {
    map_join_tc("*",swgtmf.forecast.values[,,,,"Bin",],swtf.weights,lapply_variant=lapply) %>>%
      apply(1:4,Reduce,f="+")
  }

if (!file.exists(swgt.ensemble.forecasts.file)) {
  saveRDS(swgt.ensemble.forecasts,swgt.ensemble.forecasts.file)
}

swgt.ensemble.evaluations.file = file.path(experiment.cache.dir,"swgt.ensemble.evaluations.rds")
swgt.ensemble.evaluations =
  if (file.exists(swgt.ensemble.evaluations.file)) {
    readRDS(swgt.ensemble.evaluations.file)
  } else {
    map_join_tc(get_evaluation,swgt.ensemble.forecasts,swgtm.retro.observed.values[,,,,"Bin"],
                no_join(m.forecast.types[["Bin"]]),lapply_variant=lapply)
  }
mode(swgt.ensemble.evaluations) = "numeric"

if (!file.exists(swgt.ensemble.evaluations.file)) {
  saveRDS(swgt.ensemble.evaluations,swgt.ensemble.evaluations.file)
}

swgt.ensemble.quantiles = map_join_tc(get_observed_quantile_cdc,
  swgt.ensemble.forecasts,swgtm.retro.observed.values[,,,,"Bin"])
mode(swgt.ensemble.quantiles) = "numeric"

sswgtc.ensemble.forecast.values.file = file.path(experiment.cache.dir,"sswgtc.ensemble.forecast.values.rds")
sswgtc.ensemble.forecast.values =
  if (file.exists(sswgtc.ensemble.forecast.values.file)) {
    readRDS(sswgtc.ensemble.forecast.values.file)
  } else {
    map_join_tc(function(so,s,w,g,t,cal) {
      week.idxs = w + -week.window:week.window
      week.idxs = week.idxs[1 <= week.idxs & week.idxs <= length(w.retro.model.weeks)]
      qs = swgt.ensemble.quantiles[c(-so,-s),week.idxs,,t]
      qs = qs[!is.na(qs)]
      return(cal(swgt.ensemble.forecasts[[s,w,g,t]],qs))},
      s.outer.idx.seasons,s.idx.seasons,w.idx.weeks,g.idx.groups,t.idx.targets,
      c.calibrations)
  }

if (!file.exists(sswgtfc.forecast.values.file)) {
  saveRDS(sswgtfc.forecast.values,sswgtfc.forecast.values.file)
}

gc()

sswgtc.ensemble.forecast.evaluations.file = file.path(experiment.cache.dir,"sswgtc.ensemble.forecast.evaluations.rds")
sswgtc.ensemble.forecast.evaluations =
  if (file.exists(sswgtc.ensemble.forecast.evaluations.file)) {
    readRDS(sswgtc.ensemble.forecast.evaluations.file)
  } else {
    map_join_tc(get_evaluation,sswgtc.ensemble.forecast.values,swgtm.retro.observed.values[,,,,"Bin"],
                no_join(m.forecast.types[["Bin"]]))
  }
mode(sswgtc.ensemble.forecast.evaluations) = "numeric"

if (!file.exists(sswgtc.ensemble.forecast.evaluations.file)) {
  saveRDS(sswgtc.ensemble.forecast.evaluations,sswgtc.ensemble.forecast.evaluations.file)
}

swgtc.indexer = list(all=NULL,smear=-3:3,all=NULL,each=NULL,all=NULL)

get_calibration_weights = function(sswgtc.ensemble.forecast.evaluations, indexer.list) {
  swtc.result = array(0,dim=dim(sswgtc.ensemble.forecast.evaluations)[c(1,3,5,6)])
  for (s in 1:dim(sswgtc.ensemble.forecast.evaluations)[[1]]) {
    swtc.result[s,,,] =
      cv_apply(sswgtc.ensemble.forecast.evaluations[s,-s,,,,],
        indexer.list,
        function(train, test) {
          train = R.utils::wrap(train,list(1:4,5))
          train = pmax(train[apply(train,1,function(lst) { all(!is.na(lst))}),],-8)
          degenerate.em.weights = degenerate_em_weights(exp(train))
          return(degenerate.em.weights)
        })[,,,,,1] %>>% aperm(c(2:3,1))
  }
  dimnames(swtc.result) <- dimnames(sswgtc.ensemble.forecast.evaluations)[c(2,3,5,6)]
  return(swtc.result)
}

swtc.ensemble.weights.file = file.path(experiment.cache.dir,"swtc.ensemble.weights.rds")
swtc.ensemble.weights = get_calibration_weights(sswgtc.ensemble.forecast.evaluations,swgtc.indexer)
if (!file.exists(swtc.ensemble.weights.file)) {
  saveRDS(swtc.ensemble.weights,swtc.ensemble.weights.file)
}
gc()
swgtc.ensemble.forecast.values = apply(sswgtc.ensemble.forecast.values,3:6,diag)
swgtc.ensemble.forecast.values = unlist(swgtc.ensemble.forecast.values,recursive=FALSE)
dim(swgtc.ensemble.forecast.values) = dim(sswgtc.ensemble.forecast.values)[2:6]
dimnames(swgtc.ensemble.forecast.values) = dimnames(sswgtc.ensemble.forecast.values)[2:6]
rm(sswgtc.ensemble.forecast.values)
rm(sswgtc.ensemble.forecast.evaluations)
gc()

swgt.calibrated.ensemble.forecasts.file = file.path(experiment.cache.dir,"swgt.calibrated.ensemble.forecasts.rds")
swgt.calibrated.ensemble.forecasts =
  if (file.exists(swgt.calibrated.ensemble.forecasts.file)) {
    readRDS(swgt.calibrated.ensemble.forecasts.file)
  } else {
    map_join("*",swgtc.ensemble.forecast.values,swtc.ensemble.weights,lapply_variant=lapply) %>>%
      apply(1:4,Reduce,f="+")
  }

if (!file.exists(swgt.calibrated.ensemble.forecasts.file)) {
  saveRDS(swgt.calibrated.ensemble.forecasts,swgt.calibrated.ensemble.forecasts.file)
}

swgt.calibrated.ensemble.evaluations.file = file.path(experiment.cache.dir,"swgt.calibrated.ensemble.evaluations.rds")
swgt.calibrated.ensemble.evaluations =
  if (file.exists(swgt.calibrated.ensemble.evaluations.file)) {
    readRDS(swgt.calibrated.ensemble.evaluations.file)
  } else {
    map_join_tc(get_evaluation,swgt.calibrated.ensemble.forecasts,swgtm.retro.observed.values[,,,,"Bin"],
                no_join(m.forecast.types[["Bin"]]),lapply_variant=lapply)
  }
mode(swgt.calibrated.ensemble.evaluations) = "numeric"

if (!file.exists(swgt.calibrated.ensemble.evaluations.file)) {
  saveRDS(swgt.calibrated.ensemble.evaluations,swgt.calibrated.ensemble.evaluations.file)
}


swgt.calibrated.forecasts.file = file.path(experiment.cache.dir,"swgt.ensemble.calibrated.forecasts.rds")
swgt.calibrated.forecasts = 
  if (file.exists(swgt.calibrated.forecasts.file)) {
    readRDS(swgt.calibrated.forecasts.file)
  } else {
    map_join_tc("*",swgtfc.forecast.values,swtfc.weights,lapply_variant=lapply) %>>%
      apply(1:4,Reduce,f="+")
  }

if (!file.exists(swgt.calibrated.forecasts.file)) {
  saveRDS(swgt.calibrated.forecasts,swgt.calibrated.forecasts.file)
}

swgt.calibrated.evaluations.file = file.path(experiment.cache.dir,"swgt.ensemble.calibrated.evaluations.rds")
swgt.calibrated.evaluations =
  if (file.exists(swgt.calibrated.evaluations.file)) {
    readRDS(swgt.calibrated.evaluations.file)
  } else {
    map_join_tc(get_evaluation,swgt.calibrated.forecasts,swgtm.retro.observed.values[,,,,"Bin"],
                no_join(m.forecast.types[["Bin"]]),lapply_variant=lapply)
  }
mode(swgt.calibrated.evaluations) = "numeric"

if (!file.exists(swgt.calibrated.evaluations.file)) {
  saveRDS(swgt.calibrated.evaluations,swgt.calibrated.evaluations.file)
}
