
## Observations we use for CV evaluation:
print("CV: find observed trajectories, target multivals, target values")
gc()
sg.retro.observed.trajectories = map_join(
  get_observed_trajectory,
  s.retro.seasons, g.epigroups,
  cache.prefix=file.path(epiproject.cache.dir,"sg.retro.observed.trajectories")
)

gc()
swgt.retro.observed.multivals = map_join(
  observed_multival2,
  swg.retro.voxel.data,
  target_trajectory_preprocessor,
  t.target.specs,
  sg.retro.observed.trajectories,
  shuffle=FALSE,
  cache.prefix=file.path(epiproject.cache.dir,"swgt.retro.observed.multivals")
)

gc()
swgtm.retro.observed.values = map_join(
  observed_value2,
  swg.retro.voxel.data, t.target.specs, m.forecast.types,
  swgt.retro.observed.multivals,
  shuffle=FALSE, lapply_variant=lapply,
  cache.prefix=file.path(epiproject.cache.dir,"swgtm.retro.observed.values")
)

## Tack on additional indexer_list's for CV:
e.retro.ensemble.weighting.scheme.swgtmbf.indexer.lists =
  map_join(function(wgt.indexer.list) {
    c(retro.season.indexer, # retro study on seasons
      wgt.indexer.list, # ensemble choice for weeks, epigroups, targets
      list(each=NULL), # separately for each forecast type
      list(all=NULL, all=NULL) # always group together all backcasters, forecasters
      )
  }, e.ensemble.partial.weighting.scheme.wgt.indexer.lists,
  lapply_variant=lapply, shuffle=FALSE, show.progress=FALSE
  )

## Calculate CV ensemble weights:
old.mc.cores = getOption("mc.cores")
options("mc.cores"=min(2L,old.mc.cores))
print("CV: fit weightsets")
e.retro.ensemble.weightsets = map_join(
  function(weighting.scheme.indexer.list) {
    get_ensemble_weightset(swgtmbf.retro.component.forecast.values,
                           swgtm.retro.observed.values,
                           m.forecast.types,
                           weighting.scheme.indexer.list)
  },
  e.retro.ensemble.weighting.scheme.swgtmbf.indexer.lists,
  lapply_variant=lapply,
  cache.prefix=file.path(epiproject.cache.dir,"e.retro.ensemble.weightsets")
)
options("mc.cores"=old.mc.cores)

## Calculate CV ensemble forecasts as forecast.value's
print("CV: generate ensemble forecasts")
gc()
swgtme.retro.ensemble.forecast.values.file = file.path(file.path(epiproject.cache.dir,"swgtme.retro.ensemble.forecast.values.rds"))
swgtme.retro.ensemble.forecast.values =
  if (file.exists(swgtme.retro.ensemble.forecast.values.file)) {
    readRDS(swgtme.retro.ensemble.forecast.values.file)
  } else {
    lapply(
                 e.retro.ensemble.weightsets,
                 function(weightset) {
                   map_join(
                     ## `*`, # bad if only non-NA's are 0-weighted
                     function(forecast.value, weight) {
                       if (weight == 0) {
                         weight <- NA
                       }
                       weight * forecast.value
                     },
                     swgtmbf.retro.component.forecast.values, weightset,
                     lapply_variant=lapply, show.progress=FALSE
                   ) %>>% apply(1:5, Reduce, f=function(x,y) {
                     dplyr::coalesce(x+y, x, y)
                   })
                 }) %>>%
      simplify2arrayp() %>>%
      {names(dimnames(.))[[6L]] <- dimnamesnamesp(e.retro.ensemble.weightsets); .}
  }
if (!file.exists(swgtme.retro.ensemble.forecast.values.file)) {
  saveRDS(swgtme.retro.ensemble.forecast.values, swgtme.retro.ensemble.forecast.values.file)
}
## todo in weightset application, try removing NA values rather than 0 weights? is there a way to use row/col weighted means to do the above?
## todo speed up weightset fitting&application with Rcpp?
## todo use mclapply in map_join above instead of mclapply over ensembles...
## but broke before; need to avoid memory duplication issues

## Calculate CV ensemble forecasts as target.multicasts
gc()
swge.retro.ensemble.target.multicasts.file = file.path(epiproject.cache.dir,"swge.retro.ensemble.target.multicasts.rds")
swge.retro.ensemble.target.multicasts =
  if (file.exists(swge.retro.ensemble.target.multicasts.file)) {
    readRDS(swge.retro.ensemble.target.multicasts.file)
  } else {
    apply(swgtme.retro.ensemble.forecast.values, c(1:3,6L),
          function(tm.forecast.values) {
            list(forecast.values=tm.forecast.values)
          })
  }
if (!file.exists(swge.retro.ensemble.target.multicasts.file)) {
  saveRDS(swge.retro.ensemble.target.multicasts, swge.retro.ensemble.target.multicasts.file)
}

print('Analysis: calculate CV ensemble evaluations')

## swgtme.retro.ensemble.multibin.scores = map_join(
##     get_evaluation,
##     swgtme.retro.ensemble.forecast.values[,,,,"Bin",,drop=FALSE],
##     swgt.retro.observed.multibin.values,
##     no_join(multibin.logscore.forecast.type)
## )
## mode(swgtme.retro.ensemble.multibin.scores) <- "numeric"
## saveRDS(swgtme.retro.ensemble.multibin.scores, file.path(epiproject.cache.dir,"swgtme.retro.ensemble.multibin.scores.rds"))

swgtme.retro.ensemble.evaluations = map_join(
    get_evaluation,
    swgtme.retro.ensemble.forecast.values, swgtm.retro.observed.values, m.forecast.types
)
mode(swgtme.retro.ensemble.evaluations) <- "numeric"
saveRDS(swgtme.retro.ensemble.evaluations, file.path(epiproject.cache.dir,"swgtme.retro.ensemble.evaluations.rds"))





## apply(swgtmbf.retro.component.multibin.scores, c(7L,6L), mean, na.rm=TRUE)
## apply(swgtme.retro.ensemble.multibin.scores, c(6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[,,,,"Bin",,,drop=FALSE], c(7L,6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[,,,,"Bin",,,drop=FALSE]%>>%pmax(-10), c(7L,6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[,,,,"Bin",,,drop=FALSE]%>>%pmax(-10), c(6L,7L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[names(s.retro.seasons)[s.retro.seasons>=2010L],,,,"Bin",,,drop=FALSE]%>>%pmax(-10), c(7L,6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations, c(7L,5L), mean, na.rm=TRUE)
apply(swgtme.retro.ensemble.evaluations, 6:5, mean, na.rm=TRUE)
