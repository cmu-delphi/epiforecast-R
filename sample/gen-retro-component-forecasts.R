
## CV input data
print("CV: select available input data")
gc()
swg.retro.voxel.data =
  tryCatch({
    map_join(
      get_voxel_data,
      s.retro.seasons, w.retro.model.weeks, g.epigroups,
      last.losocv.issue,
      cache.prefix=file.path(epiproject.cache.dir,"swg.retro.voxel.data"),
      use.proxy=TRUE
    )
  },
  ## issues with parallel package returning long vector results from large runs...
  error=function(e) {
    print ("Encountered error preparing voxel data in parallel.  Attempting to read cache files sequentially with no progress bar --- this make take a while.")
    map_join(
      get_voxel_data,
      s.retro.seasons, w.retro.model.weeks, g.epigroups,
      last.losocv.issue,
      lapply_variant=lapply,
      cache.prefix=file.path(epiproject.cache.dir,"swg.retro.voxel.data"),
      use.proxy=TRUE
    )
  })

## Version of the above grouped by season and model week (an array of arrays of objects):
sw.g.retro.voxel.data = map_join(
  function(swg.array, s,w) {
    swg.array[s,w,,drop=FALSE] %>>%
      select_dims(1:2, "drop") # drop s & w dimensions, but keep g dimension even if size 1
  },
  no_join(swg.retro.voxel.data),
  named_arrayvec_to_name_arrayvec(s.retro.seasons),
  named_arrayvec_to_name_arrayvec(w.retro.model.weeks),
  lapply_variant=lapply, show.progress=FALSE
)

## CV backcasts
print("CV: generate backcasts")
swgb.retro.full.dats = map_join(
  get_backcast,
  swg.retro.voxel.data, sw.g.retro.voxel.data, source.name, signal.name, b.backcasters,
  epidata_df_to_chopped_trajectory_df=epidata_df_to_chopped_trajectory_df,
  use.proxy=TRUE,
  cache.prefix=file.path(epiproject.cache.dir,"swgb.retro.full.dats")
)

## CV forecasts as target_multicast objects:
print("CV: generate component forecasts")
gc()
swgbf.retro.component.target.multicasts = map_join(
  target_multicast,
  swg.retro.voxel.data, swgb.retro.full.dats, f.forecasters,
  target_trajectory_preprocessor,
  no_join(t.target.specs),
  no_join(m.forecast.types),
  full_dat_fixup=full_dat_fixup,
  ## lapply_variant = lapply,
  cache.prefix=file.path(epiproject.cache.dir,"swgbf.retro.component.target.multicasts"),
  shuffle=FALSE,
  use.proxy=TRUE
)
## xxx loading from many cache files is slow; reduce # of cache files?

## CV forecasts as forecast.value's:
gc()
swgtmbf.retro.component.forecast.values =
  map_join(swgbf.retro.component.target.multicasts,
           f=`[[`, "forecast.values",
           shuffle=FALSE#, lapply_variant=lapply
           ) %>>%
  simplify2arrayp() %>>%
  {
    original = .
    dim(.) <- c(dim(original)[1:2], dim(swgbf.retro.component.target.multicasts))
    dimnames(.) <- c(dimnames(original)[1:2], dimnames(swgbf.retro.component.target.multicasts))
    rm(original) # somehow this reaches the global environment and takes up memory
    gc()
    .
  } %>>%
  aperm(c(3:5,1:2,6:7))

print('Analysis: calculate CV component evaluations')
old.mc.cores = getOption("mc.cores")
options("mc.cores"=min(2L,old.mc.cores))
gc()
swgtmbf.retro.component.evaluations = map_join(
  get_evaluation,
  swgtmbf.retro.component.forecast.values, swgtm.retro.observed.values, m.forecast.types
)
mode(swgtmbf.retro.component.evaluations) <- "numeric"
saveRDS(swgtmbf.retro.component.evaluations, file.path(epiproject.cache.dir,"swgtmbf.retro.component.evaluations.rds"))
options("mc.cores"=old.mc.cores)
## fixme sometimes this results in errors due to NULL's appearing in the
## evaluations, but re-running the evaluations seems to work... memory issues? gc beforehand?


swgt.retro.observed.multibin.values = map_join(
    observed_value2,
    swg.retro.voxel.data, t.target.specs, no_join(multibin.logscore.forecast.type),
    swgt.retro.observed.multivals
)

## swgtmbf.retro.component.multibin.scores = map_join(
##     get_evaluation,
##     swgtmbf.retro.component.forecast.values[,,,,"Bin",,,drop=FALSE],
##     swgt.retro.observed.multibin.values,
##     no_join(multibin.logscore.forecast.type)
## )
## mode(swgtmbf.retro.component.multibin.scores) <- "numeric"

## saveRDS(swgtmbf.retro.component.multibin.scores, file.path(epiproject.cache.dir,"swgtmbf.retro.component.multibin.scores.rds"))





## apply(swgtmbf.retro.component.multibin.scores, c(7L,6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[,,,,"Bin",,,drop=FALSE], c(7L,6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[,,,,"Bin",,,drop=FALSE]%>>%pmax(-10), c(7L,6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[,,,,"Bin",,,drop=FALSE]%>>%pmax(-10), c(6L,7L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations[names(s.retro.seasons)[s.retro.seasons>=2010L],,,,"Bin",,,drop=FALSE]%>>%pmax(-10), c(7L,6L), mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations, c(7L,5L), mean, na.rm=TRUE)
