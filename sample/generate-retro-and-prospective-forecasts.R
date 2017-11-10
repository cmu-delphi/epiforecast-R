## CV input data
print("CV: select available input data")
gc()
swg.retro.voxel.data =
  tryCatch({
    map_join(
      get_voxel_data,
      s.retro.seasons, w.retro.model.weeks, g.epigroups,
      last.losocv.issue,
      cache.prefix=file.path(epiproject.cache.dir,"swg.retro.voxel.data/swg.retro.voxel.data"),
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
      cache.prefix=file.path(epiproject.cache.dir,"swg.retro.voxel.data/swg.retro.voxel.data"),
      use.proxy=TRUE
    )
  })
## swg.retro.voxel.data = map_join(
##   get_voxel_data,
##   s.retro.seasons, w.retro.model.weeks, g.epigroups,
##   last.losocv.issue,
##   lapply_variant=pbmclapply_no_preschedule,
##   cache.prefix=file.path(epiproject.cache.dir,"swg.retro.voxel.data/swg.retro.voxel.data")
## )
## xxx this (retro voxel data) is going to consume a bunch of memory...

## CV backcasts
print("CV: generate backcasts")
swgb.retro.full.dats = map_join(
  get_backcast,
  swg.retro.voxel.data, signal.name, b.backcasters,
  cache.prefix=file.path(epiproject.cache.dir,"swgb.retro.full.dats/swgb.retro.full.dats")
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
  ## lapply_variant=pbmcapply::pbmclapply,
  ## lapply_variant=pbmclapply_no_preschedule,
  ## lapply_variant = lapply,
  cache.prefix=file.path(epiproject.cache.dir,"swgbf.retro.component.target.multicasts/swgbf.retro.component.target.multicasts"),
  use.proxy=TRUE
)
## xxx loading from many cache files is slow; reduce # of cache files?

## CV forecasts as forecast.value's:
gc()
swgtmbf.retro.component.forecast.values =
  map_join(swgbf.retro.component.target.multicasts,
           f=`[[`, "forecast.values",
           lapply_variant=lapply) %>>%
  simplify2array() %>>%
  {
    original = .
    dim(.) <- c(dim(original)[1:2], dim(swgbf.retro.component.target.multicasts))
    dimnames(.) <- c(dimnames(original)[1:2], dimnames(swgbf.retro.component.target.multicasts))
    .
  } %>>%
  aperm(c(3:5,1:2,6:7))

## Observations we use for CV evaluation:
print("CV: find observation values")
gc()
swgtm.retro.observation.values = map_join(
  get_observation_values,
  swg.retro.voxel.data,
  target_trajectory_preprocessor, t.target.specs, m.forecast.types,
  cache.prefix=file.path(epiproject.cache.dir,"swgtm.retro.observation.values/swgtm.retro.observation.values")
)

## Tack on additional indexer_list's for CV:
e.retro.ensemble.weighting.scheme.swgtmbf.indexer.lists =
  map_join(function(wgt.indexer.list) {
    c(retro.season.indexer, # retro study on seasons
      wgt.indexer.list, # ensemble choice for weeks, epigroups, targets
      list(each=NULL), # separately for each forecast type
      list(all=NULL, all=NULL) # always group together all backcasters, forecasters
      )
  }, e.ensemble.partial.weighting.scheme.wgt.indexer.lists)

## Calculate CV ensemble weights:
print("CV: fit weightsets")
e.retro.ensemble.weightsets = map_join(
  function(weighting.scheme.indexer.list) {
    get_ensemble_weightset(swgtmbf.retro.component.forecast.values,
                           swgtm.retro.observation.values,
                           m.forecast.types,
                           weighting.scheme.indexer.list)
  },
  e.retro.ensemble.weighting.scheme.swgtmbf.indexer.lists,
  lapply_variant=lapply,
  cache.prefix=file.path(epiproject.cache.dir,"e.retro.ensemble.weightsets/e.retro.ensemble.weightsets")
)

## Calculate CV ensemble forecasts as forecast.value's
print("CV: generate ensemble forecasts")
gc()
swgtme.retro.ensemble.forecast.values.file = file.path(file.path(epiproject.cache.dir,"swgtme.retro.ensemble.forecast.values.rds"))
swgtme.retro.ensemble.forecast.values =
  if (file.exists(swgtme.retro.ensemble.forecast.values.file)) {
    readRDS(swgtme.retro.ensemble.forecast.values.file)
  } else {
    parallel::mclapply(
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
      simplify2array() %>>%
      {names(dimnames(.))[[6L]] <- dimnamesnamesp(e.retro.ensemble.weightsets); .}
  }
if (!file.exists(swgtme.retro.ensemble.forecast.values.file)) {
  saveRDS(swgtme.retro.ensemble.forecast.values, swgtme.retro.ensemble.forecast.values.file)
}
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

print("Analysis: calculate CV evaluations")
gc()
swgtmbf.retro.component.evaluations = map_join(
  get_evaluation,
  swgtmbf.retro.component.forecast.values, swgtm.retro.observation.values, m.forecast.types
)
mode(swgtmbf.retro.component.evaluations) <- "numeric"
## fixme sometimes this results in errors due to NULL's appearing in the
## evaluations, but re-running the evaluations seems to work... memory issues? gc beforehand?

swgtme.retro.ensemble.evaluations = map_join(
  get_evaluation,
  swgtme.retro.ensemble.forecast.values, swgtm.retro.observation.values, m.forecast.types
)
mode(swgtme.retro.ensemble.evaluations) <- "numeric"

## ## apply(swgtmbf.retro.component.evaluations, 5:7, mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations, c(5L,7L), mean, na.rm=TRUE)
apply(swgtme.retro.ensemble.evaluations, 5:6, mean, na.rm=TRUE)

## apply(swgtmbf.retro.component.evaluations[(2003:2009)%>>%paste0("/",.+1L),,,,,,,drop=FALSE],
##       c(5L,7L), mean, na.rm=TRUE)
## apply(swgtme.retro.ensemble.evaluations[(2003:2009)%>>%paste0("/",.+1L),,,,,,drop=FALSE],
##       5:6, mean, na.rm=TRUE)

## apply(swgtmbf.retro.component.evaluations, c(7L,4:5), mean, na.rm=TRUE)
## apply(swgtme.retro.ensemble.evaluations, c(6L,4:5), mean, na.rm=TRUE)

## apply(swgtme.retro.ensemble.evaluations, c(1,4:5), mean, na.rm=TRUE)

## todo another level of stacking?  ensemble of ensembles?

## swbf.retro.forecast.spreadsheets = map_join(
##   target_multicast_epigroup_forecast_table,
##   swgbf.retro.component.target.multicasts,
##   swg.retro.voxel.data,
##   no_join(t.target.specs), no_join(m.forecast.types)
## ) %>>%
##   apply(c(1:2,4:5), dplyr::bind_rows)

## spreadsheet.template = readr::read_csv("~/long_flu_submission_template_1718.csv", col_types=readr::cols())
## attr(spreadsheet.template, "spec") <- NULL
## spreadsheet.to.check = swbf.retro.forecast.spreadsheets[["2010/2011",1L,1L,1L]]

## spreadsheet.to.check %>>% dplyr::mutate(A="a") %>>% dplyr::full_join(qwer %>>% dplyr::mutate(B="b"), by=c("Location","Target","Type","Unit","Bin_start_incl","Bin_end_notincl")) %>>% dplyr::filter(is.na(A) | is.na(B))
## class(spreadsheet.to.check[["Value"]])==class(spreadsheet.template[["Value"]])
## identical(sapply(spreadsheet.to.check, class), sapply(spreadsheet.template, class))
## identical(sapply(spreadsheet.to.check, typeof), sapply(spreadsheet.template, typeof))
## identical(sapply(spreadsheet.to.check, mode), sapply(spreadsheet.template, mode))
## identical(sapply(spreadsheet.to.check, storage.mode), sapply(spreadsheet.template, storage.mode))
## Map(all.equal,
##     spreadsheet.to.check %>>% dplyr::select(-Value),
##     spreadsheet.template %>>% dplyr::select(-Value))
## Map(identical,
##     spreadsheet.to.check %>>% dplyr::select(-Value),
##     spreadsheet.template %>>% dplyr::select(-Value))
## identical(class(spreadsheet.to.check), class(spreadsheet.template))
## identical(typeof(spreadsheet.to.check), typeof(spreadsheet.template))
## identical(mode(spreadsheet.to.check), mode(spreadsheet.template))
## identical(storage.mode(spreadsheet.to.check), storage.mode(spreadsheet.template))
## ## all.equal overload does not check attributes; if "spec" attribute set by
## ## readr is included, the all.equal check passes but the identical check fails;
## ## removing this attribute makes the identical check pass
## all.equal(spreadsheet.to.check %>>% dplyr::select(-Value),
##           spreadsheet.template %>>% dplyr::select(-Value))
## identical(spreadsheet.to.check %>>% dplyr::select(-Value),
##           spreadsheet.template %>>% dplyr::select(-Value))

gc()
s.prospective.seasons = current.issue.sw[["season"]] %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
w.prospective.model.weeks = current.issue.sw[["model.week"]] %>>%
  stats::setNames(paste0("MW",.)) %>>%
  with_dimnamesnames("Model Week")

print("Current season: select available data")
swg.prospective.voxel.data = map_join(
  get_voxel_data,
  s.prospective.seasons, w.prospective.model.weeks, g.epigroups,
  last.losocv.issue)

print("Current season: generate backcasts")
swgb.prospective.full.dats = map_join(
  get_backcast,
  swg.prospective.voxel.data, signal.name, b.backcasters,
  cache.prefix=file.path(epiproject.cache.dir,"swgb.prospective.full.dats/swgb.prospective.full.dats")
)

print("Current season: generate component forecasts")
swgbf.prospective.component.target.multicasts = map_join(
  target_multicast,
  swg.prospective.voxel.data, swgb.prospective.full.dats, f.forecasters,
  target_trajectory_preprocessor,
  no_join(t.target.specs),
  no_join(m.forecast.types),
  cache.prefix=file.path(epiproject.cache.dir,"swgbf.prospective.component.target.multicasts/swgbf.prospective.component.target.multicasts")
)

swgtmbf.prospective.component.forecast.values =
  swgbf.prospective.component.target.multicasts %>>%
  ## first, get forecast.value's in swgbf.tm format:
  map_join(f=`[[`, "forecast.values") %>>%
  ## un-nest lists to get swgbftm format:
  map_join(f=`[[`,
           named_arrayvec_to_name_arrayvec(t.target.specs),
           named_arrayvec_to_name_arrayvec(m.forecast.types)
           ) %>>%
  ## permute dimension order to get desired swgtmbf format:
  aperm(c(1:3,6:7,4:5))

## Tack on additional indexer_list's for prospective forecasting:
e.prospective.ensemble.weighting.scheme.swgtmbf.indexer.lists =
  map_join(function(wgt.indexer.list) {
    c(list(all=NULL), # use all past seasons
      wgt.indexer.list, # ensemble choice for weeks, epigroups, targets
      list(each=NULL), # separately for each forecast type
      list(all=NULL, all=NULL) # always group together all backcasters, forecasters
      )
  }, e.ensemble.partial.weighting.scheme.wgt.indexer.lists)

print("Current season: fit ensemble weightsets")
e.prospective.ensemble.weightsets = map_join(
  function(weighting.scheme.indexer.list) {
    get_ensemble_weightset(swgtmbf.retro.component.forecast.values,
                           swgtm.retro.observation.values,
                           m.forecast.types,
                           weighting.scheme.indexer.list)
  },
  e.prospective.ensemble.weighting.scheme.swgtmbf.indexer.lists,
  cache.prefix=file.path(epiproject.cache.dir,"e.prospective.ensemble.weightsets/e.prospective.ensemble.weightsets")
)

## e.prospective.ensemble.weightsets = map_join(
##   function(weighting.scheme.indexer.list) {
##     stop ("Result was not computed yet.")
##   },
##   e.prospective.ensemble.weighting.scheme.swgtmbf.indexer.lists,
##   cache.prefix=file.path(epiproject.cache.dir,"e.prospective.ensemble.weightsets/e.prospective.ensemble.weightsets")
## )

print("Current season: generate ensemble forecasts")
swgtme.prospective.ensemble.forecast.values = lapply(
  e.prospective.ensemble.weightsets,
  function(weightset) {
    map_join(
      ## `*`, # bad if only non-NA's are 0-weighted
      function(forecast.value, weight) {
        if (weight == 0) {
          weight <- NA
        }
        weight * forecast.value
      },
      swgtmbf.prospective.component.forecast.values, weightset,
      eltname.mismatch.behavior="intersect"
    ) %>>% apply(1:5, Reduce, f=function(x,y) {
      dplyr::coalesce(x+y, x, y)
    })
  }) %>>%
  simplify2array() %>>%
  {names(dimnames(.))[[6L]] <- dimnamesnamesp(e.prospective.ensemble.weightsets); .}
## todo make this work for subsets and other indexers --- get index sets from
## indexers

## Calculate CV ensemble forecasts as target.multicasts
swge.prospective.ensemble.target.multicasts =
  apply(swgtme.prospective.ensemble.forecast.values, c(1:3,6L),
        function(tm.forecast.values) {
          list(forecast.values=tm.forecast.values)
        })

## todo automatic appropriate parallel_dim_i setting
## fixme finish documenting map_join
## todo use simplified sim's to strip out unnecessary elements
## todo sw.epidata --- if will be borrowing across epigroups...
## todo add prospective forecasts...
## todo allow for subsets in ensemble... need to use cv_apply to get the forecasts again
## xxx target.settings -> voxel.settings?
## todo also need to strip out information from target forecasts if remains after using simplified sim's
## todo dots to lists in target spec functions, maybe elsewhere as well
## todo could also have save season-evaluation epidata to hopefully speed up calculating observed values
## todo reshape so that can include labels without exploding in size
## todo don't store sim's
## reverse order on the indexing?
## fixme no uniform fallback for Point; maybe make fallback the empirical distribution since it is smoothed now and has uniform pseudoweight?  or just remove since everything has uniform pseudoweight?
## fixme no Bin smoothing by default? warn about plotting target forecasts?
## todo joint epidata? voxel settings?
## joint backcast? joint forecast?
## fixme Virgin Islands --- little history
## todo document no_join
## fixme check on Season onset Point predictions of "none" vs. NA vs. ...
## fixme EB weight sum is too large
## fixme smooth sim targets with a Laplace kernel? or a spike + slab type --- can inflate bw to make up for mass on spike?
## fixme should adjust dimension ordering... for col major feel
## todo deal with parallel-related memory issues --- duplicate swg.voxel.data...
## todo work on speed getting voxel data --- possible to avoid storage?
## xxx consider just basing everything on filesystem contracts... no need to hold everything in memory
## fixme try to solve memory issues with mclapply env's? require interaction with disk?
## todo weighted cv_apply smearing schemes (boxcar kernel -> other kernels)
## todo tiny subset run to do some testing and development on
## fixme table verification
## fixme empirical distribution is using target's baseline with empirical curves again instead of pairing baselines and curves
## todo better fill-in in 2010 for states
## todo remove Season onset target for states
## fixme modularize into functions, document, ...
