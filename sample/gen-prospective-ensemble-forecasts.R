
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
                           swgtm.retro.observed.values,
                           m.forecast.types,
                           weighting.scheme.indexer.list)
  },
  e.prospective.ensemble.weighting.scheme.swgtmbf.indexer.lists
  , cache.prefix=file.path(epiproject.cache.dir,"e.prospective.ensemble.weightsets")
)

## e.prospective.ensemble.weightsets = map_join(
##   function(weighting.scheme.indexer.list) {
##     stop ("Result was not computed yet.")
##   },
##   e.prospective.ensemble.weighting.scheme.swgtmbf.indexer.lists
##   , cache.prefix=file.path(epiproject.cache.dir,"e.prospective.ensemble.weightsets")
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
  simplify2arrayp() %>>%
  {names(dimnames(.))[[6L]] <- dimnamesnamesp(e.prospective.ensemble.weightsets); .}
## todo make this work for subsets and other indexers --- get index sets from
## indexers

## Calculate CV ensemble forecasts as target.multicasts
swge.prospective.ensemble.target.multicasts =
  apply(swgtme.prospective.ensemble.forecast.values, c(1:3,6L),
        function(tm.forecast.values) {
          list(forecast.values=tm.forecast.values)
        })
