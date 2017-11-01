## author_header begin
## Copyright (C) 2017 Logan C. Brooks
##
## This file is part of epiforecast.  Algorithms included in epiforecast were developed by Logan C. Brooks, David C. Farrow, Sangwon Hyun, Shannon Gallagher, Ryan J. Tibshirani, Roni Rosenfeld, and Rob Tibshirani (Stanford University), members of the Delphi group at Carnegie Mellon University.
##
## Research reported in this publication was supported by the National Institute Of General Medical Sciences of the National Institutes of Health under Award Number U54 GM088491. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health. This material is based upon work supported by the National Science Foundation Graduate Research Fellowship Program under Grant No. DGE-1252522. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation. David C. Farrow was a predoctoral trainee supported by NIH T32 training grant T32 EB009403 as part of the HHMI-NIBIB Interfaces Initiative. Ryan J. Tibshirani was supported by NSF grant DMS-1309174.
## author_header end
## license_header begin
## epiforecast is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, version 2 of the License.
##
## epiforecast is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with epiforecast.  If not, see <http://www.gnu.org/licenses/>.
## license_header end

library("pipeR")

devtools::load_all("../epiforecast")

## Set up parallel:
options(mc.cores=parallel::detectCores()-1L)

## different location naming schemes:
fluview.location.epidata.names = c("nat", paste0("hhs",1:10))
fluview.location.spreadsheet.names = c("US National", paste0("HHS Region ",1:10))

## Load in the epidata (takes a while):
epidata.cache.dir = "~/.epiforecast-cache"
if (!dir.exists(epidata.cache.dir)) {
  dir.create(epidata.cache.dir)
}
fluview.baseline.info = fetchUpdatingResource(
  function() {
    LICENSE=RCurl::getURL("https://raw.githubusercontent.com/cdcepi/FluSight-forecasts/master/LICENSE")
    wILI_Baseline=read.csv(textConnection(RCurl::getURL("https://raw.githubusercontent.com/cdcepi/FluSight-forecasts/master/wILI_Baseline.csv")), row.names=1L, check.names=FALSE, stringsAsFactors=FALSE)
    cat("LICENSE for wILI_Baseline.csv:")
    cat(LICENSE)
    return (list(
      LICENSE=LICENSE,
      wILI_Baseline=wILI_Baseline
    ))
  },
  function(fetch.response) {
    return ()
  },
  cache.file.prefix=file.path(epidata.cache.dir,"fluview_baselines"),
  cache.invalidation.period=as.difftime(1L, units="weeks")
)
fluview.baseline.ls.mat =
  fluview.baseline.info[["wILI_Baseline"]] %>>%
  as.matrix() %>>%
  ## Adjust rownames to match 2016/2017 spreadsheet Location names:
  magrittr::set_rownames(rownames(.) %>>%
                         stringr::str_replace_all("Region", "HHS Region ") %>>%
                         dplyr::recode("National"="US National")
                         ) %>>%
  ## Re-order rows so HHS Region i is at index i and US National is at 11L:
  magrittr::extract(c(2:11,1L),) %>>%
  ## Set dimnames names:
  structure(dimnames=dimnames(.) %>>%
              setNames(c("Location", "Season"))) %>>%
  {.}
fluview.baseline.df =
  reshape2::melt(fluview.baseline.ls.mat, value.name="baseline") %>>%
  dplyr::mutate(season=Season %>>%
                  as.character() %>>%
                  stringr::str_replace_all("/.*","") %>>%
                  as.integer()) %>>%
  {.}
fluview.current.l.dfs = fluview.location.epidata.names %>>%
  setNames(fluview.location.spreadsheet.names) %>>%
  lapply(function(fluview.location.epidata.name) {
  fetchEpidataDF(
    "fluview", fluview.location.epidata.name,
    first.week.of.season=usa.flu.first.week.of.season,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_current_",fluview.location.epidata.name))
  )
})
fluview.history.l.dfs =
  fluview.location.epidata.names %>>%
  setNames(fluview.location.spreadsheet.names) %>>%
  lapply(function(fluview.location.epidata.name) {
    fetchEpidataHistoryDF(
      "fluview", fluview.location.epidata.name, 0:51,
      first.week.of.season=usa.flu.first.week.of.season,
      ## force.cache.invalidation = TRUE,
      cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_",fluview.location.epidata.name))
    )
  })

epigroup.colname = "Location"

get_voxel_data = function(season, model.week, epigroup, last.losocv.issue) {
  current.epidata.history.df = fluview.history.l.dfs[[epigroup]]
  issue = season_model.week_to_epiweek(
    season, model.week, usa.flu.first.week.of.season, 3L)
  forward.looking.history.df = mimicPastHistoryDF(
    current.epidata.history.df,
    "issue", issue,
    "epiweek", issue)
  losocv.extra.history.df = mimicPastHistoryDF(
    current.epidata.history.df,
    "issue", last.losocv.issue,
    "epiweek", last.losocv.issue) %>>%
    {.[.[["season"]] > season,]}
  epidata.history.df =
    dplyr::bind_rows(losocv.extra.history.df, forward.looking.history.df)
  forward.looking.epidata.df = mimicPastEpidataDF(
    forward.looking.history.df, issue
  )
  losocv.extra.epidata.df = mimicPastEpidataDF(
    losocv.extra.history.df, last.losocv.issue
  )
  epidata.df =
    dplyr::bind_rows(losocv.extra.epidata.df, forward.looking.epidata.df)
  baseline = fluview.baseline.df %>>%
    dplyr::filter(Location == epigroup) %>>%
    mimicPastDF("season", season, nontime.index.colnames="Location") %>>%
    magrittr::extract2("baseline")
  n.weeks.in.season = lastWeekNumber(season, 3L)
  is.inseason = usa_flu_inseason_flags(n.weeks.in.season)
  forecast.time = model_week_to_time(
    model.week, usa.flu.first.week.of.season)
  return (list(
    season = season,
    model.week = model.week,
    epigroup = epigroup,
    issue = issue,
    epidata.df = epidata.df,
    epidata.history.df = epidata.history.df,
    baseline = baseline,
    target.settings = list(
      baseline = baseline,
      is.inseason = is.inseason,
      forecast.time = forecast.time
    ),
    ## todo move to setting/task/dataset?:
    first.week.of.season = usa.flu.first.week.of.season
  ))
}

get_observation_values = function(voxel.data, target_trajectory_preprocessor, target.spec, forecast.type) {
  observation.issue = (voxel.data[["season"]]+1L)*100L + 40L
  season = voxel.data[["season"]]
  epigroup = voxel.data[["epigroup"]]
  ## todo fluview.history.l.dfs as separate argument or curried argument in this
  ## method and get_voxel_data
  ## epidata.df = mimicPastEpidataDF(
  ##   fluview.history.l.dfs[[epigroup]], observation.issue)
  ## xxx check that this doesn't change anything:
  epidata.df = fluview.current.l.dfs[[epigroup]]
  observed.trajectory = epidata.df %>>%
    {.[["wili"]][.[["season"]]==season]}
  observation.as.target.forecast = target_forecast2(
    voxel.data, target_trajectory_preprocessor, target.spec,
    match.new.dat.sim(observed.trajectory))
  observation.as.target.forecast[["method.settings"]] <-
    c(observation.as.target.forecast[["method.settings"]],
      uniform.pseudoweight.total=0,
      smooth.sim.targets=FALSE)
  observion.values =
    forecast_value2(voxel.data, target.spec, forecast.type, observation.as.target.forecast)
  return (observion.values)
}

backfill_ignorant_backsim = function(voxel.data, signal.name) {
  old.dat = voxel.data[["epidata.df"]] %>>%
    dplyr::filter(season != voxel.data[["season"]]) %>>%
    dplyr::group_by(season) %>>%
    ## todo remove hard reference to wili
    dplyr::filter(!any(is.na(wili))) %>>%
    dplyr::ungroup() %>>%
    {split(.[[signal.name]], .[["Season"]])}
  new.dat = voxel.data[["epidata.df"]] %>>%
    dplyr::filter(season == voxel.data[["season"]]) %>>%
    magrittr::extract2(signal.name)
  new.dat.sim = match.new.dat.sim(new.dat)
  voxel.Season = season_to_Season(
    voxel.data[["season"]], voxel.data[["first.week.of.season"]])
  ## concatenate new.dat.sim onto old.dat, setting the :
  full.dat = c(old.dat,
               setNames(list(new.dat.sim), voxel.Season))
  return (full.dat)
}

s.retro.seasons = c(2008:2011) %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
w.retro.model.weeks = (40:45) %>>%
  stats::setNames(paste0("MW",.)) %>>%
  with_dimnamesnames("Model Week")
g.epigroups = fluview.location.spreadsheet.names %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Location")
last.losocv.issue = 201039L
b.backcasters = list(
  ignorant=backfill_ignorant_backsim
) %>>%
  with_dimnamesnames("Backcaster")
signal.name = "wili"
f.forecasters = list(
  "Delphi_DeltaDensity_PackageDefaults"=twkde.sim,
  ## "Delphi_DeltaDensity_200Curves"=function(...) twkde.sim(..., max.n.sims=200L),
  "Basis Regression"=br.sim,
  "Empirical Bayes"=eb.sim,
  ## "Empirical"=empirical.trajectories.sim,
  "Empirical Trajectories"=empirical.trajectories.sim,
  "Empirical Futures"=empirical.futures.sim,
  "Uniform"=uniform_forecast
) %>>%
  with_dimnamesnames("Forecaster")
target_trajectory_preprocessor = flusight2016_target_trajectory_preprocessor
t.target.specs = flusight2016.target.specs %>>%
  with_dimnamesnames("Target")
f.forecast.types = flusight2016.proxy.forecast.types %>>%
  with_dimnamesnames("Type")

## CV input data
swg.retro.voxel.data = map_join(
  get_voxel_data,
  s.retro.seasons, w.retro.model.weeks, g.epigroups,
  last.losocv.issue,
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/swg.retro.voxel.data/swg.retro.voxel.data"
)
## xxx this (retro voxel data) is going to consume a bunch of memory...

## CV backcasts
swgb.retro.full.dats = map_join(
  get_backcast,
  swg.retro.voxel.data, signal.name, b.backcasters,
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/swgb.retro.full.dats/swgb.retro.full.dats"
)

## CV forecasts as target_multicast objects:
swgbf.retro.component.target.multicasts = map_join(
  target_multicast,
  swg.retro.voxel.data, swgb.retro.full.dats, f.forecasters,
  target_trajectory_preprocessor,
  no_join(t.target.specs),
  no_join(f.forecast.types),
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/swgbf.retro.component.target.multicasts/swgbf.retro.component.target.multicasts"
)

## CV forecasts as forecast.value's:
swgtmbf.retro.component.forecast.values =
  swgbf.retro.component.target.multicasts %>>%
  ## first, get forecast.value's in swgbf.tm format:
  map_join(f=`[[`, "forecast.values") %>>%
  ## un-nest lists to get swgbftm format:
  map_join(f=`[[`,
           named_arrayvec_to_name_arrayvec(t.target.specs),
           named_arrayvec_to_name_arrayvec(f.forecast.types)
           ) %>>%
  ## permute dimension order to get desired swgtmbf format:
  aperm(c(1:3,6:7,4:5))

## Observations we use for CV evaluation:
swgtm.retro.observation.values = map_join(
  get_observation_values,
  swg.retro.voxel.data,
  target_trajectory_preprocessor, t.target.specs, f.forecast.types,
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/swgtm.retro.observation.values/swgtm.retro.observation.values"
)

## Specify portions of cv_apply indexer lists corresponding to model week,
## epigroup, target:
e.ensemble.partial.weighting.scheme.wgt.indexer.lists = list(
  "constant-weights" = list(all=NULL, all=NULL, all=NULL),
  "target-based" = list(all=NULL, all=NULL, each=NULL),
  "target-3time-based" = list(smear=-1:1, all=NULL, each=NULL),
  "target-9time-based" = list(smear=-4:4, all=NULL, each=NULL)
) %>>% with_dimnamesnames("Ensemble weighting scheme")

## Use LOSOCV on seasons < 2010L (train on other seasons < 2010L), and
## oneahead on seasons >= 2010L (train on seasons < test season)
retro.season.indexer = list(loo_oneahead=match(2010L,s.retro.seasons))

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

## e.retro.ensemble.weightsets = lapply(
##   retro.ensemble.weighting.scheme.swgtmbf.indexer.lists,
##   function(weighting.scheme.indexer.list) {
##     get_ensemble_weightset(swgtmbf.retro.component.forecast.values,
##                            swgtm.retro.observation.values,
##                            f.forecast.types,
##                            weighting.scheme.indexer.list)
##   }
## ) %>>% with_dimnamesnames(dimnamesnamesp(ensemble.partial.weighting.scheme.wgt.indexer.lists))
e.retro.ensemble.weightsets = map_join(
  function(weighting.scheme.indexer.list) {
    get_ensemble_weightset(swgtmbf.retro.component.forecast.values,
                           swgtm.retro.observation.values,
                           f.forecast.types,
                           weighting.scheme.indexer.list)
  },
  e.retro.ensemble.weighting.scheme.swgtmbf.indexer.lists,
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/e.retro.ensemble.weightsets/e.retro.ensemble.weightsets"
)

## Calculate CV ensemble forecasts as forecast.value's
swgtme.retro.ensemble.forecast.values = pbmcapply::pbmclapply(
  e.retro.ensemble.weightsets,
  function(weightset) {
    map_join(
      `*`,
      swgtmbf.retro.component.forecast.values, weightset
    ) %>>% apply(1:5, Reduce, f=function(x,y) {
      dplyr::coalesce(x+y, x, y)
    })
  }) %>>%
  simplify2array() %>>%
  {names(dimnames(.))[[6L]] <- dimnamesnamesp(e.retro.ensemble.weightsets); .}

## Calculate CV ensemble forecasts as target.multicasts
swge.retro.ensemble.target.multicasts.file = "~/files/nosync/forecasts/flusight-tiny-run/swge.retro.ensemble.target.multicasts"
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

## target_multicast_percent_plot(swge.retro.ensemble.target.multicasts[[1L]],
##                               swg.retro.voxel.data[[1L]],
##                               t.target.specs,
##                               f.forecast.types)

swgtmbf.retro.component.evaluations = map_join(
  get_evaluation,
  swgtmbf.retro.component.forecast.values, swgtm.retro.observation.values, f.forecast.types
) %>>%
  {mode(.) <- "numeric"; .}

swgtme.retro.ensemble.evaluations = map_join(
  get_evaluation,
  swgtme.retro.ensemble.forecast.values, swgtm.retro.observation.values, f.forecast.types
) %>>%
  {mode(.) <- "numeric"; .}

## apply(swgtmbf.retro.component.evaluations, 5:7, mean, na.rm=TRUE)
apply(swgtmbf.retro.component.evaluations, c(5L,7L), mean, na.rm=TRUE)
apply(swgtme.retro.ensemble.evaluations, 5:6, mean, na.rm=TRUE)

## swbf.retro.forecast.spreadsheets = map_join(
##   target_multicast_epigroup_forecast_table,
##   swgbf.retro.component.target.multicasts,
##   swg.retro.voxel.data,
##   no_join(t.target.specs), no_join(f.forecast.types)
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

current.issue.sw =
  fluview.current.l.dfs[[1L]] %>>%
  dplyr::filter(season == max(season)) %>>%
  {.[!is.na(.[[signal.name]]),]} %>>%
  dplyr::filter(model.week == max(model.week)) %>>%
  dplyr::select(season, model.week)

s.prospective.seasons = current.issue.sw[["season"]] %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
w.prospective.model.weeks = current.issue.sw[["model.week"]] %>>%
  stats::setNames(paste0("MW",.)) %>>%
  with_dimnamesnames("Model Week")

swg.prospective.voxel.data = map_join(
  get_voxel_data,
  s.prospective.seasons, w.prospective.model.weeks, g.epigroups,
  last.losocv.issue)

swgb.prospective.full.dats = map_join(
  get_backcast,
  swg.prospective.voxel.data, signal.name, b.backcasters,
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/swgb.prospective.full.dats/swgb.prospective.full.dats"
)

swgbf.prospective.component.target.multicasts = map_join(
  target_multicast,
  swg.prospective.voxel.data, swgb.prospective.full.dats, f.forecasters,
  target_trajectory_preprocessor,
  no_join(t.target.specs),
  no_join(f.forecast.types),
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/swgbf.prospective.component.target.multicasts/swgbf.prospective.component.target.multicasts"
)

swgtmbf.prospective.component.forecast.values =
  swgbf.prospective.component.target.multicasts %>>%
  ## first, get forecast.value's in swgbf.tm format:
  map_join(f=`[[`, "forecast.values") %>>%
  ## un-nest lists to get swgbftm format:
  map_join(f=`[[`,
           named_arrayvec_to_name_arrayvec(t.target.specs),
           named_arrayvec_to_name_arrayvec(f.forecast.types)
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

e.prospective.ensemble.weightsets = map_join(
  function(weighting.scheme.indexer.list) {
    get_ensemble_weightset(swgtmbf.retro.component.forecast.values,
                           swgtm.retro.observation.values,
                           f.forecast.types,
                           weighting.scheme.indexer.list)
  },
  e.prospective.ensemble.weighting.scheme.swgtmbf.indexer.lists,
  cache.prefix="~/files/nosync/forecasts/flusight-tiny-run/e.prospective.ensemble.weightsets/e.prospective.ensemble.weightsets"
)

swgtme.prospective.ensemble.forecast.values = lapply(
  e.prospective.ensemble.weightsets,
  function(weightset) {
    map_join(
      `*`,
      swgtmbf.prospective.component.forecast.values, weightset,
      eltname.mismatch.behavior="intersect"
    ) %>>% apply(1:5, Reduce, f=`+`)
  }) %>>%
  simplify2array() %>>%
  {names(dimnames(.))[[6L]] <- dimnamesnamesp(e.prospective.ensemble.weightsets); .}
## todo make this work for subsets and other indexers --- get index sets from
## indexers

## Calculate CV ensemble forecasts as target.multicasts
swge.prospective.ensemble.target.multicasts.file = "~/files/nosync/forecasts/flusight-tiny-run/swge.prospective.ensemble.target.multicasts"
swge.prospective.ensemble.target.multicasts =
  if (file.exists(swge.prospective.ensemble.target.multicasts.file)) {
    readRDS(swge.prospective.ensemble.target.multicasts.file)
  } else {
    apply(swgtme.prospective.ensemble.forecast.values, c(1:3,6L),
          function(tm.forecast.values) {
          list(forecast.values=tm.forecast.values)
          })
  }
if (!file.exists(swge.prospective.ensemble.target.multicasts.file)) {
  saveRDS(swge.prospective.ensemble.target.multicasts, swge.prospective.ensemble.target.multicasts.file)
}

## target_multicast_percent_plot(swge.prospective.ensemble.target.multicasts[[1L]],
##                               swg.prospective.voxel.data[[1L]],
##                               t.target.specs,
##                               f.forecast.types)

## target_multicast_percent_plot(swgbf.prospective.component.target.multicasts[[1L,1L,1L,1L,"Uniform"]],
##                               swg.prospective.voxel.data[[1L]],
##                               t.target.specs,
##                               f.forecast.types)

## target_multicast_percent_plot(swgbf.prospective.component.target.multicasts[[1L,1L,1L,1L,"Empirical Trajectories"]],
##                               swg.prospective.voxel.data[[1L]],
##                               t.target.specs,
##                               f.forecast.types)

## target_multicast_percent_plot(swgbf.prospective.component.target.multicasts[[1L,1L,1L,1L,"Empirical Futures"]],
##                               swg.prospective.voxel.data[[1L]],
##                               t.target.specs,
##                               f.forecast.types)

## target_multicast_percent_plot(swgbf.prospective.component.target.multicasts[[1L,1L,1L,1L,"Empirical Bayes"]],
##                               swg.prospective.voxel.data[[1L]],
##                               t.target.specs,
##                               f.forecast.types)

target_multicast_percent_plot(swge.prospective.ensemble.target.multicasts[[1L,1L,1L,"constant-weights"]],
                              swg.prospective.voxel.data[[1L]],
                              t.target.specs,
                              f.forecast.types)

target_multicast_percent_plot(swge.prospective.ensemble.target.multicasts[[1L,1L,1L,"target-based"]],
                              swg.prospective.voxel.data[[1L]],
                              t.target.specs,
                              f.forecast.types)

target_multicast_percent_plot(swge.prospective.ensemble.target.multicasts[[1L,1L,1L,"target-3time-based"]],
                              swg.prospective.voxel.data[[1L]],
                              t.target.specs,
                              f.forecast.types)

target_multicast_percent_plot(swge.prospective.ensemble.target.multicasts[[1L,1L,1L,"target-9time-based"]],
                              swg.prospective.voxel.data[[1L]],
                              t.target.specs,
                              f.forecast.types)


## swe.retro.forecast.spreadsheets = map_join(
##   target_multicast_epigroup_forecast_table,
##   swge.prospective.ensemble.target.multicasts,
##   swg.prospective.voxel.data,
##   no_join(t.target.specs), no_join(f.forecast.types)
## ) %>>%
##   apply(c(1:2,4:5), dplyr::bind_rows)

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
