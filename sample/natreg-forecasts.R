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
## options(mc.cores=parallel::detectCores()-3L)

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
g.fluview.current.dfs = fluview.location.epidata.names %>>%
  setNames(fluview.location.spreadsheet.names) %>>%
  lapply(function(fluview.location.epidata.name) {
  fetchEpidataDF(
    "fluview", fluview.location.epidata.name,
    first.week.of.season=usa.flu.first.week.of.season,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_current_",fluview.location.epidata.name))
  )
})
g.fluview.history.dfs =
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
nowcast.lead = 1L
g.nowcast.current.dfs = fluview.location.epidata.names %>>%
  setNames(fluview.location.spreadsheet.names) %>>%
  lapply(function(fluview.location.epidata.name) {
  fetchEpidataDF(
    "nowcast", fluview.location.epidata.name,
    first.week.of.season=usa.flu.first.week.of.season,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("nowcast_",fluview.location.epidata.name))
  ) %>>%
    dplyr::mutate(issue = add_epiweek_integer(epiweek, -nowcast.lead)) %>>%
    dplyr::mutate(lag = subtract_epiweek_epiweek(issue, epiweek))
})

epigroup.colname = "Location"

get_voxel_data = function(season, model.week, epigroup, last.losocv.issue) {
  current.epidata.history.df = g.fluview.history.dfs[[epigroup]][
    c("season","model.week","epiweek","issue","lag","wili")
  ]
  current.nowcast.df = g.nowcast.current.dfs[[epigroup]]
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
  forward.looking.nowcast.df = mimicPastDF(
    current.nowcast.df,
    "issue", issue,
    "epiweek", add_epiweek_integer(issue, nowcast.lead)
  ) %>>%
    augmentWeeklyDF()
  losocv.extra.nowcast.df = mimicPastDF(
    current.nowcast.df,
    "issue", last.losocv.issue,
    "epiweek", add_epiweek_integer(last.losocv.issue, nowcast.lead)
  ) %>>%
    {.[.[["season"]] > season,]} %>>%
    augmentWeeklyDF()
  nowcast.df =
    dplyr::bind_rows(losocv.extra.nowcast.df, forward.looking.nowcast.df)
  baseline = fluview.baseline.df %>>%
    dplyr::filter(Location == epigroup) %>>%
    mimicPastDF("season", season, nontime.index.colnames="Location") %>>%
    magrittr::extract2("baseline")
  n.weeks.in.season = lastWeekNumber(season, 3L)
  is.inseason = usa_flu_inseason_flags(n.weeks.in.season)
  forecast.time = model_week_to_time(
    model.week, usa.flu.first.week.of.season)
  max.lag = 51L
  return (list(
    season = season,
    model.week = model.week,
    epigroup = epigroup,
    source.name = "fluview",
    signal.name = "wili",
    issue = issue,
    ## epidata.df = epidata.df,
    ## epidata.history.df = epidata.history.df,
    ## nowcast.df = nowcast.df,
    epidata.dfs = list(
      fluview=epidata.df %>>% dplyr::mutate(lag.group=pmin(lag, max.lag)),
      nowcast=nowcast.df %>>% dplyr::mutate(lag.group=pmin(lag, max.lag))
    ),
    epidata.history.dfs = list(
      fluview=epidata.history.df %>>% dplyr::mutate(lag.group=pmin(lag, max.lag)),
      nowcast=nowcast.df %>>% dplyr::mutate(lag.group=pmin(lag, max.lag))
    ),
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
  ## todo g.fluview.history.dfs as separate argument or curried argument in this
  ## method and get_voxel_data
  ## epidata.df = mimicPastEpidataDF(
  ##   g.fluview.history.dfs[[epigroup]], observation.issue)
  ## xxx check that this doesn't change anything:
  epidata.df = g.fluview.current.dfs[[epigroup]]
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

source.name = "fluview"
signal.name = "wili"

current.issue.sw =
  g.fluview.current.dfs[[1L]] %>>%
  dplyr::filter(season == max(season)) %>>%
  {.[!is.na(.[[signal.name]]),]} %>>%
  dplyr::filter(model.week == max(model.week)) %>>%
  dplyr::select(season, model.week)

s.retro.seasons = seq.int(2003L,current.issue.sw[["season"]]-1L) %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
w.retro.model.weeks = (35:78) %>>%
  stats::setNames(paste0("MW",.)) %>>%
  with_dimnamesnames("Model Week")
g.epigroups = fluview.location.spreadsheet.names %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Location")
last.losocv.issue = 201039L
b.backcasters = list(
  ignorant=backfill_ignorant_backsim,
  quantile_arx_backnowcast=quantile_arx_pancaster(TRUE, 1L),
  quantile_arx_pancast=quantile_arx_pancaster(TRUE, 53L)
) %>>%
  with_dimnamesnames("Backcaster")
f.forecasters = list(
  "Delphi_Uniform"=uniform_forecast,
  "Delphi_EmpiricalBayes_PackageDefaults"=eb.sim,
  "Delphi_EmpiricalBayes_Cond4"=function(full.dat, baseline=0, max.n.sims=2000L) {
    eb.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims,
           control.list=get_eb_control_list(max.match.length=4L))
  },
  "Delphi_BasisRegression_PackageDefaults"=br.sim,
  "Delphi_DeltaDensity_PackageDefaults"=twkde.sim,
  "Delphi_MarkovianDeltaDensity_PackageDefaults"=twkde.markovian.sim,
  "Delphi_EmpiricalFutures_PackageDefaults"=empirical.futures.sim,
  "Delphi_EmpiricalTrajectories_PackageDefaults"=empirical.trajectories.sim
) %>>%
  with_dimnamesnames("Forecaster")
target_trajectory_preprocessor = flusight2016_target_trajectory_preprocessor
t.target.specs = flusight2016.target.specs %>>%
  with_dimnamesnames("Target")
m.forecast.types = flusight2016.proxy.forecast.types %>>%
  with_dimnamesnames("Type")

## Specify portions of cv_apply indexer lists corresponding to model week,
## epigroup, target:
e.ensemble.partial.weighting.scheme.wgt.indexer.lists = list(
  ## "constant-weights" = list(all=NULL, all=NULL, all=NULL),
  "target-based" = list(all=NULL, all=NULL, each=NULL),
  "target-3time-based" = list(smear=-1:1, all=NULL, each=NULL),
  "target-9time-based" = list(smear=-4:4, all=NULL, each=NULL)
) %>>% with_dimnamesnames("Ensemble weighting scheme")

## Use LOSOCV on seasons < 2010L (train on other seasons < 2010L), and
## oneahead on seasons >= 2010L (train on seasons < test season)
retro.season.indexer = list(loo_oneahead=match(2010L,s.retro.seasons))

epiproject.cache.dir = "~/files/nosync/epiforecast-epiproject/flusight-natreg-run"

source("generate-retro-and-prospective-forecasts.R")

## Output prospective forecast spreadsheets, plots:
collab.ensemble.retro.dir = "~/files/nosync/collaborative-ensemble-submission-3"
if (!dir.exists(collab.ensemble.retro.dir)) {
  dir.create(collab.ensemble.retro.dir)
}
save_spreadsheets(swgbf.retro.component.target.multicasts,
                  swg.retro.voxel.data,
                  t.target.specs, m.forecast.types,
                  collab.ensemble.retro.dir,
                  function(swg.voxel.data,s,w,...) {
                    season = swg.voxel.data[[s,w,1L]][["season"]]
                    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
                    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
                    if (season >= 2010L && !dplyr::between(week,21L,39L)) {
                      sprintf("%s/EW%02d-%d-%s.csv", ..2, week, year, ..2)
                    }
                  })
save_spreadsheets(swge.retro.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
                  swg.retro.voxel.data,
                  t.target.specs, m.forecast.types,
                  collab.ensemble.retro.dir,
                  function(swg.voxel.data,s,w,...) {
                    season = swg.voxel.data[[s,w,1L]][["season"]]
                    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
                    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
                    if (season >= 2010L && !dplyr::between(week,21L,39L)) {
                      sprintf("%s/EW%02d-%d-%s.csv", "Delphi_Stat_FewerComponentsNoBackcastNoNowcast", week, year, "Delphi_Stat_FewerComponentsNoBackcastNoNowcast")
                    }
                  })

collab.ensemble.prospective.dir = "~/files/nosync/cdc-flusight-ensemble/model-forecasts/real-time-component-models/"
if (!dir.exists(collab.ensemble.prospective.dir)) {
  dir.create(collab.ensemble.prospective.dir)
}
save_spreadsheets(swgbf.prospective.component.target.multicasts,
                  swg.prospective.voxel.data,
                  t.target.specs, m.forecast.types,
                  collab.ensemble.prospective.dir,
                  function(swg.voxel.data,s,w,...) {
                    season = swg.voxel.data[[s,w,1L]][["season"]]
                    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
                    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
                    if (season >= 2010L && !dplyr::between(week,21L,39L)) {
                      sprintf("%s/EW%02d-%d-%s.csv", ..2, week, year, ..2)
                    }
                  })
save_spreadsheets(swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
                  swg.prospective.voxel.data,
                  t.target.specs, m.forecast.types,
                  collab.ensemble.prospective.dir,
                  function(swg.voxel.data,s,w,...) {
                    season = swg.voxel.data[[s,w,1L]][["season"]]
                    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
                    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
                    if (season >= 2010L && !dplyr::between(week,21L,39L)) {
                      sprintf("%s/EW%02d-%d-%s.csv", "Delphi_Stat_FewerComponentsNoBackcastNoNowcast", week, year, "Delphi_Stat_FewerComponentsNoBackcastNoNowcast")
                    }
                  })

save_spreadsheets(
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/linlog.plots"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/linlog.plots"
)

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/linlog.plots-percent"
)

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/stat-spreadsheets",
  function(swg.voxel.data,s,w,...) {
    season = swg.voxel.data[[s,w,1L]][["season"]]
    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
    sprintf("EW%02d-%s-%s.csv", week, "Delphi-Stat", Sys.Date())
  }
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/stat-linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/stat-linlog.plots-percent"
)

save_weighting_linlog_plots(
  e.prospective.ensemble.weightsets[["target-9time-based"]],
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/stat-weighting-plots"
)
