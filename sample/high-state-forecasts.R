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
fluview.epigroup.name.mapping =
  tibble::tibble(
            abbreviation =
              state.abb %>>%
              c("AS","MP","DC","GU","PR","VI","ORD","LAX","JFK"),
            name = state.name %>>%
              c("American Samoa", "Commonwealth of the Northern Mariana Islands", "District of Columbia", "Guam", "Puerto Rico", "Virgin Islands", "Chicago", "Los Angeles", "New York City")
          ) %>>%
  dplyr::mutate(in.spreadsheet = ! abbreviation %in% c("FL","AS","MP","GU","ORD","LAX","JFK")) %>>%
  dplyr::mutate(in.first.training = ! abbreviation %in% c("FL","AS","MP","GU","ORD","LAX")) %>>%
  ## fixme JFK shouldn't have been here; but maybe we should predict everything
  ## with data available and filter the resulting spreadsheet; e.g., add FL
  dplyr::filter(in.first.training) %>>%
  dplyr::arrange(name)
fluview.all.location.epidata.names = tolower(fluview.epigroup.name.mapping[["abbreviation"]])
fluview.all.location.spreadsheet.names = fluview.epigroup.name.mapping[["name"]]

## Load in the epidata (takes a while):
fluview.all.current.dfs = fluview.all.location.epidata.names %>>%
  setNames(fluview.all.location.spreadsheet.names) %>>%
  lapply(function(fluview.location.epidata.name) {
    print(fluview.location.epidata.name)
    get_completed_fluview_state_df(fluview.location.epidata.name, usa.flu.first.week.of.season) %>>%
      dplyr::arrange(epiweek, issue) %>>%
      {
        season = .[["season"]]
        model.week = .[["model.week"]]
        issue = .[["issue"]]
        wili = .[["wili"]]
        lag = .[["lag"]]
        first.recorded.model.week = min(model.week[season==min(season) & !is.na(wili)])
        .[["issue"]][season==min(season) & is.na(wili)] <-
          issue[season==min(season) & model.week==first.recorded.model.week][[1L]]
        .[["wili"]][season==min(season) & is.na(wili)] <-
          wili[season==min(season) & model.week==first.recorded.model.week][[1L]]
        .[["lag"]][season==min(season) & is.na(wili)] <-
          lag[season==min(season) & model.week==first.recorded.model.week][[1L]]
        ## fixme constant fill-in not reasonable, especially for week-shifting
        ## methods with low numbers of seasons
        .
      }
  })
selected.location.inds =
  sapply(fluview.all.current.dfs,
         function(df) median(df[["num_providers"]], na.rm=TRUE)) %>>%
  {which(. > median(.))}
fluview.location.epidata.names = fluview.all.location.epidata.names[selected.location.inds]
fluview.location.spreadsheet.names = fluview.all.location.spreadsheet.names[selected.location.inds]
g.fluview.current.dfs =
  fluview.all.current.dfs[fluview.location.spreadsheet.names ]
g.fluview.history.dfs =
  g.fluview.current.dfs
## fluview.location.epidata.names %>>%
  ## setNames(fluview.location.spreadsheet.names) %>>%
  ## lapply(function(fluview.location.epidata.name) {
  ##   fetchEpidataHistoryDF(
  ##     "fluview", fluview.location.epidata.name, 0:51,
  ##     first.week.of.season=usa.flu.first.week.of.season,
  ##     ## force.cache.invalidation = TRUE,
  ##     cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_",fluview.location.epidata.name))
  ##   )
  ## })

epigroup.colname = "Location"

get_voxel_data = function(season, model.week, epigroup, last.losocv.issue) {
  current.epidata.history.df = g.fluview.history.dfs[[epigroup]][
    c("season","model.week","epiweek","issue","lag","wili")
  ]
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
  baseline = 0
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

signal.name = "wili" # xxx should use "ili" rather than "wili" (above as well)

get_observed_trajectory = function(season, epigroup) {
  ## Use the current issue's version of a trajectory as the "observed" (vs. a
  ## fixed issue after the season's end):
  epidata.df = g.fluview.current.dfs[[epigroup]]
  observed.trajectory = epidata.df %>>%
    {.[[signal.name]][.[["season"]]==season]}
  return (observed.trajectory)
}

current.issue.sw =
  g.fluview.current.dfs[[1L]] %>>%
  dplyr::filter(season == max(season)) %>>%
  {.[!is.na(.[[signal.name]]),]} %>>%
  dplyr::filter(model.week == max(model.week)) %>>%
  dplyr::select(season, model.week)

s.retro.seasons = seq.int(2010L,current.issue.sw[["season"]]-1L) %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
w.retro.model.weeks = (35:78) %>>%
  stats::setNames(paste0("MW",.)) %>>%
  with_dimnamesnames("Model Week")
g.epigroups = fluview.location.spreadsheet.names %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Location")
last.losocv.issue = 201739L
b.backcasters = list(
  ignorant=backfill_ignorant_backsim
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
## fixme remove Season onset
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

## Use LOSOCV on all seasons but current
retro.season.indexer = list(loo=NULL)

epiproject.cache.dir = "~/files/nosync/epiforecast-epiproject/flusight-high-state-run"

source("generate-retro-and-prospective-forecasts.R")

save_spreadsheets(
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/linlog.plots-percent"
)

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/linlog.plots-percent"
)

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/stat-spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/stat-linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/stat-linlog.plots-percent"
)

save_weighting_linlog_plots(
  e.prospective.ensemble.weightsets[["target-9time-based"]],
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  "~/files/nosync/epiforecast-epiproject/flusight-high-state-run/stat-weighting-plots"
)
