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

## Different age naming schemes:
flusurv.age.epidata.names = paste0("rate_",c("overall",paste0("age_",0:4)))
flusurv.age.spreadsheet.names = c("Overall", "0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr", "65+ yr")

## Load in the epidata (takes a while):
epidata.cache.dir = "~/.epiforecast-cache"
if (!dir.exists(epidata.cache.dir)) {
  dir.create(epidata.cache.dir)
}

flusurv.network_all.every_age_group.current.df =
  fetchEpidataDF(
    "flusurv", "network_all",
    first.week.of.season=usa.flu.first.week.of.season,
    cache.file.prefix=file.path(epidata.cache.dir,"flusurv_network_all_current"))

g.flusurv.network_all.current.dfs =
  flusurv.age.epidata.names %>>%
  stats::setNames(flusurv.age.spreadsheet.names) %>>%
  lapply(function(flusurv.age.epidata.name) {
    flusurv.network_all.every_age_group.current.df %>>%
          magrittr::extract(c(names(.)[!grepl("^rate",names(.))], flusurv.age.epidata.name)) %>>%
          dplyr::rename_(.dots=c("rate"=flusurv.age.epidata.name))
  })

flusurv.network_all.every_age_group.history.df =
  fetchEpidataHistoryDF(
    "flusurv", "network_all", 0:51,
    first.week.of.season=usa.flu.first.week.of.season,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("flusurv_network_all"))
  )

g.flusurv.network_all.history.dfs =
  flusurv.age.epidata.names %>>%
  stats::setNames(flusurv.age.spreadsheet.names) %>>%
  lapply(function(flusurv.age.epidata.name) {
    flusurv.network_all.every_age_group.history.df %>>%
          magrittr::extract(c(names(.)[!grepl("^rate",names(.))], flusurv.age.epidata.name)) %>>%
          dplyr::rename_(.dots=c("rate"=flusurv.age.epidata.name))
  })

epigroup.colname = "Location" # calling age groups "Location"'s

get_voxel_data = function(season, model.week, epigroup, last.losocv.issue) {
  current.epidata.history.df = g.flusurv.network_all.history.dfs[[epigroup]][
    c("season","model.week","epiweek","issue","lag","rate")
  ] %>>%
    dplyr::filter(!dplyr::between(epiweek%%100L, 18L,39L))
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
  ) %>>%
    dplyr::filter(!dplyr::between(epiweek%%100L, 18L,39L))
  losocv.extra.epidata.df = mimicPastEpidataDF(
    losocv.extra.history.df, last.losocv.issue
  ) %>>%
    dplyr::filter(!dplyr::between(epiweek%%100L, 18L,39L))
  epidata.df =
    dplyr::bind_rows(losocv.extra.epidata.df, forward.looking.epidata.df)
  baseline = 0
  n.weeks.in.season = lastWeekNumber(season, 3L)
  is.inseason = flusurv_inseason_flags(n.weeks.in.season)
  forecast.time = model_week_to_time(
    model.week, usa.flu.first.week.of.season)
  max.lag = 51L
  return (list(
    season = season,
    model.week = model.week,
    epigroup = epigroup,
    source.name = "flusurv.network_all",
    signal.name = "rate",
    issue = issue,
    ## epidata.df = epidata.df,
    ## epidata.history.df = epidata.history.df,
    epidata.dfs = list(
      flusurv.network_all = epidata.df %>>% dplyr::mutate(lag.group=pmin(lag, max.lag))
    ),
    epidata.history.dfs = list(
      flusurv.network_all = epidata.history.df %>>% dplyr::mutate(lag.group=pmin(lag, max.lag))
    ),
    baseline = baseline,
    target.settings = list(
      baseline = baseline,
      age.group = epigroup,
      is.inseason = is.inseason,
      forecast.time = forecast.time
    ),
    ## todo move to setting/task/dataset?:
    first.week.of.season = usa.flu.first.week.of.season
  ))
}

source.name = "flusurv.network_all"
signal.name = "rate"

get_observed_trajectory = function(season, epigroup) {
  ## Use the current issue's version of a trajectory as the "observed" (vs. a
  ## fixed issue after the season's end):
  epidata.df = g.flusurv.network_all.current.dfs[[epigroup]]
  observed.trajectory = epidata.df %>>%
    dplyr::filter(!dplyr::between(epiweek%%100L, 18L,39L)) %>>%
    {.[[signal.name]][.[["season"]]==season]}
  return (observed.trajectory)
}

current.issue.sw =
  g.flusurv.network_all.current.dfs[[1L]] %>>%
  dplyr::filter(season == max(season)) %>>%
  {.[!is.na(.[[signal.name]]),]} %>>%
  dplyr::filter(model.week == max(model.week)) %>>%
  dplyr::select(season, model.week)

s.retro.seasons = seq.int(2010L,current.issue.sw[["season"]]-1L) %>>%
  stats::setNames(paste0(.,"/",.+1L)) %>>%
  with_dimnamesnames("Season")
## fixme expand week range
w.retro.model.weeks = (40:65) %>>%
## w.retro.model.weeks = (40:69) %>>%
  stats::setNames(paste0("MW",.)) %>>%
  with_dimnamesnames("Model Week")
g.epigroups = flusurv.age.spreadsheet.names %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Location")
last.losocv.issue = 201739L
b.backcasters = list(
  ignorant=backfill_ignorant_backsim,
  quantile_arx_backcast=quantile_arx_pancaster(FALSE, 0L),
  quantile_arx_pancast=quantile_arx_pancaster(FALSE, 53L)
) %>>%
  with_dimnamesnames("Backcaster")
f.forecasters = list(
  "Delphi_Uniform"=uniform_forecast,
  ## "Delphi_EmpiricalBayes_PackageDefaults"=eb.sim,
  ## "Delphi_EmpiricalBayes_Cond4"=function(full.dat, baseline=0, max.n.sims=2000L) {
  ##   eb.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims,
  ##          control.list=get_eb_control_list(max.match.length=4L))
  ## },
  "Delphi_BasisRegression_PackageDefaults"=br.sim,
  "Delphi_DeltaDensity_RealignWindows"=function(full.dat, baseline=0, max.n.sims=1000L) {
    twkde.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims, max.shifts=c(rep(10L,1L),10:1,rep(0L,3L),1:10,rep(10L,7L)))
  },
  "Delphi_MarkovianDeltaDensity_PackageDefaults"=twkde.markovian.sim,
  "Delphi_EmpiricalFutures_PackageDefaults"=empirical.futures.sim,
  "Delphi_EmpiricalTrajectories_PackageDefaults"=empirical.trajectories.sim
) %>>%
  with_dimnamesnames("Forecaster")
target_trajectory_preprocessor = flusurv2017_target_trajectory_preprocessor
t.target.specs = flusurv2017.target.specs %>>%
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

## Use LOSOCV on all seasons but current
retro.season.indexer = list(loo=NULL)

epiproject.cache.dir = "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run"

source("generate-retro-and-prospective-forecasts.R")

save_spreadsheets(
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/linlog.plots-percent"
)

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/linlog.plots-percent"
)

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/stat-spreadsheets"
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/stat-linlog.plots-week"
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/stat-linlog.plots-percent"
)

## save_weighting_linlog_plots(
##   e.prospective.ensemble.weightsets[["target-9time-based"]],
##   swgbf.prospective.component.target.multicasts,
##   swg.prospective.voxel.data,
##   t.target.specs, m.forecast.types,
##   "~/files/nosync/epiforecast-epiproject/flusurv-network_all-run/stat-weighting-plots"
## )
## ## fixme there appears to be a bug causing this call to generate errors

## fixme no epiweek 18--21 predictions will be made!  this will break late-season forecasts if they want those weeks, but they apparently don't.  However, the forecast date range seems to go through when epiweek 17 is released --- this is okay, backcasting only.  But it may break the forecast targets... so these forecast weeks are omitted for now...
## fixme recent rates are at a lower resolution
## xxx note bugfix of twkde windows
## xxx hospitalization challenge should include backcasting

## fixme expand week range
