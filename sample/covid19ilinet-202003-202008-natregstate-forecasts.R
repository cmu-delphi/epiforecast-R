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

## different statelike naming schemes:
fluview.epigroup.name.mapping =
  tibble::tibble(
            abbreviation =
              state.abb %>>%
              dplyr::recode('NY'='NY_MINUS_JFK') %>>%
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
fluview.all.statelike.epidata.names = tolower(fluview.epigroup.name.mapping[["abbreviation"]])
fluview.all.statelike.spreadsheet.names = fluview.epigroup.name.mapping[["name"]]

## different location naming schemes:
fluview.location.epidata.names = c("nat", paste0("hhs",1:10), fluview.all.statelike.epidata.names)
fluview.location.spreadsheet.names = c("US National", paste0("HHS Region ",1:10), fluview.all.statelike.spreadsheet.names)

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
    cache.invalidation.period=as.difftime(2L, units="days"),
    force.cache.invalidation=FALSE
)
# Manually correct for bad formatting of file
#n = names(fluview.baseline.info[["wILI_Baseline"]])[2:13]
#fluview.baseline.info[["wILI_Baseline"]] = fluview.baseline.info[["wILI_Baseline"]][,1:12]
#names(fluview.baseline.info[["wILI_Baseline"]]) = n

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

extended_season_flags = function(n.weeks.in.regular.season) {
    n.weeks.in.extended.season = n.weeks.in.regular.season + last.extended.epi.week-usa.flu.first.week.of.season+1L
    seq_len(n.weeks.in.extended.season ) %>>%
        time_to_model_week(usa.flu.first.week.of.season) %>>%
        model_week_to_epi_week(usa.flu.first.week.of.season, n.weeks.in.extended.season) %>>%
        dplyr::between(10L, 37L) %>>%
        return()
}

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
  n.weeks.in.regular.season = lastWeekNumber(season, 3L)
  ## todo change naming here
  is.inseason = extended_season_flags(n.weeks.in.regular.season)
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

source.name = "fluview"
signal.name = "wili"

## last.extended.epi.week = 37L
## last.extended.epi.week = 44L
last.extended.epi.week = 46L
epidata_df_to_chopped_trajectory_df = chop_and_extend_by_season(function(season) season*100L+usa.flu.first.week.of.season,
                                                                function(season) (season+1L)*100L+last.extended.epi.week,
                                                                function(season) season_to_Season(season, usa.flu.first.week.of.season)
                                                                )

## todo fix:
## get_observed_trajectory = function(season, epigroup) {
##   ## Use the current issue's version of a trajectory as the "observed" (vs. a
##   ## fixed issue after the season's end):
##   epidata.df = g.fluview.current.dfs[[epigroup]]
##   observed.trajectory = epidata.df %>>%
##     {.[[signal.name]][.[["season"]]==season]}
##   return (observed.trajectory)
## }

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
last.losocv.issue = 201839L
b.backcasters = list(
  "RevisionIgnorant"=backfill_ignorant_backsim,
  "RevisionIgnorantWithT2Nowcast"=backfill_ignorant_student_t2_nowcast_backcaster(n.sims=1000L),
  "QARXBackcastAllShiftsNoILINearby"=quantile_arx_pancaster(
      include.nowcast=FALSE,
      max.weeks.ahead=0L,
      lambda=1e-3, tol=1e-3,
      model.week.shift.range=NULL,
      n.sims=200L
  ),
  "QARXBacknowcastAllShiftsNoILINearby"=quantile_arx_pancaster(
      include.nowcast=FALSE,
      max.weeks.ahead=1L,
      lambda=1e-3, tol=1e-3,
      model.week.shift.range=NULL,
      n.sims=200L
  ),
  "QARXBacknowcastAllShifts"=quantile_arx_pancaster(
      include.nowcast=TRUE,
      max.weeks.ahead=1L,
      lambda=1e-3, tol=1e-3,
      model.week.shift.range=NULL,
      n.sims=200L
  ),
  "QARXPancastAllShifts"=quantile_arx_pancaster(
      include.nowcast=TRUE,
      max.weeks.ahead=53L+last.extended.epi.week-usa.flu.first.week.of.season+1L,
      lambda=1e-3, tol=1e-3,
      model.week.shift.range=NULL,
      n.sims=200L
  )
) %>>%
  with_dimnamesnames("Backcaster")
f.forecasters = list(
  "Delphi_Uniform"=uniform_forecast,
  ## "Delphi_EmpiricalBayes"=eb.sim,
  ## "Delphi_EmpiricalBayes_Cond4"=function(full.dat, baseline=0, max.n.sims=2000L) {
  ##   eb.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims,
  ##          control.list=get_eb_control_list(max.match.length=4L))
  ## },
  "Delphi_BasisRegression"=br.sim,
  "Delphi_ExtendedDeltaDensity_ExpandWindows"=function(full.dat, baseline=0, max.n.sims=1000L) {
      twkde.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims
              ## , max.shifts=c(0L,0:25,24:0)[(usa.flu.first.week.of.season+seq_len(53L+last.extended.epi.week)-1L)%%52L+1L]
              ## Holiday-aware shifts above will only use things after the holidays for modeling deltas around the end of the season; the quickest adjustment is just to include the holiday-impacted data after all:
              , max.shifts=rep(27L,53L+last.extended.epi.week-usa.flu.first.week.of.season+1L)
              ## Choose shift decay factor to weight data from the opposite time of year by 0.5:
              , shift.decay.factor=0.5^(1/27)
              ## Do not use the non-exponential sum-of-observations-so-far covariate; this is meant to relate to immunity but non-COVID-19 ILI may impart no, very weak, and/or short-term resistance to COVID-19 (although we might say the same for categories of ILI in a regular season, but to a lesser extent overall); additionally, this sum-so-far doesn't seem the most relevant when calculated over more than a season's worth of data.  Additionally, reduce the weight assigned to the exponential moving average of values covariate and tighten the weighting pattern it uses, for similar reasons.
              , tradeoff.weights=c(0.5, 0.0, 0.25, 0.25)
              , decay.factor=0.5
              , diff.decay.factor=0.5 # (same as default 0.5)
              ## Reduce re-weighting and blending heuristic weights from the default 0.1 to 0.01.  The default for the first combined with the expanded max.shifts pushes values too far toward typical values over all data.  The default for the second pushes too strongly toward typical values for the time of season.
              , uniform.weight.factor=0.01
              , y.shrink=0.01
                )
  },
  "Delphi_MarkovianDeltaDensity"=twkde.markovian.sim,
  "Delphi_EmpiricalFutures"=empirical.futures.sim,
  "Delphi_EmpiricalTrajectories"=empirical.trajectories.sim
) %>>%
  with_dimnamesnames("Forecaster")
## todo adjust:
target_trajectory_preprocessor = covid19ilinet_target_trajectory_preprocessor
full_dat_fixup = identity # avoid methods requiring fixup for now
t.target.specs = covid19ilinet.202003.202008.target.specs  %>>%
  with_dimnamesnames("Target")
m.forecast.types = covid19ilinet.forecast.types %>>%
  with_dimnamesnames("Type")

## Specify portions of cv_apply indexer lists corresponding to model week,
## epigroup, target:
e.ensemble.partial.weighting.scheme.wgt.indexer.lists = list(
  "constant-weights" = list(all=NULL, all=NULL, all=NULL),
  "target-based" = list(all=NULL, all=NULL, each=NULL)#,
  ## "target-3time-based" = list(smear=-1:1, all=NULL, each=NULL),
  ## "target-9time-based" = list(smear=-4:4, all=NULL, each=NULL)
) %>>% with_dimnamesnames("Ensemble weighting scheme")

## Use LOSOCV on all seasons but current
retro.season.indexer = list(loo=NULL)

epiproject.cache.dir = "~/files/nosync/epiforecast-epiproject/covid19ilinet-202003-202008-natregstate-run"

## source("generate-retro-and-prospective-forecasts.R")
source("gen-prospective-component-forecasts.R")




## todo natreg todos
## todo change state predictions to instead use work off of low&high state configs
