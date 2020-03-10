## author_header begin
## Copyright (C) 2017 Logan C. Brooks, Aaron Rumack
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

epidata.cache.dir = "~/.epiforecast-cache"
if (!dir.exists(epidata.cache.dir)) {
  dir.create(epidata.cache.dir)
}

## different location naming schemes:
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

## If median number of providers changes in season, a state can flip from low to high
## Use inds to avoid this issue
selected.location.inds =
  sapply(fluview.all.current.dfs,
         function(df) median(df[["num_providers"]], na.rm=TRUE)) %>>%
  {which(. > median(.))}
#inds = fluview.all.location.epidata.names %in% c("al","az","ca","ga","il","in","ks","la","me","ma","mi","ms","nv","nm","ny","jfk","nc","oh","pa","sd","tn","tx","ut","va","wv","wi")
fluview.location.epidata.names = fluview.all.location.epidata.names[selected.location.inds]
fluview.location.spreadsheet.names = fluview.all.location.spreadsheet.names[selected.location.inds]
g.fluview.current.dfs =
  fluview.all.current.dfs

r1 = c("me","ma","vt","ct","nh","ri")
r2 = c("ny_minus_jfk","nj","pr","vi","jfk")
r3 = c("de","dc","md","pa","va","wv")
r4 = c("al","ga","fl","ms","nc","sc","ky","tn")
r5 = c("mi","mn","wi","il","in","oh")
r6 = c("tx","la","ok","ar","nm")
r7 = c("ia","ks","mo","ne")
r8 = c("mt","wy","nd","sd","co","ut")
r9 = c("ca","az","nv","hi")
r10 = c("ak","wa","or","id")
reg.list = list(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)

epigroup.colname = "Location"

g.region.fluview.current.dfs = paste0("hhs",1:10) %>>%
  setNames(paste0("HHS Region ",1:10)) %>>%
  lapply(function(fluview.location.epidata.name) {
    fetchEpidataDF(
      "fluview", fluview.location.epidata.name,
      first.week.of.season=usa.flu.first.week.of.season,
      cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_current_",fluview.location.epidata.name))
    )
  })

g.fluview.history.dfs =
  fluview.all.location.epidata.names %>>%
  setNames(fluview.all.location.spreadsheet.names) %>>%
  lapply(function(fluview.location.epidata.name) {
    tryCatch( {fetchEpidataHistoryDF(
      "fluview", fluview.location.epidata.name, 0:51,
      first.week.of.season=usa.flu.first.week.of.season,
      ## force.cache.invalidation = TRUE,
      cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_",fluview.location.epidata.name))
    )}, error=function(e){
      return (g.fluview.current.dfs[fluview.location.spreadsheet.names[which(fluview.location.epidata.names == fluview.location.epidata.name)]])
    })
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

reg.history.dfs = 1:10 %>>% lapply(function(reg.num) {
  fetchEpidataHistoryDF("fluview",
                         paste0("hhs",reg.num),0:51,
                         first.week.of.season=usa.flu.first.week.of.season,
                         cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_hhs",reg.num))
  )}) %>>%
  setNames(paste0("HHS Region ",1:10))
reg.current.dfs = 1:10 %>>% lapply(function(reg.num) {
  fetchEpidataDF(
    "fluview", paste0("hhs",reg.num),
    first.week.of.season=usa.flu.first.week.of.season,
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_current_hhs",reg.num))
  )}) %>>%
  setNames(paste0("HHS Region ",1:10))

g.epigroups = fluview.all.location.spreadsheet.names %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Location")
l.lags = 1:20 %>>% stats::setNames(.) %>>% with_dimnamesnames("Lag")
gl.stats = map_join(
  function(st,lag) {
    df = g.fluview.history.dfs[[st]]
    df1 = df[df[["lag"]]==lag,c("epiweek","wili")]
    df2 = g.fluview.current.dfs[[st]][,c("epiweek","wili")]
    df = dplyr::left_join(df1,df2,by="epiweek")
    lst = df[df[["epiweek"]]>=201740,][["wili.x"]] - df[df[["epiweek"]]>=201740,][["wili.y"]]
    return(c(mean(lst,na.rm=TRUE),sd(lst,na.rm=TRUE)))
  }, g.epigroups, l.lags)

census_weights = readRDS("census_weights.rds")
census_weights = census_weights[,-10,]
sim_backfill = function(season,reg,lag) {
  reg.current = reg.current.dfs[[reg]]
  reg.history = reg.history.dfs[[reg]]
  reg.history = reg.history[reg.history[["season"]]==season & reg.history[["lag"]]==lag,c("epiweek","wili")]
  reg.current = reg.current[reg.current[["season"]]==season,c("epiweek","wili")]
  df = dplyr::left_join(reg.history,reg.current,by="epiweek")
  reg.backfill = df[["wili.x"]] - df[["wili.y"]]
  wts = census_weights[season-2009,,reg]
  wts = wts[wts!=0] / sum(wts)
  dim(wts) = length(wts)
  reg.names = fluview.all.location.spreadsheet.names[census_weights[season-2009,,reg]>0]
  state.backfill = rep(reg.backfill,length(reg.names))
  dim(state.backfill) = c(length(reg.backfill),length(reg.names))
  state.backfill = t(state.backfill)

  # See https://arxiv.org/pdf/1607.04751.pdf
  # Didn't check proof, in practice works somewhat
  mu = as.numeric(lapply(gl.stats[reg.names,1],function(lst) { lst[[1]]}))
  sig = diag(as.numeric(lapply(gl.stats[reg.names,1],function(lst) { lst[[2]]})))
  for (i in 1:length(reg.backfill)) {
    y = MASS::mvrnorm(1,mu,sig)
    alpha = solve(wts %*% sig %*% wts, reg.backfill[[i]] - wts %*% y)
    state.backfill[,i] = y + (sig %*% wts) %*% alpha
  }
  return(list(state.backfill,reg.names,df[["epiweek"]]))
}

insert_backfill = function(ili, lag, ews, reg, st) {
  reg.hist = reg.history.dfs[[reg]]
  reg.df = reg.hist[reg.hist[["epiweek"]] %in% ews & reg.hist[["lag"]]==lag,]
  st.df = reg.df
  st.df[["region"]] = fluview.all.location.epidata.names[fluview.all.location.spreadsheet.names==st]
  st.df[["ili"]] = ili
  st.df[["wili"]] = ili
  st.df[c("num_providers","num_patients","num_age_0","num_age_1","num_age_2","num_age_3","num_age_4","num_age_5","num_ili")] = NA
  g.fluview.history.dfs[[st]] <<- rbind(g.fluview.history.dfs[[st]],st.df)
}

for (s in 2010:2016) {
  for (r in 1:10) {
    for (lag in l.lags) {
      backfill = sim_backfill(s,r,lag)
      ili.change = backfill[[1]]
      st.names = backfill[[2]]
      ews = backfill[[3]]
      dimnames(ili.change) = list(st.names,1:(dim(ili.change)[[2]]))
      for (st.name in st.names) {
        st.current = g.fluview.current.dfs[[st.name]]
        st.current = st.current[st.current[["epiweek"]] %in% ews,c("epiweek","wili")]
        st.change = data.frame(ews,ili.change[st.name,])
        names(st.change) = c("epiweek","wili.change")
        joined.df = dplyr::left_join(st.current,st.change,by="epiweek")
        insert_backfill(joined.df[["wili"]]+joined.df[["wili.change"]], lag, joined.df[["epiweek"]], r, st.name)
      }
    }
  }
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
  baseline = 0
  n.weeks.in.season = lastWeekNumber(season, 3L)
  is.inseason = usa_flu_inseason_flags(n.weeks.in.season)
  forecast.time = model_week_to_time(
    model.week, usa.flu.first.week.of.season)
  max.lag = 51L
  return (list(
    season = season,
    model.week = model.week,
    epigroup = epigroup,
    issue = issue,
    source.name = "fluview",
    signal.name = "wili",
    #epidata.df = epidata.df,
    #epidata.history.df = epidata.history.df,
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
  ignorant=backfill_ignorant_backsim,
  quantile_arx_backnowcast=quantile_arx_pancaster(TRUE, 1L),
  quantile_arx_pancast=quantile_arx_pancaster(TRUE, 53L)
) %>>%
  with_dimnamesnames("Backcaster")
f.forecasters = list(
  "Delphi_Uniform"=uniform_forecast,
#  "Delphi_EmpiricalBayes_PackageDefaults"=eb.sim,
#  "Delphi_EmpiricalBayes_Cond4"=function(full.dat, baseline=0, max.n.sims=2000L) {
#    eb.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims,
#           control.list=get_eb_control_list(max.match.length=4L))
#  },
  "Delphi_BasisRegression_PackageDefaults"=br.sim,
  "Delphi_DeltaDensity_PackageDefaults"=twkde.sim,
  "Delphi_MarkovianDeltaDensity_PackageDefaults"=twkde.markovian.sim,
  "Delphi_EmpiricalFutures_PackageDefaults"=empirical.futures.sim,
  "Delphi_EmpiricalTrajectories_PackageDefaults"=empirical.trajectories.sim
) %>>%
  with_dimnamesnames("Forecaster")
target_trajectory_preprocessor = flusight2016ilinet_target_trajectory_preprocessor
full_dat_fixup = identity # avoid methods requiring fixup for now
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

epiproject.cache.dir = "../../../epiforecast-epiproject/flusight-high-state-run"

source("generate-retro-and-prospective-forecasts.R")

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  file.path(epiproject.cache.dir,"stat-spreadsheets")
)

save_spreadsheets(
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  file.path(epiproject.cache.dir,"spreadsheets")
)

save_linlog_plots(
  target_multicast_week_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"linlog.plots-week")
)

save_linlog_plots(
  target_multicast_percent_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"linlog.plots-percent")
)

save_spreadsheets(
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  file.path(epiproject.cache.dir,"spreadsheets")
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"linlog.plots-week")
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"linlog.plots-percent")
)

save_linlog_plots(
  target_multicast_week_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"stat-linlog.plots-week")
)

save_linlog_plots(
  target_multicast_percent_plot,
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"stat-linlog.plots-percent")
)

save_weighting_linlog_plots(
  e.prospective.ensemble.weightsets[["target-9time-based"]],
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"stat-weighting-plots")
)
