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

options(mc.cores=parallel::detectCores()-1L)

## Simulate installing and loading the package:
## pkgbuild::clean_dll("../epiforecast") # only need to run if ../epiforecast is shared and was last loaded on another computer or with a different version of Rcpp
devtools::document("../epiforecast") # make sure documentation is ready (side effect: loads exported and non-exported parts of the package)
## devtools::unload("../epiforecast") # unload everything to get rid of the non-exported parts of package
## devtools::load_all("../epiforecast", export_all=FALSE) # reload only exported parts of package

## devtools::unload("../epiforecast") # unload everything to get rid of the non-exported parts of package
## devtools::load_all("../epiforecast")

## fluview.location.epidata.names = c(paste0("hhs",1:10), "nat")
## fluview.location.spreadsheet.names = c(paste0("HHS Region ",1:10), "US National")
fluview.location.epidata.names = c("nat", paste0("hhs",1:10))
fluview.location.spreadsheet.names = c("US National", paste0("HHS Region ",1:10))

spreadsheet.dir = "~/files/nosync/spreadsheet-storage"
sim.methods = list(
  ## "Delphi_EmpiricalBayes_PackageDefaults"=eb.sim,
  "Delphi_EmpiricalBayes_Cond4"=function(full.dat, baseline=0, max.n.sims=2000L) {
    eb.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims,
           control.list=get_eb_control_list(max.match.length=4L))
  },
  ## "Delphi_BasisRegression_PackageDefaults"=br.sim,
  ## "Delphi_DeltaDensity_PackageDefaults"=twkde.sim
  "Delphi_MarkovianDeltaDensity_PackageDefaults"=twkde.markovian.sim,
  "Delphi_EmpiricalFutures_PackageDefaults"=empirical.futures.sim,
  "Delphi_EmpiricalTrajectories_PackageDefaults"=empirical.futures.sim
)
targets = flusight2016.targets
forecast.types = flusight2016.proxy.forecast.types
forecast.Locations = fluview.location.spreadsheet.names
forecast.epiweeks = 2010:2016 %>>%
  DatesOfSeason(usa.flu.first.week.of.season, 0L,3L) %>>%
  dplyr::combine() %>>%
  DateToYearWeekWdayDF(0L,3L) %>>%
  dplyr::filter(! week %>>% dplyr::between(21L,39L)) %>>%
  with(year*100L+week)

## Load in the epidata (takes a while):
epidata.cache.dir = "~/.epiforecast-cache"
if (!dir.exists(epidata.cache.dir)) {
  dir.create(epidata.cache.dir)
}
fluview.baseline.info = fetchUpdatingResource(
  function() {
    return (list(
      LICENSE=RCurl::getURL("https://raw.githubusercontent.com/cdcepi/FluSight-forecasts/master/LICENSE"),
      wILI_Baseline=read.csv(textConnection(RCurl::getURL("https://raw.githubusercontent.com/cdcepi/FluSight-forecasts/master/wILI_Baseline.csv")), row.names=1L, check.names=FALSE, stringsAsFactors=FALSE)
    ))
  },
  function(fetch.response) {
    return ()
  },
  cache.file.prefix=file.path(epidata.cache.dir,"fluview_baselines"),
  cache.invalidation.period=as.difftime(1L, units="weeks")
)
cat("LICENSE for wILI_Baseline.csv:")
cat(fluview.baseline.info[["LICENSE"]])
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
  dplyr::mutate(season.int=Season %>>%
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
    cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_",fluview.location.epidata.name,"_",Sys.Date()))
  )
})
fluview.history.l.dfs =
  fetchUpdatingResource(
    function() {
      return (
        fluview.location.epidata.names %>>%
        setNames(fluview.location.spreadsheet.names) %>>%
        lapply(function(fluview.location.epidata.name) {
          fetchEpidataHistoryDF(
            "fluview", fluview.location.epidata.name, 0:51,
            first.week.of.season=usa.flu.first.week.of.season,
            cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_",fluview.location.epidata.name))
          )
        })
      )
    },
    function(fetch.response) {
      return ()
    },
    cache.file.prefix=file.path(epidata.cache.dir,"fluview.history.l.dfs"),
    cache.invalidation.period=as.difftime(1L, units="weeks")
  )

flusight2016_sim = function(sim.method, forecast.epiweek, forecast.Location, max.n.sims=2000L) {
  forecast.smw = tibble::tibble(year=forecast.epiweek%/%100L,
                                week=forecast.epiweek%%100L) %>>%
    yearWeekDFToSeasonModelWeekDF(usa.flu.first.week.of.season, 3L)
  mimicked.epidata.df = mimicPastEpidataDF(fluview.history.l.dfs[[forecast.Location]], forecast.epiweek)
  mimicked.full.dat = mimicked.epidata.df %>>%
    trimPartialPastSeasons("wili", 52L) %>>%
    dplyr::mutate(Season=paste0(season,"/",season+1L)) %>>%
    {split(.[["wili"]], .[["Season"]])}
  mimicked.baseline = mimicPastDF(fluview.baseline.df,
                                  "season.int", forecast.smw[["season"]],
                                  nontime.index.colnames="Location") %>>%
    with(baseline[Location==forecast.Location])
  ##
  partial.trajectory = dplyr::last(mimicked.full.dat)
  is.inseason = usa_flu_inseason_flags(length(partial.trajectory))
  mimicked.time.of.forecast = max(0L, which(!is.na(partial.trajectory)))
  mimicked.sim = sim.method(mimicked.full.dat,
                            max.n.sims=max.n.sims,
                            baseline=mimicked.baseline) %>>%
    c(list(flusight2016.settings=list(
             baseline=mimicked.baseline,
             is.inseason=is.inseason,
             target.time.of.forecast=mimicked.time.of.forecast
           )))
  return (mimicked.sim)
}

flusight2016_subspreadsheet = function(flusight2016.sim, targets, forecast.types, ...) {
  lapply(targets, function(target) {
    target.name = target[["Target"]]
    ## print(target.name)
    flusight2016.settings = flusight2016.sim[["flusight2016.settings"]]
    target.forecast = target_forecast(
      flusight2016.sim,
      target.name=target.name,
      target.fun=target[["for_processed_trajectory"]],
      ## xxx could move preprocessing out to flusight2016_subspreadsheet:
      target_trajectory_preprocessor=flusight2016_target_trajectory_preprocessor,
      target.value.formatter=target[["to_string"]],
      baseline=flusight2016.settings[["baseline"]],
      is.inseason=flusight2016.settings[["is.inseason"]],
      target.time.of.forecast=flusight2016.settings[["target.time.of.forecast"]],
      compute.estimates=FALSE,
      ...
    )
    target.values = target.forecast[["target.values"]][[target.name]]
    target.weights = target.forecast[["target.weights"]]
    ##
    lapply(forecast.types, function(forecast.type) {
      partial_spreadsheet_from_weighted_univals(
        forecast.type, target,
        target.values, target.weights,
        baseline=flusight2016.settings[["baseline"]],
        is.inseason=flusight2016.settings[["is.inseason"]],
        target.time.of.forecast=flusight2016.settings[["target.time.of.forecast"]],
        ...
      )
    }) %>>%
      dplyr::bind_rows()
  }) %>>%
    dplyr::bind_rows()
}

## set.seed(42L)
## sample.sim = flusight2016_sim(eb.sim, 201645L, "US National")

## sample.subspreadsheet = flusight2016_subspreadsheet(
##   sample.sim, flusight2016.targets, flusight2016.proxy.forecast.types
## )

## sample.subspreadsheet %>>%
##   dplyr::filter(Type=="Bin") %>>%
##   ## ggplot2::ggplot(ggplot2::aes(Bin_start_incl, weight=Value)) +
##   ggplot2::ggplot(ggplot2::aes(factor(Bin_start_incl, unique(Bin_start_incl)), weight=Value)) +
##   ggplot2::facet_grid(factor(Target, unique(Target)) ~ Unit, scales="free_x") +
##   ggplot2::geom_bar()

## tryCatch({
##   Rprof("Rprof_file")
if (!dir.exists(spreadsheet.dir)) {
  dir.create(spreadsheet.dir)
}
for (sim.method.i in seq_along(sim.methods)) {
  sim.method.name = names(sim.methods)[[sim.method.i]]
  sim.method = sim.methods[[sim.method.i]]
  method.spreadsheet.dir = file.path(spreadsheet.dir, sim.method.name)
  if (!dir.exists(method.spreadsheet.dir)) {
    dir.create(method.spreadsheet.dir)
  }
  parallel::mclapply(seq_along(forecast.epiweeks), function(forecast.epiweek.i) {
  ## lapply(seq_along(forecast.epiweeks), function(forecast.epiweek.i) {
    forecast.epiweek = forecast.epiweeks[[forecast.epiweek.i]]
    ## print(forecast.epiweek)
    print(paste0(sim.method.name, ", ", forecast.epiweek))
    ##
    forecast.spreadsheet =
      forecast.Locations %>>%
      setNames(.) %>>%
      lapply(function(forecast.Location) {
        ## print(forecast.Location)
        flusight2016.sim = flusight2016_sim(sim.method, forecast.epiweek, forecast.Location)
        flusight2016_subspreadsheet(flusight2016.sim,
                                    flusight2016.targets, flusight2016.proxy.forecast.types)
      }) %>>%
      dplyr::bind_rows(.id="Location")
    .GlobalEnv[["g.forecast.spreadsheet"]] <- forecast.spreadsheet
    ##
    bin.sum.deviation.df =
      forecast.spreadsheet %>>%
      dplyr::filter(Type=="Bin") %>>%
      dplyr::group_by(Location, Target) %>>%
      dplyr::mutate(`Bin Sum Deviation`=sum(Value)-1) %>>%
      dplyr::filter(abs(`Bin Sum Deviation`) > sqrt(.Machine[["double.eps"]])) %>>%
      dplyr::arrange(-abs(`Bin Sum Deviation`)) %>>%
      {.}
    if (nrow(bin.sum.deviation.df) != 0L) {
      print(bin.sum.deviation.df)
      ## bin.sum.deviation.df %>>%
      ##   dplyr::filter(Location=="US National") %>>%
      ##   dplyr::filter(Target=="Season peak week") %>>%
      ##   View()
      stop ("Bin forecast not properly normalized.")
    }
    write.csv(forecast.spreadsheet,
              file.path(
                method.spreadsheet.dir,
                sprintf("EW%02d-%d-%s.csv",
                        forecast.epiweek %% 100L,
                        forecast.epiweek %/% 100L,
                        sim.method.name)),
              row.names=FALSE)
    NULL
  })
}
## }, finally = {
##   Rprof(NULL)
## })
## summaryRprof("Rprof_file")

## fixme better dataset representation... list of data sources (history df's? ilinet, fluview baselines, metadata?, in.season, ...) and auxiliary information indexed in a uniform way for location and time
