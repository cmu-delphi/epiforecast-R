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

location.forecast.values.dir = "~/files/nosync/forecast-values-storage"
spreadsheet.dir = "~/files/nosync/spreadsheet-storage"
retro.methods = list(
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
)
retro.Locations = fluview.location.spreadsheet.names
retro.epiweeks = 2010:2016 %>>%
  DatesOfSeason(usa.flu.first.week.of.season, 0L,3L) %>>%
  dplyr::combine() %>>%
  DateToYearWeekWdayDF(0L,3L) %>>%
  dplyr::filter(! week %>>% dplyr::between(21L,39L)) %>>%
  with(year*100L+week)
retro.targets = flusight2016.targets
retro.forecast.types = flusight2016.proxy.forecast.types
retro.first.noncv.epiweek = min(retro.epiweeks)

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
  fluview.location.epidata.names %>>%
  setNames(fluview.location.spreadsheet.names) %>>%
  lapply(function(fluview.location.epidata.name) {
    fetchEpidataHistoryDF(
      "fluview", fluview.location.epidata.name, 0:51,
      first.week.of.season=usa.flu.first.week.of.season,
      cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_",fluview.location.epidata.name))
    )
  })

flusight_2016_settings = function(forecast.epiweek, forecast.Location) {
  forecast.smw = yearWeekToSeasonModelWeekDF(forecast.epiweek%/%100L, forecast.epiweek%%100L,
                                             usa.flu.first.week.of.season, 3L)
  mimicked.baseline = mimicPastDF(fluview.baseline.df,
                                  "season.int", forecast.smw[["season"]],
                                  nontime.index.colnames="Location") %>>%
    with(baseline[Location==forecast.Location])
  ##
  n.weeks.in.season = lastWeekNumber(forecast.smw[["season"]], 3L)
  is.inseason = usa_flu_inseason_flags(n.weeks.in.season)
  mimicked.time.of.forecast = model_week_to_time(forecast.smw[["model.week"]],
                                                 usa.flu.first.week.of.season)
  flusight2016.settings = list(
    baseline=mimicked.baseline,
    is.inseason=is.inseason,
    target.time.of.forecast=mimicked.time.of.forecast
  )
  return (flusight2016.settings)
}

fluview_mimicked_epidata_df = function(mimicked.issue, retro.Location) {
  mimicked.fluview.cache.dir = file.path(epidata.cache.dir,"fluview_mimicked_dfs")
  if (!dir.exists(mimicked.fluview.cache.dir)) {
    dir.create(mimicked.fluview.cache.dir)
  }
  mimicked.epidata.df = fetchUpdatingResource(
    function() {
      mimicPastEpidataDF(fluview.history.l.dfs[[retro.Location]], mimicked.issue)
    },
    function(fetch.response) {
      return ()
    },
    cache.file.prefix=file.path(mimicked.fluview.cache.dir, paste0("mimicked_fluview_issue",mimicked.issue,"_",retro.Location)),
    cache.invalidation.period=as.difftime(Inf, units="weeks"),
    silent=TRUE
  )
  return (mimicked.epidata.df)
}

fluview_wili_retro_full_dat = function(retro.epiweek, retro.Location,
                                       first.noncv.retro.epiweek=retro.epiweek) {
  retro.smw = yearWeekToSeasonModelWeekDF(retro.epiweek%/%100L, retro.epiweek%%100L,
                                          usa.flu.first.week.of.season, 3L)
  retro.season = retro.smw[["season"]]
  retro.Season = paste0(retro.season,"/",retro.season+1L)
  ##
  retro.epidata.df =
    if (retro.epiweek >= first.noncv.retro.epiweek) {
      ## retro.epiweek is the same as or later than the possible fill-in week
      ## for other seasons, first.noncv.retro.epiweek; only use this later data
      fluview_mimicked_epidata_df(retro.epiweek, retro.Location) %>>%
        trimPartialPastSeasons("wili", 52L)
    } else {
      ## combine retro.epiweek's data for its season with fill-in data for other
      ## seasons from first.noncv.retro.epiweek
      first.noncv.retro.season = first.noncv.retro.epiweek %>>%
        {
          yearWeekToSeasonModelWeekDF(.%/%100L, .%%100L,
                                      usa.flu.first.week.of.season, 3L)
        } %>>%
        (season)
      dplyr::bind_rows(
               fluview_mimicked_epidata_df(first.noncv.retro.epiweek, retro.Location) %>>%
               trimPartialPastSeasons("wili", 52L) %>>%
               dplyr::filter(
                        season < first.noncv.retro.season,
                        season != retro.season
                      ),
               fluview_mimicked_epidata_df(retro.epiweek, retro.Location) %>>%
               dplyr::filter(season == retro.season)
             )
    }
  retro.full.dat =
    retro.epidata.df %>>%
    ## label and convert to list of numeric vectors:
    dplyr::mutate(Season=paste0(season,"/",season+1L)) %>>%
    {split(.[["wili"]], .[["Season"]])} %>>%
    ## reorder so that the retro.epiweek season is last:
    {c(.[names(.)!=retro.Season],
       .[names(.)==retro.Season])}
  return (retro.full.dat)
}

flusight2016retro_sim = function(sim.method, retro.epiweek, retro.Location,
                                 first.noncv.retro.epiweek=retro.epiweek,
                                 max.n.sims=2000L) {
  flusight2016.settings = flusight_2016_settings(retro.epiweek, retro.Location)
  retro.full.dat = fluview_wili_retro_full_dat(retro.epiweek, retro.Location,
                                               first.noncv.retro.epiweek=first.noncv.retro.epiweek)
  retro.sim = sim.method(retro.full.dat,
                         baseline=flusight2016.settings[["baseline"]],
                         max.n.sims=max.n.sims) %>>%
    c(list(flusight2016.settings=flusight2016.settings))
  return (retro.sim)
}

flusight2016_location_forecast_values_obj_from_sim =
  function(flusight2016.sim, targets, forecast.types, label.types=TRUE, label.bins=FALSE, ...) {
    flusight2016.settings = flusight2016.sim[["flusight2016.settings"]]
    if (label.types) {
      names(targets) <- sapply(targets, magrittr::extract2, "Target")
      names(forecast.types) <- sapply(forecast.types, magrittr::extract2, "Type")
    }
    lapply(targets, function(target) {
      target.name = target[["Target"]]
      ## print(target.name)
      flusight2016.settings = flusight2016.sim[["flusight2016.settings"]]
      target.forecast = target_forecast(
        flusight2016.sim,
        target.name=target.name,
        target.fun=target[["for_processed_trajectory"]],
        ## xxx could move preprocessing out to flusight2016_location_spreadsheet_from_sim:
        target_trajectory_preprocessor=flusight2016_target_trajectory_preprocessor,
        target.value.formatter=target[["to_string"]],
        baseline=flusight2016.settings[["baseline"]],
        is.inseason=flusight2016.settings[["is.inseason"]],
        target.time.of.forecast=flusight2016.settings[["target.time.of.forecast"]],
        compute.estimates=FALSE,
        target.info=target,
        ...
      )
      target.values = target.forecast[["target.values"]][[target.name]]
      target.weights = target.forecast[["target.weights"]]
      ##
      lapply(forecast.types, function(forecast.type) {
        forecast.value =
          forecast.type[["forecast_value_from_weighted_univals"]](
            target, target.values, target.weights,
            label.bins=label.bins,
            baseline=flusight2016.settings[["baseline"]],
            is.inseason=flusight2016.settings[["is.inseason"]],
            target.time.of.forecast=flusight2016.settings[["target.time.of.forecast"]],
            ...
          )
        forecast.value
      })
    }) %>>%
      {
        list(
          ## todo add targets and forecast types here?
          flusight2016.settings=flusight2016.settings,
          location.forecast.values=.
        )
      }
  }

flusight2016_location_spreadsheet_obj_from_location_forecast_values_obj =
  function(location.forecast.values.obj, targets, forecast.types, flusight2016.settings) {
    flusight2016.settings = location.forecast.values.obj[["flusight2016.settings"]]
    location.forecast.values = location.forecast.values.obj[["location.forecast.values"]]
    lapply(seq_along(targets), function(target.i) {
      target = targets[[target.i]]
      lapply(seq_along(forecast.types), function(forecast.type.i) {
        forecast.type = forecast.types[[forecast.type.i]]
        forecast.value = location.forecast.values[[target.i]][[forecast.type.i]]
        subspreadsheet_from_forecast_value(
          forecast.value, target, forecast.type,
          baseline=flusight2016.settings[["baseline"]],
          is.inseason=flusight2016.settings[["is.inseason"]],
          target.time.of.forecast=flusight2016.settings[["target.time.of.forecast"]]
        )
      }) %>>%
        dplyr::bind_rows()
    }) %>>%
      dplyr::bind_rows() %>>%
      {
        list(
          ## todo add targets and forecast types here?
          flusight2016.settings=flusight2016.settings,
          location.spreadsheet=.
        )
      }
  }

flusight2016_location_forecast_values_obj_from_location_spreadsheet_obj =
  function(location.spreadsheet.obj, targets, forecast.types, flusight2016.settings) {
    flusight2016.settings = location.spreadsheet.obj[["flusight2016.settings"]]
    location.spreadsheet = location.spreadsheet.obj[["location.spreadsheet"]]
    names(targets) <- sapply(targets, magrittr::extract2, "Target")
    names(forecast.types) <- sapply(forecast.types, magrittr::extract2, "Type")
    list(
      flusight2016.settings=flusight2016.settings,
      location.forecast.values=
        location.spreadsheet %>>%
        ## split rows into a list of df's based on Target column, maintaining
        ## spreadsheet's ordering of Target values:
        split(.[["Target"]] %>>% factor(unique(.))) %>>%
        lapply(function(target.df) {
          target = targets[[target.df[["Target"]][1L]]]
          target.df %>>%
            dplyr::select(-Target) %>>%
            ## split rows into a list of multiple df's based on Type column,
            ## maintaining spreadsheet's ordering of Type values:
            split(.[["Type"]] %>>% factor(unique(.))) %>>%
            lapply(function(target.type.df) {
              forecast.type = forecast.types[[target.type.df[["Type"]][1L]]]
              forecast.type[["forecast_value_from_Value"]](
                target.type.df[["Value"]], target,
                baseline=flusight2016.settings[["baseline"]],
                is.inseason=flusight2016.settings[["is.inseason"]],
                target.time.of.forecast=flusight2016.settings[["target.time.of.forecast"]]
              )
            })
        })
    )
  }

flusight2016retro_location_forecast_values_obj =
  function(sim.method.name, retro.epiweek, forecast.Location,
           targets, forecast.types,
           first.noncv.retro.epiweek=retro.epiweek,
           label.types=TRUE, label.bins=FALSE,
           force.recalc=FALSE, ...) {
    method.location.forecast.values.dir = file.path(location.forecast.values.dir, sim.method.name)
    if (!dir.exists(method.location.forecast.values.dir)) {
      dir.create(method.location.forecast.values.dir)
    }
    fetchUpdatingResource(
      function() {
        set.seed(42L)
        retro.methods[[sim.method.name]] %>>%
          flusight2016retro_sim(retro.epiweek, forecast.Location,
                                first.noncv.retro.epiweek=first.noncv.retro.epiweek) %>>%
          flusight2016_location_forecast_values_obj_from_sim(
            targets, forecast.types,
            label.types=label.types, label.bins=label.bins,
            ...
          )
      },
      function(fetch.response) {},
      cache.file.prefix=file.path(method.location.forecast.values.dir, paste0(sim.method.name,"_location.forecast.values.obj_issue",retro.epiweek,"_",forecast.Location,"_retro",first.noncv.retro.epiweek)),
      cache.invalidation.period=as.difftime(Inf, units="weeks"),
      force.cache.invalidation=force.recalc,
      silent=TRUE
    )
  }

if (!dir.exists(location.forecast.values.dir)) {
  dir.create(location.forecast.values.dir)
}

## sample.location.forecast.values.obj =
##   "Delphi_EmpiricalTrajectories_PackageDefaults" %>>%
##   ## "Delphi_EmpiricalFutures_PackageDefaults" %>>%
##   ## "Delphi_EmpiricalBayes_PackageDefaults" %>>%
##   flusight2016retro_location_forecast_values_obj(201652L, "US National", retro.targets, retro.forecast.types, first.noncv.epiweek=retro.first.noncv.epiweek, force.recalc=TRUE) %>>%
##   {.}
## sample.location.spreadsheet.obj =
##   sample.location.forecast.values.obj %>>%
##   flusight2016_location_spreadsheet_obj_from_location_forecast_values_obj(
##     flusight2016.targets, flusight2016.proxy.forecast.types
##   )
## recov.sample.location.forecast.values.obj =
##   sample.location.spreadsheet.obj %>>%
##   flusight2016_location_forecast_values_obj_from_location_spreadsheet_obj(
##     flusight2016.targets, flusight2016.proxy.forecast.types
##   )
## all.equal(sample.location.forecast.values.obj,
##           recov.sample.location.forecast.values.obj)
## all.equal(sample.location.spreadsheet.obj,
##           recov.sample.location.forecast.values.obj %>>%
##           flusight2016_location_spreadsheet_obj_from_location_forecast_values_obj(
##             flusight2016.targets, flusight2016.proxy.forecast.types
##           ))

## sample.location.spreadsheet.obj %>>%
##   (location.spreadsheet) %>>%
##   dplyr::filter(Type=="Bin") %>>%
##   ## ggplot2::ggplot(ggplot2::aes(Bin_start_incl, weight=as.numeric(Value))) +
##   ggplot2::ggplot(ggplot2::aes(factor(Bin_start_incl, unique(Bin_start_incl)), weight=as.numeric(Value))) +
##   ggplot2::facet_grid(factor(Target, unique(Target)) ~ Unit, scales="free_x") +
##   ggplot2::geom_bar()

system.time({
## tryCatch({
##   Rprof("Rprof_file")
if (!dir.exists(spreadsheet.dir)) {
  dir.create(spreadsheet.dir)
}
invisible(lapply(seq_along(retro.methods), function(sim.method.i) {
  sim.method.name = names(retro.methods)[[sim.method.i]]
  sim.method = retro.methods[[sim.method.i]]
  method.spreadsheet.dir = file.path(spreadsheet.dir, sim.method.name)
  if (!dir.exists(method.spreadsheet.dir)) {
    dir.create(method.spreadsheet.dir)
  }
  parallel::mclapply(seq_along(retro.epiweeks), function(retro.epiweek.i) {
  ## lapply(seq_along(retro.epiweeks), function(retro.epiweek.i) {
  ## lapply(seq_along(retro.epiweeks)[1:10], function(retro.epiweek.i) {
    retro.epiweek = retro.epiweeks[[retro.epiweek.i]]
    ## print(retro.epiweek)
    print(paste0(sim.method.name, ", ", retro.epiweek))
    ##
    forecast.spreadsheet =
      retro.Locations %>>%
      setNames(.) %>>%
      lapply(function(forecast.Location) {
        ## print(forecast.Location)
        sim.method.name %>>%
          flusight2016retro_location_forecast_values_obj(
            retro.epiweek, forecast.Location,
            retro.targets, retro.forecast.types,
            first.noncv.retro.epiweek=retro.first.noncv.epiweek
          ) %>>%
          flusight2016_location_spreadsheet_obj_from_location_forecast_values_obj(
            retro.targets, retro.forecast.types
          ) %>>%
          (location.spreadsheet)
      }) %>>%
      dplyr::bind_rows(.id="Location")
    ## print(forecast.spreadsheet)
    ##
    bin.sum.deviation.df =
      forecast.spreadsheet %>>%
      dplyr::filter(Type=="Bin") %>>%
      dplyr::group_by(Location, Target) %>>%
      dplyr::mutate(`Bin Sum Deviation`=sum(as.numeric(Value))-1) %>>%
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
    if (retro.epiweek >= retro.first.noncv.epiweek) {
      write.csv(forecast.spreadsheet,
                file.path(
                  method.spreadsheet.dir,
                  sprintf("EW%02d-%d-%s.csv",
                          retro.epiweek %% 100L,
                          retro.epiweek %/% 100L,
                          sim.method.name)),
                row.names=FALSE)
    }
    NULL
  })
}))
## }, finally = {
##   Rprof(NULL)
## })
## summaryRprof("Rprof_file")
})

stat.seasons = 2003:2016
stat.model.weeks =
  usa_flu_inseason_flags(53L) %>>%
  which() %>>%
  range() %>>%
  {seq.int(.[[1L]]-4L, .[[2L]]+4L)} %>>%
  time_to_model_week(usa.flu.first.week.of.season)
stat.Locations = retro.Locations
stat.targets = retro.targets
stat.forecast.types = retro.forecast.types
stat.method.names = names(retro.methods)
## xxx cache functions do not check for any stat.* versus retro.* setting
## differences...

component.ftclws.forecast.values =
  stat.seasons %>>%
  setNames(paste0(.,"/",.+1L)) %>>%
  lapply(function(stat.season) {
    ## print(paste0("stat.season: ",stat.season))
    stat.model.weeks %>>%
      setNames(sprintf("MW%02d", .)) %>>%
      parallel::mclapply(function(stat.model.week) {
      ## lapply(function(stat.model.week) {
        ## print(paste0("stat.model.week: ",stat.model.week))
        stat.epiweek = seasonModelWeekToYearWeekDF(stat.season, stat.model.week,
                                                   usa.flu.first.week.of.season, 3L) %>>%
          with(year*100L + week)
        stat.Locations %>>%
          setNames(.) %>>%
          lapply(function(stat.Location) {
            ## print(paste0("stat.Location: ",stat.Location))
            print(paste0(stat.Location, ", ", stat.epiweek))
            stat.method.names %>>%
              setNames(.) %>>%
              lapply(function(stat.method.name) {
                ## print(paste0("stat.method.name: ",stat.method.name))
                flusight2016retro_location_forecast_values_obj(
                  stat.method.name, stat.epiweek, stat.Location,
                  stat.targets, stat.forecast.types,
                  retro.first.noncv.epiweek
                )[["location.forecast.values"]] %>>% simplify2array()
                ## todo replace the simplify2array's with dplyr::combine's and dim-setting
              }) %>>% simplify2array()
          }) %>>% simplify2array()
      }) %>>% simplify2array()
  }) %>>% simplify2array() %>>%
  {
    names(dimnames(.)) <- c("Type", "Target", "Component", "Location", "Model Week","Season")
    .
  }

## tryCatch({
##   Rprof("Rprof_file")
retro.truth.ftlws.forecast.values =
  stat.seasons %>>%
  setNames(paste0(.,"/",.+1L)) %>>%
  lapply(function(stat.season) {
    ## print(paste0("stat.season: ",stat.season))
    stat.model.weeks %>>%
      setNames(sprintf("MW%02d", .)) %>>%
      parallel::mclapply(function(stat.model.week) {
      ## lapply(function(stat.model.week) {
        ## print(paste0("stat.model.week: ",stat.model.week))
        stat.epiweek = seasonModelWeekToYearWeekDF(stat.season, stat.model.week,
                                                   usa.flu.first.week.of.season, 3L) %>>%
          with(year*100L + week)
        n.weeks.next.season = lastWeekNumber(stat.season+1L, 3L)
        first.retro.forecast.time.next.season = which(usa_flu_inseason_flags(n.weeks.next.season))[[1L]]
        first.retro.forecast.mw.next.season = time_to_model_week(first.retro.forecast.time.next.season, usa.flu.first.week.of.season)
        retro.truth.epiweek = seasonModelWeekToYearWeekDF(stat.season+1L, first.retro.forecast.mw.next.season,
                                                          usa.flu.first.week.of.season, 3L) %>>%
          with(year*100L + week)
        stat.Locations %>>%
          setNames(.) %>>%
          lapply(function(stat.Location) {
            ## print(paste0("stat.Location: ",stat.Location))
            print(paste0(stat.Location, ", ", stat.epiweek))
            retro.truth.trajectory = fluview_wili_retro_full_dat(retro.truth.epiweek, stat.Location) %>>%
              magrittr::extract2(length(.)-1L)
            retro.truth.sim = list(
              ys=as.matrix(retro.truth.trajectory),
              weights=1,
              control.list=list(model="retro.truth"),
              flusight2016.settings=flusight_2016_settings(stat.epiweek, stat.Location)
            ) %>>% structure(class="sim")
            retro.truth.location.forecast.values =
              flusight2016_location_forecast_values_obj_from_sim(
                retro.truth.sim, stat.targets, stat.forecast.types,
                uniform.pseudoweight.total=0, smooth.sim.targets=FALSE
              )[["location.forecast.values"]] %>>%
              simplify2array()
            return (retro.truth.location.forecast.values)
          }) %>>% simplify2array()
      }) %>>% simplify2array()
  }) %>>% simplify2array() %>>%
  {
    names(dimnames(.)) <- c("Type", "Target", "Location", "Model Week","Season")
    .
  }
## }, finally = {
##   Rprof(NULL)
## })
## summaryRprof("Rprof_file")

## todo use more efficient representations than many files stored on disk and
## lots of lists in memory

component.lwsctf.forecast.values =
  aperm(component.ftclws.forecast.values, c(4:6,3:1))
retro.truth.lwstf.forecast.values =
  aperm(retro.truth.ftlws.forecast.values, c(3:5,2:1))

component.lwsctf.evaluations =
  seq_along(stat.forecast.types) %>>%
  setNames(names(stat.forecast.types)) %>>%
  sapply(function(stat.forecast.type.i) {
    stat.forecast.type = stat.forecast.types[[stat.forecast.type.i]]
    seq_along(stat.targets) %>>%
      setNames(names(stat.targets)) %>>%
      sapply(function(stat.target.i) {
        stat.target = stat.targets[[stat.target.i]]
        component.lwsc.forecast.values.tf = component.lwsmtf.forecast.values[,,,,stat.target.i,stat.forecast.type.i, drop=TRUE]
        Map(stat.forecast.type[["evaluate_forecast_value"]],
            component.lwsc.forecast.values.tf,
            retro.truth.lwstf.forecast.values[,,,stat.target.i,stat.forecast.type.i, drop=TRUE]
            ) %>>%
          {
            dim(.) <- dim(component.lwsc.forecast.values.tf)
            dimnames(.) <- dimnames(component.lwsc.forecast.values.tf)
            mode(.) <- "numeric"
            .
          }
      }, simplify="array")
  }, simplify="array") %>>%
  {
    names(dimnames(.))[5:6] <- c("Target", "Type")
    .
  }

## component.lwsctf.evaluations[,,match(2003:2009, stat.seasons),,,] %>>%
##   reshape2::melt(value.name="Evaluation") %>>% tibble::as_tibble() %>>%
##   dplyr::group_by(`Model Week`, Component, Target, Type) %>>%
##   ## (na.rm=TRUE for Season onset Point predictions)
##   dplyr::summarize(Evaluation = mean(Evaluation, na.rm=TRUE)) %>>%
##   dplyr::ungroup() %>>%
##   ggplot2::ggplot(ggplot2::aes(`Model Week`, Evaluation, colour=Component, group=Component)) +
##   ggplot2::facet_grid(Type ~ Target, scales="free_y") +
##   ggplot2::geom_line()

## component.lwsctf.evaluations[,,match(2003:2009, stat.seasons),,,] %>>%
##   reshape2::melt(value.name="Evaluation") %>>% tibble::as_tibble() %>>%
##   dplyr::filter(Target=="1 wk ahead", Type=="Bin") %>>%
##   dplyr::mutate(Evaluation=log(0.1/131+0.9*exp(Evaluation))) %>>%
##   dplyr::group_by(Location, `Model Week`, Season, Component) %>>%
##   ## (na.rm=TRUE for Season onset Point predictions)
##   dplyr::summarize(Evaluation = mean(Evaluation, na.rm=TRUE)) %>>%
##   dplyr::ungroup() %>>%
##   ggplot2::ggplot(ggplot2::aes(`Model Week`, Evaluation, colour=Component, group=Component)) +
##   ggplot2::facet_grid(Location ~ Season, scales="free_y") +
##   ggplot2::geom_line()

## component.lwsctf.evaluations[,,match(2003:2009, stat.seasons),,,] %>>%
##   reshape2::melt(value.name="Evaluation") %>>% tibble::as_tibble() %>>%
##   dplyr::filter(Target=="1 wk ahead", Type=="Bin") %>>%
##   dplyr::mutate(Evaluation=log(0.1/131+0.9*exp(Evaluation))) %>>%
##   ## dplyr::group_by(Location, `Model Week`, Season, Component) %>>%
##   ## ## (na.rm=TRUE for Season onset Point predictions)
##   ## dplyr::summarize(Evaluation = mean(Evaluation, na.rm=TRUE)) %>>%
##   ## dplyr::ungroup() %>>%
##   ggplot2::ggplot(ggplot2::aes(Evaluation, fill=Component, group=Component)) +
##   ggplot2::facet_wrap(~ Component) +
##   ggplot2::geom_histogram()

## ## component.ftclws.forecast.values[[
## retro.truth.ftlws.forecast.values[[
##                                     "Bin", "1 wk ahead",
##                                     ## "Delphi_BasisRegression_PackageDefaults",
##                                     ## "Delphi_DeltaDensity_PackageDefaults",
##                                     "US National", "MW50", "2016/2017"
##                                     ]] %>>%
##   plot(type="l")

stat.component.clwstf.coefs =
  ## component.lwsctf.forecast.values[1:6,1:2,,1:4,1:3,1:2,drop=FALSE] %>>%
  component.lwsctf.forecast.values %>>%
  ## component.lwsctf.forecast.values[,,,,,"Bin",drop=FALSE] %>>%
  cv_apply(list(all=NULL,
                smear=-4:4,
                loo_oneahead=
                  retro.first.noncv.epiweek %>>%
                  {yearWeekToSeasonModelWeekDF(.%/%100L, .%%100L, usa.flu.first.week.of.season, 3L)} %>>%
                  (season) %>>%
                  match(stat.seasons),
                all=NULL,
                each=NULL,
                each=NULL
                ),
           parallel_dim_i=2L,
           function(train, test) {
             Target = dimnames(train)[["Target"]]
             target = stat.targets[[Target]]
             Type = dimnames(train)[["Type"]]
             forecast.type = stat.forecast.types[[Type]]
             component.dim.i = match("Component", names(dimnames(train)))
             instance.method.forecast.values.listmat =
               train %>>% R.utils::wrap(list(seq_len(length(dim(train)))[-component.dim.i],
                                             component.dim.i))
             instance.observation.values.list =
               do.call(`[`, c(list(retro.truth.lwstf.forecast.values), dimnames(train)[-component.dim.i])) %>>%
               magrittr::extract(seq_along(.))
             forecast.type[["fit_ensemble_coefs"]](
               instance.method.forecast.values.listmat, instance.observation.values.list,
               prod(dim(train)[match(c("Location","Season"),names(dimnames(train)))]),
               "Delphi_Uniform"
             )
           }) %>>%
  {
    old.component.dim.i = which(names(dimnames(.)) == "Component")
    stopifnot(length(old.component.dim.i)==1L && old.component.dim.i!=1L)
    structure(.,
              dim=dim(.)[-old.component.dim.i],
              dimnames=c(setNames(dimnames(.)[1L], "Component"),
                         dimnames(.)[-c(1L,old.component.dim.i)])
              )
  } %>>%
  {.}

## stat.component.clwstf.coefs %>>%
##   reshape2::melt(value.name="Coefficient") %>>%
##   tibble::as_tibble() %>>%
##   dplyr::group_by(Target, `Model Week`, Component) %>>%
##   dplyr::summarize(`Mean Coefficient`=mean(Coefficient)) %>>%
##   dplyr::ungroup() %>>%
##   dplyr::mutate(`Model Week`=`Model Week` %>>%
##                   stringr::str_sub(3L) %>>%
##                   as.integer()
##                 ) %>>%
##   ggplot2::ggplot(ggplot2::aes(`Model Week`, `Mean Coefficient`, colour=Component, group=Component)) +
##   ggplot2::facet_wrap(~ Target) +
##   ggplot2::geom_line() +
##   ggplot2::scale_colour_brewer(palette="Set1") +
##   ggplot2::theme(panel.background=ggplot2::element_rect(fill="grey"))

## apply(stat.component.clwstf.coefs, c(1L,5:6), mean)

## saveRDS(component.ftclws.forecast.values, "~/files/nosync/component.ftclws.forecast.values.rds")
## saveRDS(retro.truth.ftlws.forecast.values, "~/files/nosync/retro.truth.ftlws.forecast.values.rds")
## saveRDS(stat.component.clwstf.coefs, "~/files/nosync/stat.component.clwstf.coefs.rds")

if (any(sapply(dimnames(component.lwsctf.forecast.values), is.null))) {
  stop("component.lwsctf.forecast.values must have all dimnames filled in")
}

stat.lwstf.forecast.values =
  component.lwsctf.forecast.values %>>%
  cv_apply(list(each=NULL, each=NULL, each=NULL, all=NULL, each=NULL, each=NULL),
           parallel_dim_i=2L,
           function(train, test) {
             slice.component.lwsctf.forecast.values = train
             ## stopifnot(all(dim(component.forecast.values)[names(dimnames(component.forecast.values))!="Component"]==1L))
             slice.component.clwstf.forecast.values =
               aperm(slice.component.lwsctf.forecast.values, c(4L,1:3,5:6))
             coef.inds.list = dimnames(slice.component.clwstf.forecast.values)
             for (dim.i in seq_along(coef.inds.list)) {
               dim.stat.names = dimnames(stat.component.clwstf.coefs)[[dim.i]]
               if (length(dim.stat.names)==1L && dim.stat.names=="all") {
                 coef.inds.list[[dim.i]] <- rep("all", length(coef.inds.list[[dim.i]]))
               }
             }
             slice.stat.component.clwstf.coefs = do.call(`[`, c(list(stat.component.clwstf.coefs), coef.inds.list))
             ## todo prevent these tiny negative weights from appearing earlier,
             ## and rename coefs back to weights:
             ## if (any(slice.stat.component.clwstf.coefs < 0))
             ##   print(slice.stat.component.clwstf.coefs)
             slice.stat.component.clwstf.coefs <- pmax(0, slice.stat.component.clwstf.coefs)
             slice.stat.component.clwstf.coefs <-
               slice.stat.component.clwstf.coefs %>>% magrittr::divide_by(sum(.))
             slice.component.clwstf.forecast.values %>>%
               do.call(what=rbind) %>>%
               ## fixme should sub in fallback's beforehand to avoid NA's instead
               ## of using na.rm=TRUE, so can match conditions under which the
               ## coefficients were optimized
               matrixStats::colWeightedMeans(slice.stat.component.clwstf.coefs, na.rm=TRUE) %>>%
               ## xxx this shoehorns all forecast.value's into numerics, even
               ## though other code attempts to keep the forecast.value type a
               ## function of the target and may assume this other type.
               dplyr::coalesce(NA_real_) %>>%
               list()
           }) %>>%
  structure(
    dim=dim(.)[names(dimnames(.)) != "Component"],
    dimnames=dimnames(.)[names(dimnames(.)) != "Component"]
  )

stat.ftlws.forecast.values = aperm(stat.lwstf.forecast.values, c(5:4,1:3))

ensemble.method.name = "Delphi_Stat_FewerComponentsNoBackcastNoNowcast"
method.spreadsheet.dir = file.path(spreadsheet.dir, ensemble.method.name)
if (!dir.exists(method.spreadsheet.dir)) {
  dir.create(method.spreadsheet.dir)
}
stat.ftlws.forecast.values %>>%
  cv_apply(list(all=NULL, all=NULL, all=NULL, each=NULL, each=NULL),
           parallel_dim_i=4L,
           function(stat.ftl.forecast.values.for.ws, same.thing) {
             stat.season =
               dimnames(stat.ftl.forecast.values.for.ws)[["Season"]] %>>%
               stringr::str_replace("^(\\d+)/\\d+$", "\\1") %>>%
               as.integer()
             stat.model.week =
               dimnames(stat.ftl.forecast.values.for.ws)[["Model Week"]] %>>%
               stringr::str_sub(3L) %>>%
               as.integer()
             stat.epiweek = seasonModelWeekToYearWeekDF(stat.season, stat.model.week,
                                                        usa.flu.first.week.of.season, 3L) %>>%
               with(year*100L + week)
             if (stat.epiweek %in% retro.epiweeks) {
               retro.epiweek = stat.epiweek
               print(retro.epiweek)
               forecast.spreadsheet =
                 retro.Locations %>>%
                 setNames(.) %>>%
                 lapply(function(retro.Location) {
                   flusight2016.settings = flusight_2016_settings(retro.epiweek, retro.Location)
                   ## retro.targets %>>%
                   ##   lapply(function(retro.target) {
                   ##     Target = retro.target[["Target"]]
                   ##     retro.forecast.types %>>%
                   ##       lapply(function(retro.forecast.type) {
                   ##         Type = retro.forecast.type[["Type"]]
                   ##         location.forecast.values.obj = list(
                   ##           flusight2016.settings=flusight2016.settings,
                   ##           location.forecast.values=stat.ftl.forecast.values.for.ws[[Type,Target,retro.Location,1L,1L]]
                   ##         )
                   ##       }) %>>% dplyr::bind_rows(.id="Type")
                   ##   }) %>>% dplyr::bind_rows(.id="Target")
                   location.forecast.values.obj = list(
                     flusight2016.settings=flusight2016.settings,
                     location.forecast.values=stat.ftl.forecast.values.for.ws[,,retro.Location,,] %>>%
                       ## to nested (Target-)list of (Type-)lists of forecast values:
                       apply(2L, identity)
                   )
                   location.spreadsheet.obj = flusight2016_location_spreadsheet_obj_from_location_forecast_values_obj(
                     location.forecast.values.obj,
                     retro.targets, retro.forecast.types,
                     flusight2016.settings
                   )
                   location.spreadsheet.obj[["location.spreadsheet"]]
                 }) %>>% dplyr::bind_rows(.id="Location")
               ## print(forecast.spreadsheet)
               bin.sum.deviation.df =
                 forecast.spreadsheet %>>%
                 dplyr::filter(Type=="Bin") %>>%
                 dplyr::group_by(Location, Target) %>>%
                 dplyr::mutate(`Bin Sum Deviation`=sum(as.numeric(Value))-1) %>>%
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
               ## if (any(is.na(as.numeric(forecast.spreadsheet[["Value"]])))) {
               ##   print(forecast.spreadsheet %>>% dplyr::filter(is.na(as.numeric(Value))))
               ##   stop ("All ensemble Value's should be non-NA.") # or maybe okay
               ## }
               write.csv(forecast.spreadsheet,
                         file.path(
                           method.spreadsheet.dir,
                           sprintf("EW%02d-%d-%s.csv",
                                   retro.epiweek %% 100L,
                                   retro.epiweek %/% 100L,
                                   ensemble.method.name)),
                         row.names=FALSE)
             }
           })


## todo better dataset representation... list of data sources (history df's? ilinet, fluview baselines, metadata?, in.season, ...) and auxiliary information indexed in a uniform way for location and time
## xxx exclude 2009/2010 (test & training) as is typical?
## xxx lost opportunities/experimentation on earlier seasons?
## xxx Method vs. Component
