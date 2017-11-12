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
fluview.first.two.current.dfs = fluview.all.location.epidata.names %>>%
  setNames(fluview.all.location.spreadsheet.names) %>>%
  magrittr::extract(1:2) %>>%
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
fluview.location.epidata.names = fluview.all.location.epidata.names[1:2]
fluview.location.spreadsheet.names = fluview.all.location.spreadsheet.names[1:2]
g.fluview.current.dfs =
  fluview.first.two.current.dfs[fluview.location.spreadsheet.names ]
g.fluview.history.dfs =
  g.fluview.current.dfs

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

current.issue.sw =
  fluview.first.two.current.dfs[[1L]] %>>%
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
g.epigroups = fluview.all.location.spreadsheet.names[1:2] %>>%
  stats::setNames(.) %>>%
  with_dimnamesnames("Location")
last.losocv.issue = 201739L
t.target.specs = flusight2016.target.specs %>>%
  with_dimnamesnames("Target")
m.forecast.types = flusight2016.proxy.forecast.types %>>%
  with_dimnamesnames("Type")

desired.weighting.scheme.name = "target-9time-based"
desired.Targets = names(t.target.specs) %>>% magrittr::extract(. != "Season onset")
desired.Types = names(m.forecast.types)

swg.prospective.voxel.data = map_join(
  get_voxel_data,
  s.prospective.seasons, w.prospective.model.weeks, g.epigroups,
  last.losocv.issue)

## map.join.df.result = map_join_(
map_join_(
  arraylike.args=c(
    list(no_join(swg.prospective.voxel.data)),
    named_array_to_name_arrayvecs(swg.prospective.voxel.data)[-3L]
  ),
  f=function(swg.voxel.data, s,w) {
    in.filename = paste(s, w, desired.weighting.scheme.name, sep=".") %>>%
      stringr::str_replace_all("/","-") %>>%
      paste0(".csv")
    low.file = file.path("~/files/nosync/epiforecast-epiproject/flusight-low-state-run/stat-spreadsheets",in.filename)
    high.file = file.path("~/files/nosync/epiforecast-epiproject/flusight-high-state-run/stat-spreadsheets",in.filename)
    low.df = readr::read_csv(low.file)
    high.df = readr::read_csv(high.file)
    combined.df = dplyr::bind_rows(low.df, high.df) %>>%
      dplyr::filter(Target %in% desired.Targets) %>>%
      dplyr::filter(Type %in% desired.Types) %>>%
      dplyr::right_join(fluview.epigroup.name.mapping %>>%
                        dplyr::filter(in.spreadsheet) %>>%
                        dplyr::transmute(Location=name),
                        by="Location") %>>%
      dplyr::arrange(Location)
    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
    out.dir = "~/files/nosync/epiforecast-epiproject/flusight-state-run/stat-spreadsheets"
    if (!dir.exists(out.dir)) {
      dir.create(out.dir, recursive=TRUE)
    }
    out.filename = sprintf("EW%02d-%s-StateILI-%s.csv", week, "Delphi-Stat", Sys.Date())
    readr::write_csv(combined.df, file.path(out.dir, out.filename))
    NULL
  },
  lapply_variant=lapply, shuffle=FALSE, show.progress=FALSE)

## spreadsheet.to.check = map.join.df.result[[1L]] %>>% readr::type_convert()
## spreadsheet.template = readr::read_csv("~/files/nosync/epiforecast-epiproject/stateili_submission_template_1718_v3.csv")
## class(spreadsheet.to.check[["Value"]])==class(spreadsheet.template[["Value"]])
## identical(sapply(readr::type_convert(spreadsheet.to.check), class), sapply(spreadsheet.template, class))
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
