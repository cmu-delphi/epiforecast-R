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

source("natreg-config.R")

source("generate-retro-and-prospective-forecasts.R")

## Output prospective forecast spreadsheets, plots:
collab.ensemble.retro.dir = "../../../collaborative-ensemble-potential-submission-4"
if (!dir.exists(collab.ensemble.retro.dir)) {
  dir.create(collab.ensemble.retro.dir)
}
save_spreadsheets(swgbf.retro.component.target.multicasts[,,,"quantile_arx_backnowcast",,drop=FALSE],
                  swg.retro.voxel.data,
                  t.target.specs, m.forecast.types,
                  epigroup.colname,
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
                  epigroup.colname,
                  collab.ensemble.retro.dir,
                  function(swg.voxel.data,s,w,...) {
                    season = swg.voxel.data[[s,w,1L]][["season"]]
                    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
                    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
                    if (season >= 2010L && !dplyr::between(week,21L,39L)) {
                      sprintf("%s/EW%02d-%d-%s.csv", "Delphi_Stat_FewerComponentsNoBackcastNoNowcast", week, year, "Delphi_Stat_FewerComponentsNoBackcastNoNowcast")
                    }
                  })
save_spreadsheets(swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
                  swg.prospective.voxel.data,
                  t.target.specs, m.forecast.types,
                  epigroup.colname,
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
  swge.prospective.ensemble.target.multicasts[,,,"target-9time-based",drop=FALSE],
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  epigroup.colname,
  file.path(epiproject.cache.dir,"stat-spreadsheets"),
  function(swg.voxel.data,s,w,...) {
    season = swg.voxel.data[[s,w,1L]][["season"]]
    year = swg.voxel.data[[s,w,1L]][["issue"]] %/% 100L
    week = swg.voxel.data[[s,w,1L]][["issue"]] %% 100L
    sprintf("EW%02d-%s-%s.csv", week, "Delphi-Stat", Sys.Date())
  }
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
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"linlog.plots")
)

save_linlog_plots(
  target_multicast_percent_plot,
  swgbf.prospective.component.target.multicasts,
  swg.prospective.voxel.data,
  t.target.specs, m.forecast.types,
  file.path(epiproject.cache.dir,"linlog.plots")
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
