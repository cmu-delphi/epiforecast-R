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

source("covidplot.R")

source("covid19ilinet-202003-202008-natreg-config.R")
source("covid19ilinet-templates.R")

## source("generate-retro-and-prospective-forecasts.R")



## source("gen-retro-component-forecasts.R")
## save_spreadsheets(
##     swgbf.retro.component.target.multicasts,
##     swg.retro.voxel.data,
##     t.target.specs, m.forecast.types,
##     epigroup.colname,
##     file.path(epiproject.cache.dir,"test001-component-spreadsheets"),
##     format_spreadsheet = function(spreadsheet) {
##         reformat_to_predx_v2_spreadsheet(spreadsheet, covid19ilinet.natreg.spreadsheet.template)
##     },
##     post_action = function(spreadsheet,filepath,spreadsheet.dir,subpath,swg.voxel.data,s,w,...) {
##         plot.parent.dir = file.path(epiproject.cache.dir, "test001-component-plots")
##         if (!dir.exists(plot.parent.dir)) {
##             dir.create(plot.parent.dir)
##         }
##         output.dir = file.path(plot.parent.dir, gsub("\\.csv$","",subpath))
##         if (!dir.exists(output.dir)) {
##             dir.create(output.dir)
##         }
##         plot.covid.forecast(filepath, output.dir)
##     }
## )



source("gen-prospective-component-forecasts.R")
save_spreadsheets(
    swgbf.prospective.component.target.multicasts,
    swg.prospective.voxel.data,
    t.target.specs, m.forecast.types,
    epigroup.colname,
    file.path(epiproject.cache.dir,"test001-component-spreadsheets"),
    format_spreadsheet = function(spreadsheet) {
        reformat_to_predx_v2_spreadsheet(spreadsheet, covid19ilinet.natreg.spreadsheet.template)
    },
    post_action = function(spreadsheet,filepath,spreadsheet.dir,subpath,swg.voxel.data,s,w,...) {
        plot.parent.dir = file.path(epiproject.cache.dir, "test001-component-plots")
        if (!dir.exists(plot.parent.dir)) {
            dir.create(plot.parent.dir)
        }
        output.dir = file.path(plot.parent.dir, gsub("\\.csv$","",subpath))
        if (!dir.exists(output.dir)) {
            dir.create(output.dir)
        }
        plot.covid.forecast(filepath, output.dir)
    }
)



## todo sin-cos-intercept-holidayindicatorinteraction pancaster
## todo use different version of nowcast
## todo other locations
## todo spreadsheet format
## todo plotting and sanity checks
## todo exclude methods that don't make sense
