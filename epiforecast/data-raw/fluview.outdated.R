## author_header begin
## Copyright (C) 2016 Logan C. Brooks
##
## This file is part of epiforecast.  Algorithms included in epiforecast were developed by Logan C. Brooks, David C. Farrow, Sangwon Hyun, Shannon Gallagher, Ryan J. Tibshirani, Roni Rosenfeld, and Rob Tibshirani (Stanford University), members of the Delphi group at Carnegie Mellon University.
##
## Research reported in this publication was supported by the National Institute Of General Medical Sciences of the National Institutes of Health under Award Number U54 GM088491. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health. This material is based upon work supported by the National Science Foundation Graduate Research Fellowship Program under Grant No. DGE-1252522. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation. David C. Farrow was a predoctoral trainee supported by NIH T32 training grant T32 EB009403 as part of the HHMI-NIBIB Interfaces Initiative. Ryan J. Tibshirani was supported by NSF grant DMS-1309174.
## author_header end
## license_header begin
## epiforecast is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## epiforecast is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with epiforecast.  If not, see <http://www.gnu.org/licenses/>.
## license_header end

## source("R/utils.R", chdir=TRUE)
## source("R/match.R", chdir=TRUE)
## source("R/weeks.R", chdir=TRUE)
## source("R/delphi_epidata.R", chdir=TRUE)
## source("R/loaders.R", chdir=TRUE)

## devtools::load_code()

devtools::load_all()

fluview.first.model.week = 21L
fluview.area.names = c(sprintf("hhs%d",1:10),"nat")
## fluview.2003on.dfs = structure(lapply(fluview.area.names, function(area.name) {
##   trimPartialPastSeasons(fetchEpidataDF("fluview", area.name,
##                                         first.week.of.season=fluview.first.model.week,
##                                         cache.file.prefix=sprintf("fluview_%s_fetch.Rdata", area.name)),
##                          "wili", 52)
## }), names=fluview.area.names)
## fluview.2003on.full.dats = lapply(fluview.2003on.dfs, function(df) {
##   full.dat = split(df$wili, df$season) # historical seasons + current season
##   names(full.dat) <- sprintf("S%s", names(full.dat))
##   full.dat <- full.dat[names(full.dat)!="S2009"]
##   full.dat
## })
fluview.2003on.full.dats = setNames(lapply(fluview.area.names, function(area.name) {
  fetchEpidataFullDat("fluview", area.name, "wili",
                      min.points.in.season=52L, first.week.of.season=fluview.first.model.week,
                      cache.file.prefix=sprintf("fluview_%s_fetch.Rdata", area.name))
}), fluview.area.names)
fluview.2003on.dats = lapply(fluview.2003on.full.dats, head, n=-1L)
fluview.2003on.new.dats = lapply(fluview.2003on.full.dats, function(full.dat) tail(full.dat, 1L)[[1]])

## todo document

fluview.2003on.outdated.dats = fluview.2003on.dats
fluview.2003on.outdated.new.dats = fluview.2003on.new.dats

devtools::use_data(fluview.2003on.outdated.dats, fluview.2003on.outdated.new.dats,
                   overwrite = TRUE)

## fixme should this just be a method? get.fluview.2003on.dats()? or something
## formatted like fetchEpidataDF but returning a dat and new.dat, or a full.dat?
