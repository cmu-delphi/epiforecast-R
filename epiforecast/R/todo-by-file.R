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

######################
## interface.R #######
######################

## todo length checks --- time.of.forecast not only with dat, but with fit
## todo remove n.out arg --- determine by min fit length

## ## source("../sample/sample_config_flu_1516.R")
## source("loaders.R")
## source("fitters.R")
## source("nidss/nidss_fetch_data.R", chdir=TRUE)

## nidss.datfit = fetchNIDSSDatFit("flu", "nationwide")

## olddatplus = nidss.datfit$olddatplus
## oldfit = nidss.datfit$fit
## include.ss = nidss.datfit$include.ss
## current.s = nidss.datfit$current.s
## first.year = nidss.datfit$first.year
## first.model.week = nidss.datfit$first.model.week

## fit.ss = include.ss[include.ss < current.s]
## exclude.2009.pandemic.season=TRUE
## if (exclude.2009.pandemic.season) {
##     fit.ss <- fit.ss[-2]
##     oldfit <- list(f=oldfit$f[,-2], tau=oldfit$tau[-2])
## }
## train.ss = fit.ss
## test.s = max(train.ss)+1

## newdat = olddatplus.to.newdat(olddatplus)
## newdat.attributes = attributes(newdat)
## newdat <- newdat[match(fit.ss, include.ss)]
## newdat.attributes$names <- names(newdat)
## attributes(newdat) <- newdat.attributes

## qwer = fit.eb.control.list(oldfit.to.newfit(oldfit), get.eb.control.list())
## asdf = eb.createForecasts(newdat, olddatplus$wili[olddatplus$season==test.s], oldfit.to.newfit(oldfit), 0)
## asdf = eb.createForecasts(newdat, olddatplus$wili[olddatplus$season==test.s], oldfit.to.newfit(oldfit), 1L)
## source("plotters.R")
## newfit = smooth.curves.to.newfit(eb.fitSmoothCurves(newdat))
## matplot.newdat(newdat)
## matplot.newfit(newdat, newfit)
## seriesplot.newfit(newdat, smooth.curves.to.newfit(eb.fitSmoothCurves(newdat)))

## xxx instead of n.out, allow NA's in the future trajectories, just fill in all; use !is.na as another ii.match mask?
## todo explicitly make object that represents a distribution of curves, corresponding fitting functions, then the conditioning method?
## todo rename forecast time to something with "ind"?
## todo documentation
## todo imports
## todo examples

######################
## loaders.R #########
######################

## xxx consider fetching by issue instead in fetchEpidataHistoryDF; at least for
## fluview, the set of all issues should be a subset of the set of all epiweeks
## from the current data frame

## todo version of mimicPastEpidataDF that never uses future data, instead
## taking the seasonally-expected change from the last available data point or
## stopping if there is no available data point beforehand (will need to handle finalized versions inputted later with lags outside the =lags= range... override their lag with the max lag in =lag= and update =issue= accordingly?)

## ## todo turn into test
## history.dt = fetchEpidataHistoryDT("fluview", "hhs1", 0:51,
##                            first.week.of.season = 31L,
##                            cache.file.prefix="~/.epiforecast-cache/fluview_hhs1")
## list(mimicPastEpidataDF1, mimicPastEpidataDF2) %>>%
##   lapply(function(mimicPastEpidataDFn) {
##     ## mimicPastEpidataDFn(history.dt, 201540L) %>>%
##     mimicPastEpidataDFn(history.dt, 201040L) %>>%
##       dplyr::arrange(-epiweek) %>>%
##       dplyr::select(epiweek, issue, forecast.epiweek, wili)
##   }) %>>%
##   do.call(what=identical) %>>%
##   {.}

## todo forecast.sim: rather than randomly sampling a value from a multival, use all and weight?
