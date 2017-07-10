## author_header begin
## Copyright (C) 2016 Sangwon Hyun
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

library(epiforecast)

## Make dummy csv
area.name = "hhs1"
hhs1.dat = fetchEpidataFullDat("fluview", area.name, "wili",
                               min.points.in.season=52L,
                               first.week.of.season = 21L,
                               cache.file.prefix=sprintf("fluview_%s_fetch", area.name))

## Create a csv in the correct, desirable format
filename="./correct.csv"
hhs1.dat = lapply(hhs1.dat, function(myvec){ if(length(myvec)<53) return(c(myvec,NA))  else {return(myvec)}})
correct.hhs1.dat = do.call(cbind, hhs1.dat)
correct.hhs1.dat = correct.hhs1.dat[-53,]
write.csv(correct.hhs1.dat, file = filename,row.names=FALSE)
full.dat = read.from.file(filename)


## Try EB with basic options
mysim = eb.sim(full.dat, n.sims=1000)
plot(mysim)
targets = forecast(mysim,'pht')


## Try EB with more simulation settings
control.list = get_eb_control_list(sd.option="prior",max.match.length=5)
mysim = eb.sim(full.dat, n.sims=100, control.list=control.list)
## epiforecast:::plot.sim(mysim)
plot(mysim)
targets = forecast(mysim,'pwk')


## Try BR
mysim = br.sim(full.dat, n.sims=100, bootstrap=T)
plot(mysim)
print(mysim,verbose=TRUE)
targets = forecast(mysim, "pht")

## Try BR with more simulation settings
control.list = get_br_control_list(df=5, cv.rule="1se")
mysim = br.sim(full.dat, n.sims=100, bootstrap=T, control.list=control.list)
plot(mysim)
print(mysim,verbose=TRUE)
targets = forecast(mysim, "pht", plot.hist=TRUE)
targets = forecast(mysim, "pwk", plot.hist=TRUE)

## Try twkde
mysim = twkde.sim(full.dat)
plot(mysim)
