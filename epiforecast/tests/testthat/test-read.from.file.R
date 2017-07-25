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

context("Testing the read.from.file() function..")

## ## Setup
## dummytable = matrix(c(runif(n=9,1,3),NA),ncol=5)
## colnames(dummytable) = c(1997:2001)


## ## Tests:
## test_that("Return type is list of vector", {
##     write.table(dummytable,
##                 file = "a.csv",
##                 col.names = T,
##                 row.names=F,
##                 sep=",")
##     mydat  = read.from.file("a.csv")
##     mydat
##     expect_equal("list", typeof(mydat))
##     file.remove("a.csv")
## })


## Fetch some data
area.name = "hhs1"
hhs1.dat = fetchEpidataFullDat("fluview", area.name, "wili",
                               min.points.in.season=52L,
                               first.week.of.season = 21L,
                               cache.file=sprintf("fluview_%s_fetch.Rdata", area.name))
alt.names1=paste("season",1:length(names(hhs1.dat)))
alt.names2 = alt.names1;    alt.names2[5] = "" 
alt.names3 = rep("", length(names(hhs1.dat)))

## When one label is missing
    

## Create csv with /no/ names
test_that("When alls labels are missing, error is thrown", {
    filename="./all.names.missing.file.csv"
    all.names.missing.hhs1.dat = hhs1.dat
    names(all.names.missing.hhs1.dat) = alt.names3
    all.names.missing.hhs1.dat = do.call(cbind, all.names.missing.hhs1.dat)
    all.names.missing.hhs1.dat = all.names.missing.hhs1.dat[-53,]
    write.csv(all.names.missing.hhs1.dat, file = filename, row.names=FALSE)
    expect_error(full.dat = read.from.file(filename))
})

## Create csv with one name missing
test_that("When one label is missing, error is thrown", {
    filename="./one.name.missing.file.csv"
    one.name.missing.hhs1.dat = hhs1.dat
    names(one.name.missing.hhs1.dat) = alt.names2
    one.name.missing.hhs1.dat = do.call(cbind, one.name.missing.hhs1.dat)
    one.name.missing.hhs1.dat = one.name.missing.hhs1.dat[-53,]
    write.csv(one.name.missing.hhs1.dat, file = filename, row.names=FALSE)
    expect_error(full.dat = read.from.file(filename))
})

## Last season is not partially observed
test_that("When last season is not partially observed, error is thrown", {
    filename="./last.column.weird.file.csv"
    last.column.weird.hhs1.dat = do.call(cbind, hhs1.dat)
    write.csv(last.column.weird.hhs1.dat, file = filename,row.names=FALSE)
    expect_error(full.dat = read.from.file(filename))
})

## CSV is correct, so everything should be fine here.
test_that("When last season is not partially observed, error is thrown", {
    filename="./correct.csv"
    correct.hhs1.dat = do.call(cbind, hhs1.dat)
    correct.hhs1.dat = correct.hhs1.dat[-53,]
    write.csv(correct.hhs1.dat, file = filename,row.names=FALSE)
    full.dat = read.from.file(filename)
    expect_equal(class(full.dat), "list")
})
