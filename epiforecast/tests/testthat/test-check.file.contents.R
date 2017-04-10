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

context("Testing the check.file.contents() function..")

## Setup
dummytable = matrix(c(runif(n=9,1,3),NA),ncol=5)
dummytable.full = matrix(runif(n=10,1,3),ncol=5)
colnames(dummytable.full) = c(1997:2001)
library(s

## Tests:

test_that("Missing column names returns error.", {
    wrong.names = c("a","b","c","","d")
    colnames(dummytable) = wrong.names
    write.table(dummytable,
                file = "a.csv",
                col.names = T,
                row.names=F,
                sep=",")
    expect_error(check.file.contents("a.csv"))
    file.remove("a.csv")
})


test_that("Column names containing NA returns error.", {
    wrong.names = c(NA,1:4)
    colnames(dummytable) = wrong.names
    write.table(dummytable,
                file = "a.csv",
                col.names = T,
                row.names= F,
                sep=",")
    expect_error(check.file.contents("a.csv"))
    file.remove("a.csv")
})



test_that("Last column being full returns error.", {
    write.table(dummytable.full,
                file = "a.csv",
                col.names = T,
                row.names=F,
                sep=",")
    print(fread("a.csv",header=TRUE))
    expect_error(check.file.contents("a.csv"))
    file.remove("a.csv")
})
