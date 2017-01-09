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

parent_function = function(ch1=letters[1:5], ch2=c("AAA","AAB","BBB"),
                           num1=1:5, int1=6L:10L,
                           list1=list(1:2,3:4), list2=list(NULL,"a",c(two=2),1:5)) {
  return (list(
    ch1=match.arg.else.default(ch1),
    ch2=match.arg.else.default(ch2),
    num1=match.arg.else.default(num1),
    int1=match.arg.else.default(int1),
    list1=match.arg.else.default(list1),
    list2=match.arg.else.default(list2)
  ))
}

## Return default on missing:
expect_equal(parent_function(),
             list(ch1="a", ch2="AAA", num1=1, int1=6L, list1=1:2, list2=NULL))

## Return default on NULL:
expect_equal(parent_function(NULL, NULL, NULL, NULL, NULL, NULL),
             list(ch1="a", ch2="AAA", num1=1, int1=6L, list1=1:2, list2=NULL))

## Allow partial matches, all.equal ignoring attributes:
expect_equal(parent_function("b", c(extraneous.name="B"), 3L, 8.00000000001, c(p=3,q=4), 2),
             list(ch1="b", ch2="BBB", num1=3.0, int1=8L, list1=3:4, list2=c(two=2)))

## Return default with warning on mismatched inputs:
expect_equal(suppressWarnings(parent_function("q", "A", "nonnumeric", 11L, 1:4, c("A","B","C"))),
             list(ch1="a", ch2="AAA", num1=1, int1=6L, list1=1:2, list2=NULL))
expect_warning(parent_function("q"))
expect_warning(parent_function(,"A"))
expect_warning(parent_function(,,"nonnumeric"))
expect_warning(parent_function(,,,11L))
expect_warning(parent_function(,,,,1:4))
expect_warning(parent_function(,,,,,c("A","B","C")))

## Produce error on inappropriate inputs:
expect_error(parent_function(3), "length-1 character")
expect_error(parent_function(letters[1:2]), "length-1 character")

## todo produce error on inappropriate =choices=
