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

##' Test if Date is Christmas Day (vectorized)
##'
##' @param Date \code{Date} vector to test
##' @return logical vector
##' 
##' @export
is_christmas = function(Date) {
  posixlt = as.POSIXlt(Date)
  return (posixlt$mon + 1L == 12L & posixlt$mday == 25L)
}

##' Test if Date is (Gregorian) New Year's Day (vectorized)
##'
##' @param Date \code{Date} vector to test
##' @return logical vector
##'
##' @export
is_newyear = function(Date) {
  posixlt = as.POSIXlt(Date)
  return (posixlt$mon + 1L == 1L & posixlt$mday == 1L)
}

##' Test if Date is Thanksgiving Day (vectorized)
##'
##' @param Date \code{Date} vector to test
##' @return logical vector
##'
##' @export
is_thanksgiving = function(Date) {
  posixlt = as.POSIXlt(Date)
  return (posixlt$mon + 1L == 11L & (posixlt$mday - 1L) %/% 7L + 1L == 4L & posixlt$wday == 4L)
}

## todo: Chinese New Year; cannot use seasonal::cny since it is GPL3.0;
## investigate seasonal-cited source,
## http://www.chinesenewyears.info/chinese-new-year-calendar.php (license does
## not mention reuse) or ConvCalendar package (projectpluto.com calendars)

## alternative: timeDate::holiday / chron::is.holiday
