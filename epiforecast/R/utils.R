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


##' A \code{seq} variant that produces a 0-length vector when \code{!(from <=
##' to)}.
##'
##' @param from starting number (or other object compatible with \code{seq})
##' @param to ending number (or other object compatible with \code{seq})
##' @param ... arguments to forward to \code{seq}
##'
##' @export
Seq = function(from, to, ...) {
  if (from <= to) {
    ## Just delegate in the "nice" case:
    return (seq(from, to, ...))
  } else {
    ## Return 0 elements of a reasonable type by taking 0 elements from
    ## =seq(from,from)=.
    return (seq(from, from)[integer(0)])
    ## We could do =seq(from, to, ...)=, but that could potentially hurt
    ## performance and generate bugs for certain dispatches.
  }
}

dat.to.matrix1 = function(dat, n) {
  dat <- match.dat(dat)
  n <- match.single.nonna.integer(n)
  sapply(dat, `[`, seq_len(n))
}
dat.to.matrix2 = function(dat, n) {
  dat <- match.dat(dat)
  n <- match.single.nonna.integer(n)
  result = matrix(NA_real_, n, length(dat))
  for (trajectory.i in seq_along(dat)) {
    trajectory = dat[[trajectory.i]]
    result[,trajectory.i] <- trajectory[seq_len(n)]
  }
  dimnames(result)[[2]] <- names(dat)
  return (result)
}
dat.to.matrix3 = function(dat, n) {
  dat <- match.dat(dat)
  n <- match.single.nonna.integer(n)
  result = matrix(NA_real_, n, length(dat))
  inds = seq_len(n)
  for (trajectory.i in seq_along(dat)) {
    trajectory = dat[[trajectory.i]]
    result[,trajectory.i] <- trajectory[inds]
  }
  dimnames(result)[[2]] <- names(dat)
  return (result)
}
## all.equal(dat.to.matrix1(dat, 52), dat.to.matrix2(dat, 52), dat.to.matrix3(dat, 52))
## benchmark(dat.to.matrix1(dat,52),dat.to.matrix2(dat,52),dat.to.matrix3(dat,52),replications=100000)
## 1 dat.to.matrix1(dat, 52)       100000  35.548    4.814    35.564    0.004
## 2 dat.to.matrix2(dat, 52)       100000   7.948    1.076     7.952    0.000
## 3 dat.to.matrix3(dat, 52)       100000   7.384    1.000     7.388    0.000

##' Numeric matrix of the first \code{n} elements of each numeric vector in
##' \code{dat}.
##'
##' A more efficient implementation of \code{sapply(dat, `[`, seq_len(n))}. Any
##' vectors in \code{dat} with length less than \code{n} are extended with
##' \code{NA_real_}'s at the end.
##'
##' @param dat a list of numeric vectors
##' @param n a single integer: the number of elements to take from each vector
##' @return a \code{n}-by-\code{length(dat)} numeric matrix
##'
##' @examples
##' dat = list(11:15, 21:26)
##' dat.to.matrix(dat, 5) # (5x2: dat[[2]] is cut off)
##' dat.to.matrix(dat, 6) # (6x2: dat[[1]] is extended with NA_real_)
##' n = 3
##' identical(c(n, length(dat)), dim(dat.to.matrix(dat, n)))
##'
##' @export
dat.to.matrix = dat.to.matrix3
