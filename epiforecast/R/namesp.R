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

namesp = function(x) {
  maybe.namesp.x = names(x)
  ## require non-NULL result:
  if (is.null(maybe.namesp.x)) {
    maybe.namesp.x <- rep("", length(x))
  }
  return (maybe.namesp.x)
}

dimp = function(x, invent.scalars=TRUE) {
  maybe.dimp.x = dim(x)
  ## for dim-less things, treat as 1-D or scalar:
  if (is.null(maybe.dimp.x)) {
    if (invent.scalars && length(x)==1L && namesp(x)=="") {
      maybe.dimp.x <- integer(0L)
    } else {
      maybe.dimp.x <- length(x)
    }
  }
  return (maybe.dimp.x)
}

ndimp = function(x, invent.scalars=TRUE) {
  return (length(dimp(x, invent.scalars=invent.scalars)))
}

dimnamesp = function(x, invent.scalars=TRUE) {
  dimp.x = dimp(x, invent.scalars=invent.scalars)
  maybe.dimnamesp.x = dimnames(x)
  ## replace NULL dimnames with singleton list of 1-D names or list of NULL's:
  if (is.null(maybe.dimnamesp.x)) {
    if (length(dimp.x)==1L) {
      maybe.dimnamesp.x <- list(names(x))
    } else {
      maybe.dimnamesp.x <- rep(list(NULL), length(dimp.x))
    }
  }
  ## replace any NULL's within dimnames with vector of ""'s
  nameless.dim.is = which(vapply(maybe.dimnamesp.x, is.null, TRUE))
  maybe.dimnamesp.x[nameless.dim.is] <-
    lapply(nameless.dim.is, function(dim.i) {
      rep("", dimp.x[[dim.i]])
    })
  ## require non-NULL result names:
  if (is.null(names(maybe.dimnamesp.x))) {
    names(maybe.dimnamesp.x) <- rep("", length(maybe.dimnamesp.x))
  }
  return (maybe.dimnamesp.x)
}

dimnamesnamesp = function(x, invent.scalars=TRUE) {
  return (names(dimnamesp(x, invent.scalars=invent.scalars)))
}

vector_as_named_array_ = function(vector, dimension.name, entry.names) {
  result = as.array(vector)
  dimnames(result) <- stats::setNames(list(entry.names), dimension.name)
  return (result)
}

vector_as_named_array = function(vector,
                                 dimension.name=NULL,
                                 entry.names=NULL) {
  if (is.null(dimension.name)) {
    dimension.name <- capture.output(print(match.call()[["vector"]]))
  }
  if (is.null(entry.names)) {
    if (is.null(names(vector)) && is.character(vector)) {
      entry.names <- vector
    } else {
      entry.names <- namesp(vector)
    }
  }
  return (vector_as_named_array_(vector, dimension.name, entry.names))
}

with_dimnamesnames = function(arraylike, dimnamesnames) {
  arraylike <- as.array(arraylike)
  ## replace NULL dimnames (but not NULL entries in dimnames):
  if (is.null(dimnames(arraylike))) {
    dimnames(arraylike) <- rep(list(NULL), length(dim(arraylike)))
  }
  names(dimnames(arraylike)) <- dimnamesnames
  return (arraylike)
}

with_dimnames = function(arraylike, dimnames) {
  arraylike <- as.array(arraylike)
  dimnames(arraylike) <- dimnames
  return (arraylike)
}

named_arrayvec_to_name_arrayvec = function(arrayvec) {
  with_dimnames(names(arrayvec), dimnames(arrayvec))
}

named_array_to_name_arrayvecs = function(named.array) {
  dimnamesp(named.array) %>>%
    Map(f=vector_as_named_array, names(.)) %>>%
    unname()
}
