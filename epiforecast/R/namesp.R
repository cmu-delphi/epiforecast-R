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

##' names, but always outputting character vector of natural length
##'
##' If object does not have names, outputs a vector of empty strings of matching
##' length instead of NULL.
##'
##' @param x the object
##' @template param_invent.scalars
##'
##' @export
namesp = function(x) {
  maybe.namesp.x = names(x)
  ## require non-NULL result:
  if (is.null(maybe.namesp.x)) {
    maybe.namesp.x <- rep("", length(x))
  }
  return (maybe.namesp.x)
}

##' dim, but always outputting integer vector of natural length
##'
##' @param x the object
##' @template param_invent.scalars
##'
##' @export
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

##' The number of dimensions of an object; length of dimp
##'
##' @param x the object
##' @template param_invent.scalars
##'
##' @export
ndimp = function(x, invent.scalars=TRUE) {
  return (length(dimp(x, invent.scalars=invent.scalars)))
}

##' dimnames, but always as a list of lists with natural lengths
##'
##' Operates like dimnames, but always produces a list of lists with the "natural" lengths.  Fixes up cases where the dimnames or any element of the dimnames would be NULL.  Any element of dimnames (i.e., vector of names associated with a particular dimension) that would be NULL is replaced with a vector of empty strings.
##'
##' @param x array or array-like object of which to take and format the dimnames
##' @template param_invent.scalars
##'
##' @export
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

##' Get integer or character index sets for each dimension of an object
##'
##' Favors character-type index sets, using  the entries of \code{dimnamesp} for
##' every dimension with names that are not all empty strings.
##'
##' @param x an array-like object
##' @template param_invent.scalars
##'
##' @examples
##' library("pipeR")
##' array(1:120, 2:5) %>>%
##'   with_dimnames(list(
##'     "First dimension" = NULL,
##'     "Second dimension" = rep("", 3),
##'     "Third dimension" = c("","","c","d"),
##'     "Fourth dimension" = letters[1:5]
##'   )) %>>%
##' dimnames_or_inds()
##'
##' @export
dimnames_or_inds = function(x, invent.scalars=TRUE) {
  dnp.x = dimnamesp(x, invent.scalars=invent.scalars)
  dimnamesp_to_dimindices(dnp.x)
}

##' Post-process dimnamesp, changing all-\code{""} dimensions to \code{seq_along}'s
##'
##' @param x object of which to get the dimnames names
##' @template param_invent.scalars
##'
##' @examples
##' ## The first two dimensions are re-index with integers; the last two are not:
##' library("pipeR")
##' array(1:120, 2:5) %>>%
##'   with_dimnames(list(
##'     "First dimension" = NULL,
##'     "Second dimension" = rep("", 3),
##'     "Third dimension" = c("","","c","d"),
##'     "Fourth dimension" = letters[1:5]
##'   )) %>>%
##'   dimnamesp() %>>%
##'   dimnamesp_to_dimindices()
##'
##' @export
dimnamesp_to_dimindices = function(dnp) {
  lapply(dnp, function(dimension.np) {
    if (all(dimension.np=="")) {
      seq_along(dimension.np)
    } else {
      dimension.np
    }
  })
}

##' Pipe-friendly, robust dimnames names getter, as a character vector of natural length
##'
##' A    pipe-friendly,     robustified,    more    consistent     version    of
##' \code{names(dimnames(x))}.  Accepts array-like  objects, including  vectors (although these never have dimnames names).
##' Always  returns a  character  vector  with length  equal  to  the number  of
##' dimensions identified for \code{x} (never \code{NULL}).
##'
##' @param x array-like object
##' @template param_invent.scalars
##'
##' @examples
##' library("pipeR")
##' array(1:24, 2:4) %>>%
##'     with_dimnamesnames(c("a","","b")) %>>%
##'     dimnamesnamesp()
##' array(1:24, 2:4) %>>%
##'   dimnamesnamesp()
##' array(1:24, 2:4) %>>%
##'     with_dimnames(list(letters[1:2], letters[1:3], NULL)) %>>%
##'     {names(dimnames(.))[[2L]] <- "Second dimension"; .} %>>%
##'     dimnamesnamesp()
##' dimnamesnamesp(1:5)
##' dimnamesnamesp(1)
##' dimnamesnamesp(1, invent.scalars=FALSE)
##'
##' @export
dimnamesnamesp = function(x, invent.scalars=TRUE) {
  return (names(dimnamesp(x, invent.scalars=invent.scalars)))
}

##' Standard-evaluation, no-NULL-input version of \code{vector_as_named_array}
##'
##' @export
vector_as_named_array_ = function(vector, dimension.name, entry.names) {
  result = as.array(vector)
  dimnames(result) <- stats::setNames(list(entry.names), dimension.name)
  return (result)
}

##' Convert vector to 1-D array with non-NULL names & named dimnames
##'
##' @param vector the vector to convert
##' @param dimension.name length-1 character vector or \code{NULL} (default); if
##'     the former,  then the  \code{names(dimnames(output))[[1]]} will  be this
##'     value; otherwise, \code{names(dimnames(output))[[1]]}  will be set using
##'     nonstandard  evaluation   to  grab  the  expression   provided  for  the
##'     \code{vector} argument and convert it to a string
##'  @param entry.names  character  vector (expected  to be  of  same length  as
##'     \code{vector}) or  \code{NULL}; if the  former, the names of  the output
##'     are set to this value; if the latter, the names of the output are set to
##'     a character vector of the same  length: \code{vector} itself if it is an
##'     unnamed character vector, or \code{namesp(vector)} otherwise.
##'
##' @export
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

##' Pipe-friendly, robust way of setting names(dimnames(object))
##'
##' @examples
##' library("pipeR")
##' array(1:24, 2:4) %>>%
##'     with_dimnames(list(c("a","b"),
##'                        c("c","d","e"),
##'                        c("f","g","h","i"))) %>>%
##'     with_dimnamesnames(c("First dimension", "Second dimension", "Third dimension"))
##' ## Edge case with null dimnames is handled:
##' array(1:24, 2:4) %>>%
##'     with_dimnamesnames(c("First dimension", "Second dimension", "Third dimension"))
##' matrix(1:6, 2,3) %>>%
##'     with_dimnamesnames(c("First dimension", "Second dimension"))
##' ## Vectors are converted to arrays by with_dimnamesnames; the dimnames names
##' ## are set properly whether or not the vector is already has names, but the
##' ## print function will not show the dimnames names for the result if the
##' ## vector was originally unnamed:
##' from.named.vector = c(a=1,b=2) %>>% with_dimnamesnames("Only dimension")
##' class(from.named.vector) # "array"
##' print(from.named.vector) # shows dimnames names
##' names(dimnames(from.named.vector)) # "Only dimension"
##' from.unnamed.vector = c(1,2) %>>% with_dimnamesnames("Only dimension")
##' class(from.unnamed.vector) # "array"
##' print(from.unnamed.vector) # doesn't show dimnames names
##' names(dimnames(from.unnamed.vector)) # "Only dimension"
##'
##' @export
with_dimnamesnames = function(arraylike, dimnamesnames) {
  arraylike <- as.array(arraylike)
  ## replace NULL dimnames (but not NULL entries in dimnames):
  if (is.null(dimnames(arraylike))) {
    dimnames(arraylike) <- rep(list(NULL), length(dim(arraylike)))
  }
  names(dimnames(arraylike)) <- dimnamesnames
  return (arraylike)
}

##' Pipe-friendly, robust way of setting dimnames(object)
##'
##' @examples
##' library("pipeR")
##' array(1:24, 2:4) %>>%
##'     with_dimnames(list(c("a","b"),
##'                        c("c","d","e"),
##'                        c("f","g","h","i")))
##' ## Two ways of setting named dimnames:
##' array(1:24, 2:4) %>>%
##'     with_dimnames(list("First dimension"=c("a","b"),
##'                        "Second dimension"=c("c","d","e"),
##'                        "Third dimension"=c("f","g","h","i")))
##' array(1:24, 2:4) %>>%
##'     with_dimnames(list(c("a","b"),
##'                        c("c","d","e"),
##'                        c("f","g","h","i"))) %>>%
##'     with_dimnamesnames(c("First dimension", "Second dimension", "Third dimension"))
##' ## Vectors are converted to arrays by with_dimnames:
##' (1:5) %>>% with_dimnames(list(letters[1:5]))
##' (1:5) %>>% with_dimnames(list("Only dimension" = letters[1:5]))
##'
##' @export
with_dimnames = function(arraylike, dimnames) {
  arraylike <- as.array(arraylike)
  dimnames(arraylike) <- dimnames
  return (arraylike)
}

##' Names of a 1-D array, as a character array with the same dimnames(names)
##'
##' Given a 1-D  array or other object  with non-NULL names, turn it  into a 1-D
##' character  array with  values  given  by the  input's  names,  and the  same
##' dimnames as the input (if any), including the names of the dimnames. Objects
##' are converted  to arrays first  for more consistent behavior  (e.g., vectors
##' can have names but not dimnames; converting to an array gives dimnames equal
##' to \code{list(names(vec))}).
##'
##' @examples
##' library("pipeR")
##' array(1:5, 5) %>>%
##'   with_dimnames(list(letters[1:5])) %>>%
##'   with_dimnamesnames("Only dimension") %>>%
##'   named_arrayvec_to_name_arrayvec()
##' array(1:5, 5) %>>%
##'     with_dimnames(list(letters[1:5])) %>>%
##'     named_arrayvec_to_name_arrayvec()
##' ## Conversion to arrays is attempted for non-array inputs, e.g., the two
##' ## expressions below should be exchangeable:
##' stats::setNames(1:5, letters[1:5]) %>>%
##'     as.array() %>>%
##'     named_arrayvec_to_name_arrayvec()
##' stats::setNames(1:5, letters[1:5]) %>>%
##'     named_arrayvec_to_name_arrayvec()
##'
##' @export
named_arrayvec_to_name_arrayvec =
    function(arrayvec) {
        arrayvec <- as.array(arrayvec)
        if (is.null(names(arrayvec))) {
            stop ("Input must have non-NULL names.")
        }
        with_dimnames(names(arrayvec), dimnames(arrayvec))
    }

##' dimnamesp of an array, as a list of named lists
##'
##' @examples
##' library("pipeR")
##' array(1:24, 2:4) %>>%
##'     with_dimnames(list(c("a","b"),
##'                        c("c","d","e"),
##'                        c("f","g","h","i"))) %>>%
##'     named_array_to_name_arrayvecs()
##' array(1:120, 2:5) %>>%
##'     with_dimnames(list("First dimension"=c("a","b"),
##'                        "Second dimension"=NULL,
##'                        c("f","g","h","i"),
##'                        NULL)) %>>%
##'     named_array_to_name_arrayvecs()
##' ## Some unusual inputs may give surprising output that may not have the exact
##' ## same properties as for more common cases shown above. Behavior on these cases
##' ## likely may change in subsequent versions of this package. In particular,
##' ## behavior for arrays with names assigned to elements of dimnames may be
##' ## surprising and should be considered unstable functionality:
##' array(1:24, 2:4) %>>%
##'     with_dimnames(list("Unusual dimension"=c(A="a",B="b"),
##'                        c("c","d","e"),
##'                        c("f","g","h","i"))) %>>%
##'     named_array_to_name_arrayvecs()
##'
##' @export
named_array_to_name_arrayvecs = function(named.array) {
  dimnamesp(named.array) %>>%
    Map(f=vector_as_named_array, names(.)) %>>%
    unname()
}
