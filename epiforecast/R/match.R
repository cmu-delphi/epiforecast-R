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

##' Match integer-valued input
##'
##' Returns a possibly-named integer-class vector version of the input, or
##' produces an error if the input seems inappropriate.
##'
##' @param inp supposed to be a possibly-named numeric object with integer/NA
##'   values.
##'
##' @return \code{inp} as an integer vector
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.integer = function(inp) {
  if (is.numeric(inp) && all(is.na(inp) | as.integer(inp)==inp)) {
    storage.mode(inp) <- "integer" # like as.integer, but doesn't strip attrs like names
    return (inp)
  } else {
    stop (sprintf("Argument =%s= must be a possibly-named numeric vector containing integral/NA values.",
                  paste(deparse(match.call()$inp), collapse="\n")))
  }
}

##' Match length-1 numeric input
##'
##' Returns a possibly-named length-1 possibly-NA numeric vector version of the
##' input, or produces an error if the input seems inappropriate.
##'
##' @param x supposed to be a possibly-named length-1 numeric object
##'
##' @return \code{x} as a length-1 numeric vector
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.single.na.or.numeric = function(x) {
  if (is.numeric(x) && length(x)==1L) {
    return (as.vector(x))
  } else {
    stop (sprintf("Argument =%s= must be a length 1 numeric vector representing an integral value.", paste(deparse(match.call()$n), collapse="\n")))
  }
}

##' Match length-1 non-\code{NA} numeric-valued input
##'
##' Returns a possibly-named length-1 non-NA numeric-class vector version of the
##' input, or produces an error if the input seems inappropriate.
##'
##' @param x supposed to be a possibly-named length-1 non-NA numeric-valued
##'   numeric object
##'
##' @return \code{x} as a length-1 numeric-class vector
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.single.nonna.numeric = function(x) {
  if (is.numeric(x) && length(x)==1L && !is.na(x)) {
    x <- as.vector(x)
    storage.mode(x) <- "numeric"
    return (x)
  } else {
    stop (sprintf("Argument =%s= must be a length-1 numeric object.", paste(deparse(match.call()$x), collapse="\n")))
  }
}

##' Match length-1 non-\code{NA} integer-valued input
##'
##' Returns a possibly-named length-1 non-NA integer-class vector version of the
##' input, or produces an error if the input seems inappropriate.
##'
##' @param n supposed to be a possibly-named length-1 non-NA integer-valued
##'   is.numeric object
##'
##' @return \code{n} as a length-1 integer-class vector
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.single.nonna.integer = function(n) {
  if (is.numeric(n) && length(n)==1L && !is.na(n) && as.integer(n)==n) {
    n <- as.vector(n)
    storage.mode(n) <- "integer"
    return (n)
  } else {
    stop (sprintf("Argument =%s= must be a length-1 integer-valued numeric object.", paste(deparse(match.call()$n), collapse="\n")))
  }
}

##' Match length-1 non-\code{NA} integer-valued input or \code{NULL}
##'
##' Returns a possibly-named length-1 non-NA integer-class vector version of the
##' input if non-\code{NULL}, \code{NULL} if input is \code{NULL}, or produces
##' an error if the input seems inappropriate.
##'
##' @param n supposed to be a possibly-named length-1 non-\code{NA}
##'   integer-valued numeric object or \code{NULL}
##'
##' @return \code{n} as a length-1 integer-class vector
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.single.nonna.integer.or.null = function(n) {
  if (is.null(n)) {
    return (NULL)
  } else if (is.numeric(n) && length(n)==1L && !is.na(n) && as.integer(n)==n) {
    n <- as.vector(n)
    storage.mode(n) <- "integer"
    return (n)
  } else {
    stop (sprintf("Argument =%s= must be a length-1 integer-valued numeric object.", paste(deparse(match.call()$n), collapse="\n")))
  }
}

##' Match all non-negative --- within \code{eps} --- \code{is.numeric} input
##'
##' Returns the input with any nearly-non-negative entries replaced with 0 if
##' the original input \code{is.numeric} and all of its entries are within
##' \code{eps} of being non-negative. Otherwise raises an error.
##'
##' @param x object that \code{is.numeric}; object to check
##'
##' @param eps single non-NA non-negative \code{is.numeric}; threshold of tolerance for negative values
##'
##' @return \code{x} with negative (nearly non-negative) entries replaced with 0
##'   if \code{x} meets the nearly non-negative criterion; otherwise an error is
##'   raised
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.nonnegative.numeric = function(x, eps=sqrt(.Machine[["double.eps"]])) {
  eps <- match.single.nonna.numeric(eps)
  if (eps < 0) {
    stop ("eps must be >= 0")
  }
  if (is.numeric(x) && all(x >= -eps)) {
    return (pmax(x, 0))
  } else {
    stop (sprintf("All entries in argument =%s= should be non-negative.", paste(deparse(match.call()$n), collapse="\n")))
  }
}

##' Match length-1 non-\code{NA} numeric-valued input
##'
##' Returns a possibly-named length-1 non-NA numeric-class vector version of the
##' input, or produces an error if the input seems inappropriate.
##'
##' @param x supposed to be a possibly-named length-1 non-NA numeric-valued
##'   numeric object
##'
##' @return \code{x} as a length-1 numeric-class vector
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.single.nonna.numeric = function(x) {
  if (is.numeric(x) && length(x)==1L && !is.na(x)) {
    x <- as.vector(x)
    storage.mode(x) <- "numeric"
    return (x)
  } else {
    stop (sprintf("Argument =%s= must be a length-1 numeric object.", paste(deparse(match.call()$x), collapse="\n")))
  }
}

##' Match dat object input
##'
##' Returns a list of possibly-named numeric-class vectors given a list of
##' possibly-named (is.)numeric vectors as input, or generates an error if the
##' input seems inappropriate.
##'
##' @param dat supposed to be a list of possibly-named numeric vectors
##'
##' @return \code{dat} as a list of numeric-class vectors
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.dat = function(dat) {
  if (is.list(dat)
      && all(sapply(dat, is.numeric))
      && all(sapply(dat, is.vector))) {
    return (lapply(dat, as.numeric))
  } else {
    stop (sprintf("Argument =%s= must be a list of numeric vectors.", paste(deparse(match.call()$dat), collapse="\n")))
  }
}

##' Convert trajectory, trajectory matrix, or sim to sim; error otherwise
##'
##' *.sim methods should eventually all support sim objects as input for new.dat
##' rather than just single trajectories, but, for the convenience of the user,
##' allow other types of input as well. Specifically, we should accept:
##' \itemize{
##'   \item{Trajectory: }{a numeric vector;}
##'   \item{Trajectory matrix: }{a numeric matrix with each column a trajectory;
##'   and}
##'   \item{Sim: }{a list with \code{$ys} a #times by #trajectories numeric
##'   matrix with #trajectories >= 1 and each row either all \code{NA} or all
##'   non-\code{NA}, and \code{$weights} a #trajectories-length numeric matrix
##'   with entries all >= 0.}
##' }
##'
##' This method checks that it receives such an input and outputs a
##' corresponding sim object.
##'
##' @param new.dat.sim trajectory / trajectory matrix / sim object
##'
##' @return sim object
##'
##' @examples
##' match.new.dat.sim(1:5)
##' match.new.dat.sim(as.matrix(1:5)[,rep(1,10)])
##' match.new.dat.sim(list(ys=as.matrix(1:5)[,rep(1,10)], weights=0:9))
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
match.new.dat.sim = function(new.dat.sim) {
  if (inherits(new.dat.sim, "sim")) {
    return (new.dat.sim)
  } else if (is.numeric(new.dat.sim) && is.vector(new.dat.sim)) {
    ys = as.matrix(new.dat.sim)
    weights = 1
    sim = structure(
      list(ys=ys, weights=weights,
           control.list=list(
             model="single trajectory as sim"
           )),
      class="sim")
    return (sim)
  } else if (is.numeric(new.dat.sim) && is.matrix(new.dat.sim)) {
    ys = new.dat.sim
    weights = rep(1/ncol(ys), ncol(ys)) # normalization isn't necessary though
    sim = structure(
      list(ys=ys, weights=weights,
           control.list=list(
             model="trajectory matrix as sim"
           )),
      class="sim")
    return (sim)
  } else if (is.list(new.dat.sim)
             && is.numeric(new.dat.sim$ys) && is.matrix(new.dat.sim$ys)
             && ncol(new.dat.sim$ys) >= 1
             && all(is.na(new.dat.sim$ys) == is.na(new.dat.sim$ys[,1]))
             && is.numeric(new.dat.sim$weights) && is.vector(new.dat.sim$weights)
             && all(new.dat.sim$weights >= 0)
             && ncol(new.dat.sim$ys) == length(new.dat.sim$weights)) {
    ## xxx consider requiring sum(new.dat.sim$weights) > 0
    if (sum(new.dat.sim$weights) < .Machine$double.eps^0.5) {
      warning("sum(new.dat.sim$weights) < .Machine$double.eps^0.5")
    }
    sim = new.dat.sim
    sim = structure(
      c(new.dat.sim,
        list(control.list=list(
               model="list with ys & weights as sim"
             ))
        ),
      class="sim")
    return (sim)
  } else {
    stop (sprintf("Argument =%s= must be (a) a numeric vector (a single trajectory), (b) a numeric matrix (cbound trajectories), or (c) a list with $ys a #times by #trajectories matrix with #trajectories >= 1 and each row either all NA or non-NA, and $weights a #trajectories-length numeric vector with entries all >=0.", paste(deparse(match.call()$new.dat.sim), collapse="\n")))
  }
}

##' \code{match.arg} variant replacing unmatched args with \code{choices[[1]]}, allowing non-character choices
##'
##' Assumes this usage: parent_fun = function(parent_arg[=parent_choices]) {
##' ... match.arg.forgiving(parent_arg) ... } (with or without
##' "=parent_choices").
##'
##' If \code{arg} is \code{NULL}, returns \code{choices[[1]]}.
##'
##' If \code{choices} is a character vector, this performs partial matches;
##' otherwise, it checks for \code{arg}'s that are \code{all.equal} with
##' \code{check.attributes=FALSE}.
##' 
##' @param arg the argument to match to a choice; should
##'
##' @param choices a positive-length vector; if it contains \code{NULL}, first
##'   choice should should be \code{NULL} to avoid ambiguity
##'
##' @return \code{arg}, the corresponding match in \code{choices}, or
##'   \code{choices[[1]]} with a warning (when \code{arg} fails to match a
##'   choice)
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @examples
##' library(testthat)
##' @example /tests/testthat/test_match.arg.or.default.R
##' 
##' @export
match.arg.else.default = function(arg, choices) {
  ## Assume this usage: parent_fun = function(parent_arg[=parent_choices]) {
  ## ... match.arg.else.default(parent_arg) ... } (with or without
  ## "=parent_choices").
  parent_frame = sys.frame(sys.parent())
  parent_function = sys.function(sys.parent())
  parent_match_call = match.call(parent_function, sys.call(sys.parent()))
  this_match_call = match.call() 
  parent_arg_name = as.character(this_match_call$arg)
  parent_choices_expr = formals(parent_function)[[parent_arg_name]]
  if (missing(choices)) {
    ## Assume parent_choices was specified; fill in choices with parent_choices:
    parent_choices = eval(parent_choices_expr, parent_frame)
    choices <- parent_choices
  }
  if (!is.vector(choices) || length(choices)==0L) {
    stop ("Problem with specification of choices (not arg): choices must be a positive-length vector.")
  } else if (any(is.null(choices)) && !is.null(choices[[1L]])) {
    stop ("Problem with specification of choices (not arg): if choices contains NULL, the first choice must be NULL, so that passing NULL as an arg is not ambiguous.")
  } else if (is.null(parent_match_call[[parent_arg_name]])) {
    ## parent_arg was missing in parent call; return first choice:
    return (choices[[1L]])
  } else if (is.null(arg)) {
    ## parent_arg was NULL; return first choice:
    return (choices[[1L]])
  } else {
    if (is.character(choices)) {
      if (!is.character(arg) || length(arg) != 1L) {
        stop (sprintf('%s must be NULL or a length-1 character with a partial match in %s.',
                      parent_arg_name, deparse(parent_choices_expr)))
      }
      ind = pmatch(arg, choices, nomatch=0L, duplicates.ok=TRUE) # last two for speed xxx check this
      ## fixme if there is a nonunique partial match, an error should be
      ## generated, not a warning&default
      if (ind==0L) {
        ## arg was not a partial match to a choice; return first choice with warning
        warning (sprintf('%s should have a partial match in %s; defaulting to %s.',
                         parent_arg_name, deparse(parent_choices_expr), choices[[1L]]))
        return (choices[[1L]])
      } else {
        ## arg was (partial) match to a choice; return that choice
        return (choices[[ind]])
      }
    } else { ## choices is not a character vector
      ## ind = Position(function(choice) idential(arg, choice), choices, nomatch=0L)
      ind = Position(function(choice) isTRUE(all.equal(arg, choice, check.attributes=FALSE)), choices, nomatch=0L) # be a bit forgiving
      if (ind == 0L) {
        ## arg was not in choices, return first choice with warning
        warning (sprintf('%s should be all.equal to something in choice vector %s; defaulting to first choice.',
                         parent_arg_name, deparse(parent_choices_expr)))
        return (choices[[1L]])
      } else {
        ## arg was all.equal to some choice: return that choice
        return (choices[[ind]])
      }
    }
  }
}
## todo testing on match.arg.else.default
