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

## Do nothing if \code{trueOrStopMessage} is \code{TRUE}; otherwise,
## output the error message given by \code{trueOrStopMessage}.  Used
## with \code{is.or.why.not.*} functions.
ensureTRUE = function(trueOrStopMessage) {
    if (!isTRUE(trueOrStopMessage))
        stop(trueOrStopMessage)
}

## Do nothing if \code{trueOrStopMessage} is \code{TRUE}; otherwise,
## output the error message given by \code{trueOrStopMessage} with a
## tag indicating that it is a post-condition check failure
## (indicating a bug).  Used with \code{is.or.why.not.*} functions.
ensureTRUEpostcondition = function(trueOrStopMessage) {
    if (!isTRUE(trueOrStopMessage))
        stop(paste("Post-condition check failed: ",trueOrStopMessage))
}

## Output \code{TRUE} if \code{smooth.dat} is list of numeric vectors
## with same length as \code{dat}, but each with length 53 and no
## \code{NA}'s (should be a smoothed version of \code{dat}).
## Otherwise, output a string to be used in an error message.
is.or.why.not.smooth.dat = function(dat, smooth.dat) {
    if (length(dat) != length(smooth.dat)) {
        return("length(dat) != length(smooth.dat)")
    }
    ## if (length(smooth.dat) != 0 &&
    ##     !identical(53L, unique(sapply(smooth.dat, length)))) {
    ##     return("=smooth.dat= should contain only 53-length vectors")
    ## }
    if (!identical(unname(lengths(dat)), unname(lengths(smooth.dat))))
        stop("lengths of vectors in =dat= and =smooth.dat= do not match")
    ## xxx additional checks for numeric type
    return (TRUE)
}
## xxx remove reference to 53, maybe change around length requirements to match =eb.fitSmoothCurves=.

## Output \code{TRUE} if \code{curve.models} is a list of lists; one
## list per season in \code{dat}, each list containing three elements
## \code{`f`}, the corresponding smoothed curve, \code{`tau`}, the
## estimate of the sd under the iid Gaussian noise assumption, and
## \code{`type`}, a string ("Gaussian") indicating what noise model
## was used.  Otherwise, output a string to be used in an error
## message.
is.or.why.not.curve.models = function(dat, curve.models) {
    if (length(dat) != length(curve.models))
        return ("length(dat) != length(curve.models)")
    if (!inherits(curve.models, "list") ||
          !all(vapply(curve.models, FUN.VALUE = logical(1L), inherits, "list")))
        return ("curve.models must be a list of lists")
    if (length(curve.models)>0 && !identical(sapply(curve.models, names), matrix(rep(c("f","tau","type"), length(curve.models)),3)))
        return ("list elt's should be f, tau, type")
    ## xxx additional checks on f's, tau's, type's
    return (TRUE)
}

## xxx use the =match= paradigm and forget about / use try-catch to alter post-condition error checks?
## xxx use OO system to ensure predicates with types
