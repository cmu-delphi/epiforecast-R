
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
    if (class(curve.models) != 'list' || !identical(rep('list',length(curve.models)), as.character(sapply(curve.models, class))))
        return ("curve.models must be a list of lists")
    if (length(curve.models)>0 && !identical(sapply(curve.models, names), matrix(rep(c("f","tau","type"), length(curve.models)),3)))
        return ("list elt's should be f, tau, type")
    ## xxx additional checks on f's, tau's, type's
    return (TRUE)
}

## xxx use the =match= paradigm and forget about / use try-catch to alter post-condition error checks?
## xxx use OO system to ensure predicates with types
