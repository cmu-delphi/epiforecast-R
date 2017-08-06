##' @template param_full.dat
##' @template param_baseline
##' @param max.n.sims single non-\code{NA} integer value or \code{NULL}: the
##'   number of curves to sample from the inferred distribution
##' @return a sim object --- a list with two components:
##'
##' \code{ys}: a numeric matrix, typically with multiple columns; each column is
##' a different possible trajectory for the current season, with NA's in the
##' input for the current season filled in with random draws from the forecasted
##' distribution, and non-\code{NA}'s (observed data) filled in with an imagined
##' resampling of noise based on the model (for some models, the non-\code{NA}
##' values will remain unchanged).
##'
##' \code{weights}: a numeric vector; assigns a weight to each column of
##' \code{ys}, which is used by methods relying on importance sampling.
