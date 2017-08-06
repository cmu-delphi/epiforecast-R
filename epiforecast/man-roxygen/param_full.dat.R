##' @param full.dat list of (a) numeric vectors, one per past season, containing
##'   historical trajectories, followed by (b) a numeric vector (trajectory),
##'   numeric matrix (cbound trajectories), or sim object (list with $ys a
##'   numeric matrix (cbound trajectories) and $weights a numeric vector
##'   (associated weights)), with \code{NA}'s for all future or missing data
##'   points to forecast or infer; currently only supports \code{NA}'s at future
##'   points, not mixed in between non-\code{NA} data
