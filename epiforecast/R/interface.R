## author_header begin
## Copyright (C) 2016 Logan C. Brooks, Ryan J. Tibshirani
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
##' @include match.R
##' @include loaders.R
##' @import glmgen
NULL

seasonalGaussian = function(dat, smooth.dat) {
    ensureTRUE(is.or.why.not.smooth.dat(dat, smooth.dat))

    curve.models = lapply(seq_along(dat), function(season.i) {
        trajectory = dat[[season.i]]
        curve = smooth.dat[[season.i]]
        trim.curve = curve[seq_along(trajectory)]
        tau = sqrt(mean((trajectory-trim.curve)^2, na.rm=TRUE))
        return (list(f=curve, tau=tau, type="Gaussian"))
    })

    ensureTRUEpostcondition(is.or.why.not.curve.models(dat, curve.models))
    return (curve.models)
}

##' Function for fitting smooth curves to past seasons' data
##'
##' Arguments:
##' @param dat.obj assumed to be a list, of length equal to number of past
##'   seasons. Each item here is itself a list, each component containing a
##'   vector of "signals" for that seasons
##' @param method one of "ss" and "tf". The former uses R's built-in smoothing
##'   spline method; the latter uses the glmgen package
##' @param cv if TRUE, uses cross-validation to find the smoothing parameter; if
##'   FALSE, uses generalized cross-validation (more efficient, less accurate)
##' @param cv.rule one of "min", "1se", or "gcv": the rule for selecting the
##'   smoothing parameter, where "min" gives the CV usual rule, and "1se" uses
##'   the CV 1-standard-error rule, and "gcv" uses generalized cross-validation
##'   (more efficient, less accurate); the "ss" method accepts only "min" and
##'   "gcv"
##' @param tf.ord the order of the piecewise polynomial fit by trend filtering.
##'   Default is 2
##' @param verbose logical; if TRUE, progress information will be printed out to
##'   the terminal.
##'
##' @return a list with components
##'
##' \code{smooth.obj}: a list of the same dimension as \code{dat.obj}, except
##'   all observed signal values have all been replaced by smoothed values
##'
##' \code{sigma.hat}: a vector of length equal to the number of seasons, each
##'   component giving an estimate of the standard deviation of the noise in
##'   that season's data
##'
##' @export
eb.fitSmoothCurves = function(dat.obj, method=c("ss","tf"),
                              cv.rule=c("min","1se","gcv"), tf.ord=2,
                              verbose=FALSE) {
  method <- match.arg(method)
  cv.rule <- match.arg.else.default(cv.rule)
  ns = length(dat.obj)
  if (ns==0) stop("dat.obj must be of length >= 1")

  smooth.obj = vector(mode="list",length=ns)
  sigma.hat = rep(0,ns)

  for (i in 1:ns) {
    if (verbose) cat(sprintf("Season %i (of %i)\n",i,ns))

    y = dat.obj[[i]]
    x = seq_along(y)
    obs = which(!is.na(y))
    yobs = y[obs]
    xobs = x[obs]
    nobs = length(obs)

    if (nobs < 2) {
      cat(sprintf("Too little data available for season %i;",
                  "returning unsmoothed values ...\n",i))
      smooth.obj[[i]] = y
    } else {
      if (method=="ss") {
        ## fixme check: 1se not allowed for spline
        if (cv.rule == "1se") {
          stop ("The 1se cv.rule is not supported by the ss method.")
        }
        if (verbose) cat("\tFitting smoothing spline object ...\n")
        obj = stats::smooth.spline(xobs,yobs, cv=switch(cv.rule, min=TRUE, gcv=FALSE))
        yhat.obs = stats::predict(obj, xobs)$y
        yhat = stats::predict(obj, x)$y
        df = obj$df
      } else {
        if (verbose) cat("\tFitting trend filtering object ...\n")
        if (cv.rule=="gcv") a = trendfilter.gcv(xobs,yobs,k=tf.ord)
        else a = trendfilter.cv(xobs,yobs,k=tf.ord,cv.rule=cv.rule,verbose=verbose)
        yhat.obs = stats::predict(a$obj, x.new=xobs, lambda=a$lambda.hat)
        yhat = stats::predict(a$obj, x.new=x, lambda=a$lambda.hat)
        df = a$df.hat
      }
      smooth.obj[[i]] = yhat
      sigma.hat[i] = sqrt(sum((yobs-yhat.obs)^2)/(nobs-df))
      if (is.na(sigma.hat[i]) || is.nan(sigma.hat[i])) sigma.hat[i] = sd(yobs)
    }
  }

  return(list(smooth.obj=smooth.obj, sigma.hat=sigma.hat))
}

trendfilter.gcv = function(x, y, k=2, ntimes=3, ...) {
  n = length(y)
  obj = glmgen::trendfilter(x=x, y=y, k=k, ...)
  ## xxx should dispatch to glmnet::predict.glmnet; perhaps it would be better to call directly:
  fits = stats::predict(obj, x.new=x)
  ymat = matrix(rep(y,ncol(fits)),nrow=n)
  crit = colSums((ymat - fits)^2)/(1-obj$df/n)^2

  ## Usual way, find the minimizer, doesn't work well imin = which.min(crit)

  ## Strange version of GCV where we require that it demonstrate ntimes
  ## successive increases before declaring a min point. This works better!
  count = 0
  for (i in 1:(length(crit)-1)) {
    if (crit[i+1] <= crit[i]) count = 0
    else count = count+1
    if (count == ntimes) break
  }
  if (count < ntimes) imin = length(crit)
  else imin = i-ntimes+1

  lambda.hat = obj$lambda[imin]
  df.hat = sum(coef(obj,lambda=lambda.hat)!=0)
  return(list(obj=obj,lambda.hat=obj$lambda[imin],df.hat=obj$df[imin]))
}

trendfilter.cv = function(x, y, k=2, cv.rule=c("min","1se"), nfolds=5,
  verbose=FALSE, ...) {
  cv.rule <- match.arg.else.default(cv.rule)

  obj = glmgen::trendfilter(x=x,y=y,k=k,...)
  lambda = obj$lambda
  nlambda = length(lambda)

  ## In case some of the x,y pairs were binned away (i.e., the case of nearly
  ## identical x values), take those stored by the trendfilter obj
  x = obj$x
  y = obj$y
  n = length(y)

  if (n < nfolds*2+2) stop(sprintf("Need >= %i data points for %i-fold-cross-validation",
                                  nfolds*2+2,nfolds))

  foldid = c(0,rep(1:nfolds,length=n-2),0)
  cvall = matrix(0,nfolds,nlambda)

  for (i in 1:nfolds) {
    if (verbose) cat(sprintf("\tCV fold %i (of %i) ...\n",i,nfolds))

    ite = which(foldid==i)
    itr = which(foldid!=i)
    ytr = y[itr]
    xtr = x[itr]
    yte = y[ite]
    xte = x[ite]

    tmp = glmgen::trendfilter(x=xtr,y=ytr,k=k,lambda=lambda,...)
    prd = stats::predict(tmp, x.new=xte, lambda=lambda)
    cvall[i,] = colMeans((matrix(rep(yte,nlambda),ncol=nlambda) - prd)^2)
  }

  cverr = colMeans(cvall)
  cvse = apply(cvall,2,sd)/sqrt(nfolds)

  imin = which.min(cverr)
  if (cv.rule=="min") {
    i0 = imin
  } else {
    i0 = min(which(cverr<=cverr[imin]+cvse[imin]))
  }
  return(list(obj=obj,lambda.hat=lambda[i0],df.hat=obj$df[i0]))
}

##' Generates a control list for \code{\link{eb.createForecasts}}.
##'
##' With no arguments, returns the default control list.  Optional
##' arguments provided to the function will override these default
##' values (currently with no validation checks).
##'
##' @param parent
##'
##' @param max.n.sims the number of simulated curves to generate in a
##'     forecast
##'
##' @param peak.time.dist the distribution of smoothed-curve peak
##'     times that the prior should follow.  If enabled, each smoothed
##'     curve will be x-shifted to have a peak time which is drawn
##'     from this distribution.  The default setting is to disable
##'     this transformation.  The default enabled distribution is a
##'     discrete uniform distribution fitted to the peak times of the
##'     smoothed curves provided to \code{\link{eb.createForecasts}}
##'     in the argument \code{fit.obj}.
##'
##' @param x.shift.dist the distribution of x-shifts to apply (after
##'     any x-shift from \code{peak.time.dist}).  The default setting
##'     is to enable this transformation.  The default enabled
##'     distribution is a discrete uniform distribution centered at
##'     zero with width equal to twice the bin width of a histogram of
##'     the \code{fit.obj} peak times, using Sturges' rule.
##'
##' @param x.scale.dist
##'
##' @param y.scale.baseline a single numeric value.  Any y-scale
##'     transforms will only transform about and above this baseline
##'     value; for example, for a baseline of 4 and scaling factor of
##'     2, the y-value 1 will not be scaled, since 1<4, and the
##'     y-value 5 will be scale to 5+(5-4)*2=6.  The default is 0,
##'     which, for non-negative smoothed curves, corresponds to simply
##'     multiplying y-values by the scaling factor.
##'
##' @param peak.height.dist the distribution of smoothed-curve peak
##'     heights that the prior should follow.  If enabled, each
##'     smoothed curve will be y-scaled to have a peak height which is
##'     drawn from this distribution.  The default setting is to
##'     disable this transformation.  The default enabled distribution
##'     is a uniform distribution fitted to the \code{fit.obj} peak
##'     heights.  If a smoothed curve remains completely below
##'     \code{y.scale.baseline} the entire time, then y-scaling will
##'     have no effect on that curve, and the peak height will remain
##'     at its original value.  If a peak height selected from the
##'     distribution is lower than \code{y.scale.baseline}, parts of
##'     the curve above the baseline will be scaled by a negative
##'     factor so that the original peak is mapped to the drawn peak
##'     height value; however, this inversion is likely undesirable,
##'     and the resulting peak height may be any value between the
##'     drawn peak height value and the baseline value.
##'
##' @param y.scale.dist the distribution of y-scales to apply (after
##'     any y-scale from \code{peak.height.dist}).  The default
##'     setting is to enable this transformation.  The default enabled
##'     transformation is a log-uniform distribution centered at 0 in
##'     the log-scale with log-scale width equal to twice the bin
##'     width of a histogram of the logarithms of the \code{fit.obj}
##'     peak heights, using Sturges' rule.  Note that this default
##'     behavior can significantly bias the mean of the prior for the
##'     peak heights, but does not significantly affect the median of
##'     the prior for the peak heights.
##'
##' @param sd.option one of \code{"match"}, \code{"scale"}, or
##'     \code{"prior"}, which controls the assignment of a single
##'     noise level, or distribution of possible noise levels, to a
##'     transformed curve: \code{"match"} chooses the \code{sigma.hat}
##'     associated with the selected smooth curve; \code{"scale"} does
##'     the same, but scales this \code{sigma.hat} by the y-scale
##'     factor given by the transformations selected from
##'     \code{peak.height.dist} and \code{y.scale.dist};
##'     \code{"prior"} selects a noise level uniformly from the
##'     \code{sigma.hat}'s of all smoothed curves fed into the EB
##'     method, not just the one corresponding to the current
##'     transformed curve.
##'
##' @param sd.prior controls the distribution of noise levels used
##'     when \code{sd.option} is \code{"prior"}; currently, the only
##'     choice is \code{"uniform"}.
##'
##' @param sd.scale.dist controls the distribution of noise level
##'     scaling factors, which are applied after \code{sd.option} and
##'     \code{sd.prior} are used to select an initial noise level.
##'     The default setting is to disable this transformation.  There
##'     is no default enabled distribution.
##'
##' @param reasonable.future.weight controls the coefficient of the
##'     "reasonable future" term added to the conditional
##'     log-likelihood of the observed values given a transformed
##'     curve and noise level when calculating importance weights.
##'     The default value is 0, which disables this feature.
##'
##' @param n.future.neighbors controls the number of neighbors to use when
##'   determining the "reasonable future" term: for a given transformed curve
##'   and noise parameter, the neighbors are the \code{n.future.neighbors}
##'   historical noisy curves from the \code{dat} argument in
##'   \code{\link{eb.createForecasts}} with the highest log-likelihoods in
##'   future weeks (after \code{time.of.forecast}; the reasonable future term is
##'   the average across these neighbors of the log-likelihood in future weeks.
##'
##' @param bias_peaktime_weight is the coefficient of a "peak time bias" term
##'   added to the conditional log-likelihood of the observed values given a
##'   transformed curve and noise level when calculating importance weights. The
##'   default value is \code{NA}, which disables this feature. When it is
##'   enabled, the peak time bias term is the log-likelihood of the underlying
##'   peak time (before adding noise) of a curve, according to a Gaussian
##'   distribution with mean \code{bias_peaktime_mean} and standard deviation
##'   \code{bias_peaktime_sd}, multiplied by \code{bias_peaktime_weight}.
##'
##' @param bias_peaktime_mean the mean of the peak time bias Gaussian
##'   distribution
##' 
##' @param bias_peaktime_sd the standard deviation of the peak time bias
##'   distribution
##'
##' @param bias_peakheight_weight is the coefficient of a "peak height bias"
##'   term added to the conditional log-likelihood of the observed values given
##'   a transformed curve and noise level when calculating importance weights.
##'   The default value is \code{NA}, which disables this feature. When it is
##'   enabled, the peak height bias term is the log-likelihood of the underlying
##'   peak height (before adding noise) of a curve, according to a Gaussian
##'   distribution with mean \code{bias_peakheight_mean} and standard deviation
##'   \code{bias_peakheight_sd}, multiplied by \code{bias_peakheight_weight}.
##' 
##' @param bias_peakheight_mean the mean of the peak height bias Gaussian
##'   distribution
##' 
##' @param bias_peakheight_sd the standard deviation of the peak height bias
##'   distribution
##'
##' @param inactive.seasons is currently ignored.
##'
##' @param n.out is the number of observations that each outputted
##'     noisy curve should contain; it should be less than or equal to
##'     the length of the shortest smooth curve.  For weekly data and
##'     year-long seasons, this should be 52 or 53.
##'
##' @param ii.match.mask is a vector of indices in
##'     \code{seq_len(n.out)}; only observations at these times will
##'     be considered when computing the likelihood of observations
##'     and assigning "reasonable future" terms to transformed curves
##'     and noise levels.
##'
##' @param max.match.length is a single integer controlling the
##'     maximum number of observations to use when computing the
##'     log-likelihood of \code{new.dat} given a transformed curve and
##'     noise level; if more than \code{max.match.length} observations
##'     are available at \code{time.of.forecast} after applying
##'     \code{ii.match.mask}, only the \code{max.match.length} most
##'     recent observations are used in the likelihood calculation.
##'
##' @param n.unpinned.observations is a single integer controlling
##'     what values in the noisy transformed curves in the posterior
##'     are "pinned" to the observed values in \code{new.dat}; any
##'     observations after
##'     \code{time.of.forecast-n.unpinned.observations} are not
##'     pinned.
##'
##' @details Most settings are single integers, single reals, or
##'     character vectors where the first entry holds the desired
##'     value.  Transformation distribution settings, on the other
##'     hand, can be one of several options, with the following
##'     associated meanings:
##'     \itemize{
##'       \item{\code{NULL}: }{Use the default setting for this
##'       transformation, either disabling the transformation or using
##'       the default distribution}
##'       \item{\code{TRUE}: }{Enable this transformation, and use the
##'       default distribution for this transformation}
##'       \item{\code{FALSE}: }{Disable this transformation}
##'       \item{Single integer: }{Enable this transformation and use the
##'       default distribution for this transformation, but break the
##'       distribution into the specified number of bins (rather than
##'       the default) when applying the grid importance sampling
##'       algorithm used by \code{\link{eb.createForecasts}}}
##'       \item{Distribution with bins:}{}
##'       \item{Function from curve.obj to distribution with bins:}{}
##'     }
##'
##' @return a list of parameter settings used by
##'     \code{\link{eb.createForecasts}}
##'
##' @examples
##'   default.control.list get_eb_control_list()
##'   with.less.sims = get_eb_control_list(max.n.sims = 10000L)
##'   with.less.sims.another.way = get_eb_control_list(default.control.list, max.n.sims = 10000L)
##'   with.less.sims.and.sd.option.scale = get_eb_control_list(with.less.sims, sd.option="scale")
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
get_eb_control_list = function(parent = NULL,
                               max.n.sims = 2000L,
                               peak.time.dist = NULL,
                               x.shift.dist = NULL,
                               x.scale.dist = NULL,
                               y.scale.baseline = 0,
                               peak.height.dist = NULL,
                               y.scale.dist = NULL,
                               sd.option = c("match", "scale", "prior"),
                               sd.prior = "uniform", # xxx alternatives?
                               sd.scale.dist = NULL,
                               reasonable.future.weight = 0,
                               n.future.neighbors = 3L,
                               bias_peaktime_weight = NA_real_,
                               bias_peaktime_mean = NA_real_,
                               bias_peaktime_sd = NA_real_,
                               bias_peakheight_weight = NA_real_,
                               bias_peakheight_mean = NA_real_,
                               bias_peakheight_sd = NA_real_,
                               ## xxx weight.future.neighbors = 1, ???
                               inactive.seasons = NULL, # xxx impl
                               n.out=53L,
                               ii.match.mask = NULL, # xxx should also allow boolean vec, or perhaps numeric for partial weighting (shrinking), and check for usage in futweight
                               max.match.length = NULL, # xxx NULL or integer (or numeric if using partial weighting)
                               n.unpinned.observations = 0L,
                               model = "Empirical Bayes"
                               ) {
    fresh.bindings = as.list(if (is.null(parent)) {
                                 environment()
                             } else {
                                 match.call()[-1]
                             })
    fresh.bindings <- fresh.bindings[-match("parent", names(fresh.bindings))]
    control.list = parent
    for (binding.i in seq_along(fresh.bindings)) {
        ## todo validation checks, do most of the input.eb.control.list here
        binding.name = names(fresh.bindings)[binding.i]
        binding.value = fresh.bindings[[binding.i]]
        control.list[binding.name] <- list(eval.parent(binding.value)) # stores NULL's too
    }
    return (control.list)
}





##' Returns a possibly-named length-1 non-NA integer-class vector
##' version of the input, or produces an error if the input seems
##' inappropriate.  Designed for input processing and validation
##' within another function.
match.single.integer = function(n) {
    if (is.numeric(n) && length(n)==1L && !is.na(n) && floor(n)==n)
        return (as.integer(n))
    else
        stop(sprintf("Argument %s must be a length 1 numeric vector representing an integral value.", as.character(match.call()$n)))
}

##' Returns a distribution that can be divided into buckets --- a
##' named list containing \code{`n`}, a single integer representing
##' the number of buckets into which the distribution can be broken
##' down; \code{`choices`}, an \code{n}-length vector containing a
##' representative element from each bucket; \code{`probs`}, an
##' \code{n}-length vector containing the probability mass assigned to
##' each bucket, and \code{`sampler`}, a function from vectors of
##' bucket indices (in \code{seq_len(n)}) to randomly sampled elements
##' within the corresponding buckets.  Designed for input processing
##' and validation within another function.
##'
##' @param curve.obj output of \code{\link{eb.fitSmoothCurves}}
##'
##' @param dist one of the following: (a) \code{NULL}, (b) a single boolean, (c) a
##'     single integer, (d) a function that outputs a distribution given a
##'     \code{curve.obj}, or (e) a distribution.
##'
##' @param null.replacement one of (b)--(e), which is used to replace
##'     \code{NULL} inputs
##'
##' @param true.replacement one of (c)--(e), which is used to replace
##'     \code{TRUE} inputs
##'
##' @param false.replacement one of (c)--(e), which is used to replace
##'     \code{FALSE} inputs
##'
##' @param integer.replacement.fn a function from a single integer to
##'     either (d) or (e), called on integer inputs to produce a
##'     replacement value
##'
##' @return a distribution, which incorporates any fitting procedure
##'     to \code{curve.obj};
match.dist = function(curve.obj, dist,
                      null.replacement,
                      true.replacement, false.replacement,
                      integer.replacement.fn
                      ) {
    if (is.null(dist))
        dist <- null.replacement

    if (isTRUE(dist))
        dist <- true.replacement
    else if(identical(dist,FALSE))
        dist <- false.replacement

    if (is.numeric(dist) && length(dist)==1L && floor(dist)==dist) {
        ## dist <- {integer.replacement.fn(as.integer(dist))}
        n.bins = as.integer(dist) # rename seems required
        dist <- integer.replacement.fn(n.bins)
    }

    ## xxx match.fun
    if (is.function(dist))
        dist <- dist(curve.obj)

    if (is.list(dist) && all(c("n","choices","probs","sampler")%in%names(dist)))
        return (dist)

    stop(sprintf("Distribution argument %s must be NULL (default), a single boolean (enable/disable), a single integer (bin count), distribution (list of n,choices,probs,sampler), or a function from a curve.obj to a distribution.", as.character(match.call()$dist)))
}

fit.eb.control.list.old = function(curve.obj, control.list) {
    return (with(control.list, { # access inputted value of `control.list$thing` with `thing`
        ## Almost every statement here should have the form "control.list$thing <- expression(curve.obj, thing)", which takes the input `thing`, validates input, and assigns a fit and standardized version to the control list (with `thing` usually retaining its inputted value).

        ## =n.out= and =max.n.sims= should be single-length integers:
        control.list$n.out <- n.out <- match.single.integer(n.out) # also overwrite local =n.out=
        control.list$max.n.sims <- match.single.integer(max.n.sims)

        ## Fill in defaults and fit distributions according to the specification in =get_eb_control_list=:
        ## Peak time dist: default: disabled; enabled default: fit a discrete uniform; disabled default: uniform distribution over =NA=, interpreted specially by the EB procedure; enabled default with specified number of buckets: not yet implemented.
        control.list$peak.time.dist <- match.dist(
            curve.obj, peak.time.dist,
            FALSE, function(newfit) {
                unifLocGridPrior(newfit.to.oldfit(newfit))
            }, unifChoicePrior(NA), stop("Specifying n.bins with integer-domain distributions is not yet implemented.")) # todo implement
        ## x-shift dist: default: enabled; enabled default: uniform with width twice the Sturges histogram bin width; disabled default: uniform distribution over =NA=, interpreted specially by the EB procedure; enabled default with specified number of buckets: evenly divide the domain of the distribution into buckets.
        control.list$x.shift.dist <- match.dist(
            curve.obj, x.shift.dist,
            TRUE, 31L, unifChoicePrior(0), function(n.bins) function(newfit) {
                ## =diff(hist(...))= gives widths of histogram bins, which may be different due to floating-point arithmetic; use =mean= rather than =unique= to get a single bin-width:
                x.shift.end = round(mean(diff(hist(sapply(newfit, function(season.fit) {
                    which.max(season.fit$f)}),plot=FALSE)$breaks)))
                x.shift.begin = -x.shift.end
                unifChoicePrior(x.shift.begin:x.shift.end)
            })
        control.list$x.scale.dist <- match.dist(
            curve.obj, x.scale.dist,
            FALSE, 11L, unifChoicePrior(1), function(n.bins) function(newfit) {
                logUnifDistrLogUnifGridPrior(0.75, 1.25, n.bins)
            })
        control.list$y.scale.baseline <-
            if (is.numeric(y.scale.baseline) && length(y.scale.baseline)%in%c(1L,n.out)) {
                if (length(y.scale.baseline) != 1L) {
                    stop("Length>1 baselines not implemented yet.") # todo implement
                } else y.scale.baseline
            } else stop("y.scale.baseline must be a numeric of length 1 or length n.out.")
        control.list$peak.height.dist <- match.dist(
            curve.obj, peak.height.dist,
            FALSE, 61L, unifChoicePrior(NA), function(n.bins) function(newfit) {
                unifScaleUnifGridPrior(newfit.to.oldfit(newfit), n.bins)
            })
        control.list$y.scale.dist <- match.dist(
            curve.obj, y.scale.dist,
            TRUE, 31L, unifChoicePrior(1), function(n.bins) function(newfit) {
                ## todo include option for mean correction? nonlogscale?
                y.scale.end = exp(mean(diff(hist(log(sapply(newfit, function(season.fit) {
                    max(season.fit$f)})),plot=FALSE)$breaks)))
                y.scale.begin = 1/y.scale.end
                logUnifDistrLogUnifGridPrior(y.scale.begin, y.scale.end, n.bins)
            })
        control.list$sd.option <-
            if (is.character(sd.option) && length(sd.option)>=1L
                && sd.option[1L]%in%get_eb_control_list()$sd.option) {
                sd.option[1L]
            } else stop("sd.option must be a length>=1 character vector with first entry in get_eb_control_list()$sd.option")
        if (!identical(sd.prior, "uniform")) stop('sd.prior must be "uniform"')
        control.list$sd.scale.dist <- match.dist(
            curve.obj, sd.scale.dist,
            FALSE, -1L, unifChoicePrior(1), stop("There is no default enabled sd.scale dist."))

        if (!is.numeric(reasonable.future.weight) || length(reasonable.future.weight)!=1L || reasonable.future.weight<0 || reasonable.future.weight>1)
            stop("reasonable.future.weight must be a length-1 numeric between 0 and 1.")
        control.list$n.future.neighbors <- match.single.integer(n.future.neighbors)

        control.list$bias_peaktime_mean = match.single.na.or.numeric(bias_peaktime_mean)-1L
        control.list$bias_peaktime_sd = match.single.na.or.numeric(bias_peaktime_sd)
        control.list$bias_peaktime_weight = match.single.na.or.numeric(bias_peaktime_weight)
        control.list$bias_peakheight_mean = match.single.na.or.numeric(bias_peakheight_mean)
        control.list$bias_peakheight_sd = match.single.na.or.numeric(bias_peakheight_sd)
        control.list$bias_peakheight_weight = match.single.na.or.numeric(bias_peakheight_weight)

        if (!is.null(inactive.seasons)) stop("inactive seasons not implemented yet") # todo impl

        control.list$ii.match.mask <-
            if (is.null(ii.match.mask)) {
                seq_len(n.out)
            } else if (is.numeric(ii.match.mask) && all(floor(ii.match.mask)==ii.match.mask)
                       && all(ii.match.mask>=1L & ii.match.mask<=n.out)) {
                sort(as.integer(ii.match.mask))
            } else stop("ii.match.mask must be a numeric vector with values in seq_len(n.out).")
        control.list$max.match.length <-
            if (is.null(max.match.length)) { n.out
            } else match.single.integer(max.match.length) # xxx want NULL possibility in error message
        control.list$n.unpinned.observations <- match.single.integer(n.unpinned.observations)
        control.list
    }))
}
## xxx transformation.ss

##' Takes an EB control list containing arguments that may require fitting to a
##' curve object, performs any fitting, and outputs a static EB control list
##' containing the results of the fitting procedure. The contents of the EB
##' control list are also validated and standardized to a more rigid form, e.g,
##' replacing some \code{NULL} values with defaults and some non-integer-class
##' integral input with integer-class versions.
fit.eb.control.list = function(dat, new.dat, fit.obj, time.of.forecast, control.list) {
    fit.control.list = fit.eb.control.list.old(fit.obj, control.list)
    ii.match.candidates = seq_len(time.of.forecast)
    ii.match.candidates <- ii.match.candidates[!is.na(new.dat[ii.match.candidates])]
    ii.match = tail(intersect(ii.match.candidates, fit.control.list$ii.match.mask),
                    fit.control.list$max.match.length)
    ii.future.candidates = Seq(time.of.forecast+1, min(control.list$n.out, sapply(dat, length)))
    if (length(ii.future.candidates) > 0) { ## length 0 gives undesired simplify2array list output
        ii.future.candidates <- ii.future.candidates[!rowSums(is.na(sapply(dat,`[`,ii.future.candidates)))]
    }
    ii.future = intersect(ii.future.candidates, fit.control.list$ii.match.mask)

    fit.control.list$ii.match <- ii.match
    fit.control.list$ii.future <- ii.future
    fit.control.list$curve.index.dist <- unifChoicePrior(seq_along(fit.obj))
    ## fit.control.list$sd.dist <- unifChoicePrior(sapply(fit.obj,`[[`,"tau"))
    fit.control.list$sd.dist <- unifChoicePrior(NA)

    return (fit.control.list)
}



##' Function for making forecasts with the empirical Bayes method.
##'
##' @template sim.method_template
##'
##' @param fit.obj optional argument to forward to
##'   \code{\link{eb.createForecasts}}
##' @param control.list optional argument to forward to
##'   \code{\link{eb.createForecasts}}; \code{n.out} is overridden with
##'   length(new.dat)
##' @param max.n.sims maximum number of (weighted) simulated curves to produce.
##'   Defaults to 1000.
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
eb.sim = function(full.dat, baseline=0, time.of.forecast = NULL, max.n.sims = 1000L,
                  fit.obj=NULL, control.list=get_eb_control_list()) {

    ## Check input
    ## check.list.format(full.dat)

    ## Split into old dat (list) and new dat (vector)
    old.dat = head(full.dat, -1L)
    new.dat.sim = match.new.dat.sim(tail(full.dat, 1L)[[1]])
    ## xxx Using multiple trajectories from a new.dat.sim is too slow (without
    ## revisiting the Rcpp code); collapse new.dat.sim into a single new.dat
    ## trajectory:
    new.dat = matrixStats::rowWeightedMeans(new.dat.sim[["ys"]], new.dat.sim[["weights"]])
    old.season.labels = head(names(full.dat), -1L)
    new.season.label = tail(names(full.dat), 1L)

    ## Make forecast
    if(is.null(time.of.forecast)) time.of.forecast = get.latest.time(new.dat)
    control.list = get_eb_control_list(parent = control.list,
                                       n.out=length(new.dat),
                                       y.scale.baseline = baseline,
                                       max.n.sims = max.n.sims)
    sim = eb.createForecasts(old.dat, new.dat, fit.obj, time.of.forecast,
                             control.list = control.list)

    ## Append the control parameter list
    sim = c(sim, control.list = list(control.list))

    ## Append old and new dat
    sim = c(sim,
            old.dat = list(old.dat),
            new.dat = list(new.dat),
            old.season.labels = list(old.season.labels),
            new.season.label = list(new.season.label))

    ## Cast it as class |sim| before returning.
    class(sim) <- "sim"

    return (sim)
}



##' Generates control list for BR
##' @param smooth logical; if TRUE, past observations and future
##'   "pseudo-observations" (predictions) will be smoothed; if FALSE, the
##'   observations and pseudo-observations will be returned unsmoothed.
##' @param model.noisy logical; if TRUE and bootstrapping is enabled, will
##'   inject noise into future values and pin past observations to the observed
##'   values; if FALSE or TRUE but used when bootstrapping is disabled, br.sim
##'   will return the "non-noisy" regression curve
##' @param basis type of basis to use. So far only "bs" (B-splines) are
##'   implemented.
##' @param df the degrees of freedom for the basis. Default is 10.
##' @param w the mixing weight between the two loss terms, as in: sum over obs
##'   times (yobs - f)^2 + w * sum over unobs times (ypast - f)^2, where yobs is
##'   the current season's observed data, ypast is the past season's data,
##'   suitably transformed, and f is the function to be estimated.
##' @param scale.method whether and how to scale past seasons to match data from
##'   the \code{cur.season}th trajectory: "none" performs no scaling; "max"
##'   scales the maximum of each other season's trajectory --- restricted to
##'   times which correspond to non-\code{NA} values in the \code{cur.season}th
##'   trajectory --- so that it matches the maximum of the \code{cur.season}th
##'   trajectory; and "last" performs the same scaling using data at the time
##'   corresponding to the latest observation in the \code{cur.season}th
##'   trajectory
##' @param max.scale.factor single numeric: a limit on the amount of scaling
##'   performed by the scaling method: scale factors over
##'   \code{max.scale.factor} and under \code{1/max.scale.factor} will be
##'   clipped.
##' @param max.match.length the maximum number of past data points to which the
##'   spline is fitted. The default, NULL, indicates to use all past data points
##'   when fitting the spline.
##' @param cv.rule one of "min" or "1se", where "min" gives the usual rule, and
##'   "1se" uses the 1-standard-error rule.
##' @param baseline the anchoring point used for scaling past season's data;
##'   data above the baseline are scaled about the baseline. The default value,
##'   NA, indicates to scale about 0 (regardless of sign).
##' @export
get_br_control_list = function(parent = NULL,
                               max.n.sims = 100L,
                               n.out=53L,
                               model = "Basis Regression",
                               max.match.length=NULL,
                               df = 10,
                               w = 1,
                               smooth=TRUE,
                               model.noisy=TRUE,
                               baseline=0,
                               basis="bs", ## df.to.weight.ratio=1/5,
                               max.scale.factor=3,
                               cv.rule=c("min","1se"),
                               scale.method=c("none","max","last")) {

    ## Fetching + bundling all function arguments into a list.
    fresh.bindings = as.list(if (is.null(parent)) {
                                 environment()
                             } else {
                                 match.call()[-1]
                             })
    fresh.bindings <- fresh.bindings[-match("parent", names(fresh.bindings))]
    control.list = parent
    for (binding.i in seq_along(fresh.bindings)) {
        binding.name = names(fresh.bindings)[binding.i]
        binding.value = fresh.bindings[[binding.i]]
        control.list[binding.name] <- list(eval.parent(binding.value)) # stores NULL's too
    }

    ## Manually doing 'match.arg()' on cv.rule and scale.method
    if(length(control.list$cv.rule)>1) control.list$cv.rule = (control.list$cv.rule)[1]
    if(length(control.list$scale.method)>1) control.list$scale.method = (control.list$scale.method)[1]
    return (control.list)
}
