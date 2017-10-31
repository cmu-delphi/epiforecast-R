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

##' @include match.R
##' @include loaders.R
##' @include interface.R
##' @include eb_dists.R
##' @import Rcpp
##' @useDynLib epiforecast
NULL

obs.vec.to.trajectory = function(obs.vec) {
  data.frame(time=seq_along(obs.vec)-1, observation=obs.vec)
}

## fixme documentation

## Produce log-weights of hypergrid of buckets of parameters; this forms a
## histogram-like proposal distribution for importance sampling.
eb.log.bucket.weights = function(dat, new.dat, fit.obj, time.of.forecast, fit.control.list) {
  ## indices of observations to condition on:
  ii.match = fit.control.list$ii.match
  ## indices used in "reasonable future" subobjective:
  ii.future = fit.control.list$ii.future
  ## Log CONDITIONAL bucket weights given the data (not weighted by prior):
  log.cond.bucket.weights =
    CartesianProductLogWeights(
      obs.vec.to.trajectory(new.dat),
      lapply(dat, obs.vec.to.trajectory),
      data.frame(time=ii.match-1, weight=rep(1,length(ii.match))),
      data.frame(time=ii.future-1, weight=rep(fit.control.list$reasonable.future.weight,length(ii.future))),
      fit.control.list$n.future.neighbors,
      fit.obj,
      fit.control.list$y.scale.baseline,
      fit.control.list$curve.index.dist$choices - 1,
      fit.control.list$peak.time.dist$choices,
      fit.control.list$x.shift.dist$choices,
      fit.control.list$x.scale.dist$choices,
      fit.control.list$sd.dist$choices,
      fit.control.list$sd.scale.dist$choices,
      fit.control.list$peak.height.dist$choices,
      fit.control.list$y.scale.dist$choices,
      fit.control.list$bias_peaktime_mean,
      fit.control.list$bias_peaktime_sd,
      fit.control.list$bias_peaktime_weight,
      fit.control.list$bias_peakheight_mean,
      fit.control.list$bias_peakheight_sd,
      fit.control.list$bias_peakheight_weight
    )
  ## set friendly dimnames:
  dn = rev(list(
    ## y.scale.baseline=fit.control.list$y.scale.baseline,
    curve.index.dist=fit.control.list$curve.index.dist$choices - 1,
    peak.time.dist=fit.control.list$peak.time.dist$choices,
    x.shift.dist=fit.control.list$x.shift.dist$choices,
    x.scale.dist=fit.control.list$x.scale.dist$choices,
    sd.dist=fit.control.list$sd.dist$choices,
    sd.scale.dist=fit.control.list$sd.scale.dist$choices,
    peak.height.dist=fit.control.list$peak.height.dist$choices,
    y.scale.dist=fit.control.list$y.scale.dist$choices
  ))
  dim(log.cond.bucket.weights) <- sapply(dn,length)
  dimnames(log.cond.bucket.weights) <- dn

  ## Prior over parameters:
  col.major.dists = rev(list(
    curve.i=unifChoicePrior(seq_along(fit.obj)),
    peak.time.dist=fit.control.list$peak.time.dist,
    x.shift.dist=fit.control.list$x.shift.dist,
    x.scale.dist=fit.control.list$x.scale.dist,
    sd.dist=fit.control.list$sd.dist,
    sd.scale.dist=fit.control.list$sd.scale.dist,
    peak.height.dist=fit.control.list$peak.height.dist,
    y.scale.dist=fit.control.list$y.scale.dist
  ))

  ## Log PRIOR bucket likelihoods:
  log.prior.bucket.lkhds = Reduce(function(subgrid.bucket.weights, parm.dist)
    structure(kronecker(log(parm.dist$probs), subgrid.bucket.weights, FUN=`+`),
              dim=c(dim(subgrid.bucket.weights), parm.dist$n)),
    col.major.dists, 0)

  ## Combine PRIOR + CONDITIONAL:
  log.bucket.weights = log.prior.bucket.lkhds + log.cond.bucket.weights
  return (log.bucket.weights)
}

## Infer EB parameters, mean+sd curves; gives list of mean curves, sd curves,
## hyperparemeters, and importance sampling weights.
##
## Uses hypergrid histogram importance sampling:
## 1. Find weights for all parameter combinations in a coarse hypergrid.
## 2. Interpret these weights as "heights" in a histogram-shaped proposal
##    distribution over parameter configurations.
## 3. Sample proposed parameter buckets of the histogram from the proposal histogram.
## 4. Sample proposed parameter values within the proposed parameter buckets.
## 5. Build curves associated with the proposed parameter values.
## 6. Find corresponding importance weights.
eb.infer = function(dat, new.dat,
                    fit.obj = NULL,
                    time.of.forecast = NULL,
                    control.list = get.eb.control.list()
                    ) {
  ## Fit smooth curves to trajectories in dat unless =fit.obj= has been passed
  ## in:
  if (is.null(fit.obj))
    fit.obj <- smooth.curves.to.fit(eb.fitSmoothCurves(dat))
  ## Assuming all NA's contiguous, the time of forecast for a partial trajectory
  ## is the last non-NA measurement or 0L if all measurements are missing:
  if (is.null(time.of.forecast))
    time.of.forecast <- sum(!is.na(new.dat))

  ## Check input:
  ensureTRUE(is.or.why.not.curve.models(dat, fit.obj))
  time.of.forecast <- match.single.integer(time.of.forecast)
  ## xxx combine error msgs
  if (time.of.forecast < 0L || time.of.forecast > length(new.dat))
    stop("time.of.forecast must be in c(0, seq_along(new.dat))")

  ## Fit/choose control list parameters based on data:
  fit.control.list <- fit.eb.control.list(dat, new.dat, fit.obj, time.of.forecast, control.list)

  ## Extract parameters from fit control list:
  ii.match = fit.control.list$ii.match
  ii.future = fit.control.list$ii.future
  ## List of prior distributions over "parameters":
  col.major.dists = rev(list(
    curve.i=unifChoicePrior(seq_along(fit.obj)),
    peak.time.dist=fit.control.list$peak.time.dist,
    x.shift.dist=fit.control.list$x.shift.dist,
    x.scale.dist=fit.control.list$x.scale.dist,
    sd.dist=fit.control.list$sd.dist,
    sd.scale.dist=fit.control.list$sd.scale.dist,
    peak.height.dist=fit.control.list$peak.height.dist,
    y.scale.dist=fit.control.list$y.scale.dist
  ))

  ## Log-weights of hypergrid of buckets of parameters; this forms a
  ## histogram-like proposal distribution for importance sampling:
  log.bucket.weights = eb.log.bucket.weights(dat, new.dat, fit.obj, time.of.forecast, fit.control.list)

  ## Scale proposal weights so we can bring them out of the log domain:
  scaled.bucket.weights = exp(log.bucket.weights - max(log.bucket.weights))
  ## Sample bucket indices from histogram proposal distribution:
  sampled.bucket.inds = sample(seq_along(scaled.bucket.weights), fit.control.list$max.n.sims, prob=scaled.bucket.weights, replace=TRUE)
  ## xxx should sort inds and optimize for repeated discrete parm configs in C++ code
  ## Convert vector (1-D) indices to array (#parameter-D) indices:
  sampled.bucket.rayinds = arrayInd(sampled.bucket.inds, dim(log.bucket.weights))
  ## Sample within the selected parameter BUCKETS to get the proposed
  ## parameter VALUES:
  sampled.row.major.parms = rev(structure(lapply(seq_along(col.major.dists), function(dist.i)
    col.major.dists[[dist.i]]$sampler(sampled.bucket.rayinds[,dist.i])
    ), names=names(col.major.dists)))

  ## Given n proposed hyperparemeter configurations, build n data trajectories
  ## and n associated importance sampling weights:
  sampled.curves.and.log.weights = ZipProductCurvesAndLogWeightsp(
    seq_len(fit.control.list$n.out)-1,
    obs.vec.to.trajectory(new.dat),
    lapply(dat, obs.vec.to.trajectory),
    data.frame(time=ii.match-1, weight=rep(1,length(ii.match))),
    data.frame(time=ii.future-1, weight=rep(fit.control.list$reasonable.future.weight,length(ii.future))),
    fit.control.list$n.future.neighbors,
    fit.obj,
    fit.control.list$y.scale.baseline,
    sampled.row.major.parms$curve.i - 1,
    sampled.row.major.parms$peak.time.dist,
    sampled.row.major.parms$x.shift.dist,
    sampled.row.major.parms$x.scale.dist,
    sampled.row.major.parms$sd.dist,
    sampled.row.major.parms$sd.scale.dist,
    sampled.row.major.parms$peak.height.dist,
    sampled.row.major.parms$y.scale.dist,
    fit.control.list$bias_peaktime_mean,
    fit.control.list$bias_peaktime_sd,
    fit.control.list$bias_peaktime_weight,
    fit.control.list$bias_peakheight_mean,
    fit.control.list$bias_peakheight_sd,
    fit.control.list$bias_peakheight_weight
  )

  ## Return list of n mean curves, n sd curves, n log-importance-weights, and n
  ## parameter configurations:
  return (list(
    mean=sampled.curves.and.log.weights[[1]],
    sd=sampled.curves.and.log.weights[[2]],
    log.weights=sampled.curves.and.log.weights[[3]],
    sampled.row.major.parms=sampled.row.major.parms
  ))
  ## return (sampled.curves.and.log.weights)
}

##' Function for making forecasts with the empirical Bayes method.
##'
##' @param dat a list of numeric vectors, one per past season, containing
##'   historical trajectories.
##' @param new.dat a single numeric vector containing the observations for the
##'   current season so far, and possibly future data points as well (when
##'   performing retrospective analysis); should not contain any NA's.
##' @param fit.obj a collection of fit curves and noise level estimates to use
##'   when forming the prior; defaults to
##'   \code{smooth.curves.to.fit(eb.fitSmoothCurves(dat))}; while the smoothing
##'   method is quite fast, repeated calls to \code{eb.createForecasts} may
##'   benefit from caching the smoothed curves and feeding them in each time.
##' @param time.of.forecast integer in [0..\code{length(new.dat.partial)}]; if
##'   specified, the forecast is prepared as if
##'   \code{new.dat.partial[seq_len(time.of.forecast)]} was fed in.
##' @param control.list optional control list to forward to
##'   \code{\link{eb.createForecasts}}.
##'
##' @return a list with two components:
##'
##' \code{ys}: a numeric matrix; in most other methods, each column is a different
##' possible trajectory for the current season, with NA's in new.dat filled in
##' with random draws from the forecasted distribution, and non-NA's (observed
##' data) filled in with an imagined resampling of noise based on the model.
##'
##' \code{weights}: a numeric vector; assigns a weight to each column of
##' \code{ys}, which is used by methods relying on importance sampling.
##'
##' @export
eb.createForecasts = function(dat, new.dat,
                              fit.obj = NULL,
                              time.of.forecast = NULL,
                              control.list = get.eb.control.list()
                              ) {
  if (is.null(time.of.forecast))
    time.of.forecast <- length(new.dat)

  sampled.curves.and.log.weights = eb.infer(dat, new.dat, fit.obj, time.of.forecast, control.list)

  ys = structure(stats::rnorm(control.list$n.out*control.list$max.n.sims,
                              sampled.curves.and.log.weights[[1L]],
                              sampled.curves.and.log.weights[[2L]]),
                 dim=c(control.list$n.out, control.list$max.n.sims))
  ys[seq_len(time.of.forecast),] <- new.dat[seq_len(time.of.forecast)]

  log.weights = sampled.curves.and.log.weights[[3L]]
  scaled.weights = exp(log.weights-max(log.weights))

  sim = list(ys=ys,
             ## ## make total weight equal to length(dat):
             ## weights=scaled.weights*length(dat)/sum(scaled.weights),
             ## make total weight equal to length(weights):
             weights = scaled.weights/mean(scaled.weights),
             sampled.row.major.parms=sampled.curves.and.log.weights[[4L]])

  return (sim)
}

eb.cv.future.data = function(dat, fit.obj = NULL,
                             time.of.forecast,
                             control.list = get.eb.control.list()) {
  if (is.null(fit.obj))
    fit.obj <- smooth.curves.to.fit(eb.fitSmoothCurves(dat))
  if (length(dat) != length(fit.obj))
    stop("dat trajectories must correspond to fit.obj curves")
  eval.ii = Seq(time.of.forecast+1, min(control.list$n.out, sapply(dat, length)))
  fold.eval.log.weights = sapply(seq_along(dat), function(s) {
    new.dat = dat[[s]]
    sampled.curves.and.log.weights = eb.infer(dat[-s], new.dat, fit.obj[-s], time.of.forecast, control.list)
    curve.eval.log.weights = colSums(stats::dnorm(new.dat[eval.ii],
                                                  sampled.curves.and.log.weights$mean,
                                                  sampled.curves.and.log.weights$sd,
                                                  log=TRUE))
    curve.weights = exp(sampled.curves.and.log.weights$log.weights
                        - max(sampled.curves.and.log.weights$log.weights))
    fold.eval.log.weight = stats::weighted.mean(curve.eval.log.weights,
                                                curve.weights)
    return (fold.eval.log.weight)
  })
  return (mean(fold.eval.log.weights))
}

## todo ii.match shrinkage lambda, decay factor
## todo rcpp cv data logscore
## todo model logscore dist
