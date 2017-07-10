## author_header begin
## Copyright (C) 2016 Logan C. Brooks, Sangwon Hyun
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

##' @include utils.R
##' @include match.R
NULL

##' bw.SJ with bw.nrd0 as fallback
##'
##' At least sometimes when calling bw.SJ on backfill updates
##' (lag.info$residual), an error is generated
##' ("sample is too sparse to find TD"). bw.SJnrd0 uses bw.SJ if it succeeds,
##' and falls back to bw.nrd0 if it generates any error.
##'
##' @param x numeric vector: the observations
##'
##' @return single numeric: the bandwidth selection
bw.SJnrd0 = function(x) {
  tryCatch(stats::bw.SJ(x),
           error=function(e) stats::bw.nrd0(x))
}

##' Time-parameterized kernel density estimation sim method, Markovian version
##'
##' Function for making forecasts with the basic time-parameterized kernel
##' density estimation method. This method estimates \code{diff(new.dat)[t-1]}
##' (used to produce \code{dat[t]}) for a trajectory \code{new.dat} based on
##' weighted kernel density estimation using values of \code{dat} at the
##' corresponding time of season (a Markov process). The weights are based on
##' \code{new.dat[t-1]} and the corresponding values in \code{dat}; the
##' weighting function is a Gaussian kernel with width determined by
##' \code{bw.SJnrd0}.
##'
##' @param dat a list of numeric vectors, one per past season, containing
##'   historical trajectories.
##' @param new.dat.sim a numeric vector (trajectory), numeric matrix (cbound
##'   trajectories), or sim object (list with $ys a numeric matrix (cbound
##'   trajectories) and $weights a numeric vector (associated weights)), with
##'   \code{NA}'s for all future or missing data points to forecast or infer;
##'   currently only supports \code{NA}'s at future points, not mixed in between
##'   non-\code{NA} data
##' @param baseline a "baseline level" for this dataset; roughly speaking, data
##'   below this level does not grow like an epidemic; currently ignored, but
##'   can be used as the \code{y.scale.baseline} by passing it through the
##'   \code{control.list} argument.
##' @param n.sims single non-\code{NA} integer value or \code{NULL}: the number
##'   of curves to sample from the inferred distribution, or \code{NULL} to
##'   match the number of trajectories in \code{new.dat.sim}
##'
##' @return a sim object (list with two components: \itemize{
##' \item{\code{ys}: }{a numeric matrix; each column is a different possible
##' trajectory for the current season, with \code{NA}'s in \code{new.dat.sim}
##' filled in with random draws from the forecasted distribution, and non-NA's
##' (observed data) filled in with an imagined resampling of noise based on the
##' model.}
##' \item{\code{weights}: }{a numeric vector; assigns a weight to each column of
##' \code{ys}, which is used by other methods that rely on importance
##' sampling).}
##' }
##'
##' @examples
##' fluview.nat.recent.df =
##'    trimPartialPastSeasons(fetchEpidataDF("fluview", "nat",
##'                           first.week.of.season=21L,
##'                           cache.file.prefix="fluview_nat_allfetch"),
##'            "wili", min.points.in.season=52L)
##' ## Recent historical seasons + current season, minus 2009 (nonseasonal
##' ## pandemic) season:
##' full.dat = split(fluview.nat.recent.df$wili, fluview.nat.recent.df$season)
##' names(full.dat) <- sprintf("S%s", names(full.dat))
##' full.dat <- full.dat[names(full.dat)!="S2009"]
##' ## Recent historical seasons minus 2009:
##' dat = head(full.dat, -1L)
##' ## Current season:
##' new.dat = tail(full.dat, 1L)[[1]]
##' ## Sample from conditional curve distribution estimate using CDC's 2015
##' ## national %wILI onset threshold baseline of 2.1:
##' sim = twkde.markovian.sim(dat, new.dat, 2.1, n.sims=100)
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
twkde.markovian.sim = function(dat, new.dat.sim, baseline=NA_real_, n.sims=2000) {
  dat <- match.dat(dat)
  new.dat.sim <- match.new.dat.sim(new.dat.sim)
  baseline <- match.single.na.or.numeric(baseline) # (ignored by twkde though)
  if (is.null(n.sims)) {
    n.sims <- ncol(new.dat.sim$ys)
  }
  if (ncol(new.dat.sim$ys)==n.sims) {
    ys = new.dat.sim$ys
    sim.weights = new.dat.sim$weights
  } else {
    inds = sample(seq_along(new.dat.sim$weights), n.sims, prob=new.dat.sim$weights, replace=TRUE)
    ys = new.dat.sim$ys[,inds]
    sim.weights = rep(1, n.sims)
  }
  min.n.out = min(sapply(dat, length))
  obs.mat = sapply(dat, `[`, seq_len(min.n.out))
  obs.bws = sapply(seq_len(min.n.out), function(time.of.obs)
    bw.SJnrd0(obs.mat[time.of.obs,]))
  diff.obs.mat = diff(obs.mat)
  diff.obs.bws = sapply(seq_len(min.n.out-1), function(time.of.obs1m)
    bw.SJnrd0(diff.obs.mat[time.of.obs1m,]))
  for (time.of.obs in seq_len(nrow(ys))) {
    if (time.of.obs > min.n.out) {
      model.time.of.obs = min.n.out
    } else {
      model.time.of.obs = time.of.obs
    }
    if (!is.na(ys[time.of.obs,1])) {
    } else {
      for (sim.i in seq_len(n.sims)) {
        last.obs.log.weight = dnorm(ys[time.of.obs-1,sim.i],
                                    obs.mat[model.time.of.obs-1,],
                                    obs.bws[model.time.of.obs-1], log=TRUE)
        log.weights = last.obs.log.weight
        weights = exp(log.weights-max(log.weights))
        ys[time.of.obs,sim.i] <- ys[time.of.obs-1,sim.i] +
          sample(diff.obs.mat[model.time.of.obs-1,], 1,
                 prob=weights,
                 replace=TRUE) +
          rnorm(1, 0, diff.obs.bws[model.time.of.obs-1])
      }
    }
  }
  ## weights = rep(1, n.sims)
  ## sim = list(ys=ys, weights=weights)
  sim = list(ys=ys, weights=sim.weights)
  return (sim)
}

##' Time-parameterized kernel density estimation sim method with heuristic adjustments
##'
##' Function for making forecasts with the time-parameterized kernel density
##' estimation method with tweaks. This method estimates \code{diff(new.dat)}
##' based on weighted kernel density estimation. Weights are based on the time
##' of season and four functions of \code{new.dat}: (a) the last observed value
##' in \code{new.dat}, (b) the sum of observed values in \code{new.dat}, (c) an
##' exponential moving average of the observed values in \code{new.dat}, and (d)
##' an exponential moving average of the changes in observed values in
##' \code{new.dat} (i.e., in \code{diff(new.dat)}). The weighting function is
##' separable, and consists of two components: a highly weighted "base"
##' weighting function and a lowly weighted boxcar weighting function. The base
##' weighting function is the product of an integral Laplacian kernel with
##' respect to time of season, and Gaussian kernels with respect to the four
##' \code{new.dat}-based components (with bandwidths selected by the
##' \code{\link{bw.SJnrd0}} method, and relative weighting controlled by
##' \code{tradeoff.weights}). Each time a difference is drawn, simulating
##' \code{diff(new.dat)[t-1]}, the corresponding result for \code{new.dat[t]} is
##' linearly mixed with a randomly selected value from historical curves around
##' that time.
##'
##' @param dat a list of numeric vectors, one per past season, containing
##'   historical trajectories.
##' @param new.dat.sim a numeric vector (trajectory), numeric matrix (cbound
##'   trajectories), or sim object (list with $ys a numeric matrix (cbound
##'   trajectories) and $weights a numeric vector (associated weights)), with
##'   \code{NA}'s for all future or missing data points to forecast or infer;
##'   currently only supports \code{NA}'s at future points, not mixed in between
##'   non-\code{NA} data
##' @param baseline a "baseline level" for this dataset; roughly speaking, data
##'   below this level does not grow like an epidemic; currently ignored, but
##'   can be used as the \code{y.scale.baseline} by passing it through the
##'   \code{control.list} argument.
##' @param n.sims single non-\code{NA} integer value or \code{NULL}: the number
##'   of curves to sample from the inferred distribution, or \code{NULL} to
##'   match the number of trajectories in \code{new.dat.sim}
##' @param decay.factor decay factor for the exponential moving average of
##'   covariate.
##' @param diff.decay.factor decay factor for the exponential moving average of
##'   differences covariate.
##' @param max.shifts numeric vector with length matching the trajectory length
##'   in \code{new.dat.sim}; specifies the width of the time-of-season kernel as
##'   a function of the time of season.
##' @param shift.decay.factor decay factor for the time-of-season Laplacian
##'   kernel component.
##' @param tradeoff.weights log-scale weighting factors for the four
##'   non-time-based kernel components (last observed value, sum of observed
##'   values, exponential moving average of values, exponential moving average
##'   of differences).
##'
##' @return a sim object (list with two components: \itemize{
##' \item{\code{ys}: }{a numeric matrix; each column is a different possible
##' trajectory for the current season, with \code{NA}'s in \code{new.dat.sim} filled in with
##' random draws from the forecasted distribution, and non-\code{NA}'s (observed data)
##' filled in with an imagined resampling of noise based on the model.}
##' \item{\code{weights}: }{a numeric vector; assigns a weight to each column of
##' \code{ys}, which is used by other methods that rely on importance
##' sampling).}
##' }
##'
##' @examples
##' fluview.nat.recent.df =
##'    trimPartialPastSeasons(fetchEpidataDF("fluview", "nat",
##'                           first.week.of.season=21L,
##'                           cache.file.prefix="fluview_nat_allfetch"),
##'            "wili", min.points.in.season=52L)
##' ## Recent historical seasons + current season, minus 2009 (nonseasonal
##' ## pandemic) season:
##' full.dat = split(fluview.nat.recent.df$wili, fluview.nat.recent.df$season)
##' names(full.dat) <- sprintf("S%s", names(full.dat))
##' full.dat <- full.dat[names(full.dat)!="S2009"]
##' ## Recent historical seasons minus 2009:
##' dat = head(full.dat, -1L)
##' ## Current season:
##' new.dat = tail(full.dat, 1L)[[1]]
##' ## Sample from conditional curve distribution estimate using CDC's 2015
##' ## national %wILI onset threshold baseline of 2.1:
##' sim = twkde.sim(dat, new.dat, 2.1, n.sims=50)
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
twkde.sim = function(## dat, new.dat.sim
                     full.dat,
                     baseline=NA_real_,
                     n.sims=2000,
                     decay.factor=0.7,
                     diff.decay.factor=0.5,
                     max.shifts=c(rep(10,20),10:1,rep(0,3),1:10,rep(10,10)),
                     shift.decay.factor = 0.7,
                     tradeoff.weights=c(0.5, 0.25, 0.25, 0.5)) {

  ## extract historical data and future data from full.dat
  dat = head(full.dat, -1L)
  dat <- match.dat(dat)
  new.dat.sim = tail(full.dat, 1L)[[1]]
  new.dat.sim <- match.new.dat.sim(new.dat.sim)
  old.season.labels = head(names(full.dat), -1L)
  new.season.label = tail(names(full.dat), 1L)

  baseline <- match.single.na.or.numeric(baseline) # (ignored by twkde though)
  n.sims <- match.single.nonna.integer.or.null(n.sims)


  if (is.null(n.sims)) {
    n.sims <- ncol(new.dat.sim$ys)
  }
  if (ncol(new.dat.sim$ys)==n.sims) {
    ys = new.dat.sim$ys
    sim.weights = new.dat.sim$weights
  } else {
    inds = sample(seq_along(new.dat.sim$weights), n.sims, prob=new.dat.sim$weights, replace=TRUE)
    ys = new.dat.sim$ys[,inds]
    sim.weights = rep(1, n.sims)
  }
  orig.min.n.out = min(sapply(dat, length))
  ## orig.obs.mat = sapply(dat, `[`, seq_len(orig.min.n.out)) # xxx use dat.to.matrix
  orig.obs.mat = dat.to.matrix(dat, orig.min.n.out)
  for (time.of.obs in seq_len(nrow(new.dat.sim$ys))) {
    if (!is.na(new.dat.sim$ys[time.of.obs,1])) {
    } else {
      ## Include data from other times with discounted weights; do not consider
      ## any data from more than =max.shifts[time.of.obs]= weeks away.
      max.shift = max.shifts[time.of.obs]
      ## "Model times" are the times for the historical data we are using; they
      ## should be no more than =max.shift= away from =time.of.obs=, >=2 so we
      ## can calculate a delta, and <= the season length.
      min.model.time = max(2L, time.of.obs-max.shift)
      max.model.time = min(orig.min.n.out, time.of.obs+max.shift)
      ## Check that we will have any valid model times (one example where we
      ## don't is when predicting week 53 with max.shift=0 and no other seasons
      ## in the history with 53 weeks):
      if (min.model.time > max.model.time) stop("not enough shifting to have valid points")
      model.times = min.model.time:max.model.time
      ## Make a new obs.mat like the one used in the Markovian simulation but
      ## containing shifted versions of itself according to min.model.time and
      ## max.model.time.
      obs.mat = do.call(cbind, lapply(model.times, function(model.time) {
        ## We want all histories to have the same length so "area under curve",
        ## etc., are comparable. Select =min.model.time=-length segments of each
        ## historical curves, with time shift model.time-min.model.time so that
        ## the max time included is =model.time=.
        orig.obs.mat[seq_len(min.model.time)+model.time-min.model.time,]
      }))
      shift.decay = rep(shift.decay.factor^(abs(model.times-time.of.obs)), each=ncol(orig.obs.mat))
      min.n.out = nrow(obs.mat)
      ## Since we selected =min.model.time=-length segments from the history to
      ## populate =obs.mat=, it's as if we're predicting at =min.model.time=
      ## when looking at =obs.mat=, regardless of the actual =model.time=.
      model.time.of.obs = min.model.time
      obs.bws = sapply(seq_len(min.n.out), function(time.of.obs)
        bw.SJnrd0(obs.mat[time.of.obs,]))
      diff.obs.mat = diff(obs.mat)
      diff.obs.bws = sapply(seq_len(min.n.out-1), function(time.of.obs1m)
        bw.SJnrd0(diff.obs.mat[time.of.obs1m,]))
      last.obs.hist = obs.mat[model.time.of.obs-1,]
      last.obs.bw = obs.bws[model.time.of.obs-1]
      sum.obs.hist = colSums(obs.mat[seq_len(model.time.of.obs-1),,drop=FALSE])
      sum.obs.bw = sum(obs.bws[seq_len(model.time.of.obs-1)])
      ## sum.obs.bw = bw.SJnrd0(sum.obs.hist)
      decay.pattern = rev(decay.factor^seq_len(model.time.of.obs-1))
      decay.obs.hist = colSums(decay.pattern*obs.mat[seq_len(model.time.of.obs-1),,drop=FALSE])
      decay.obs.bw = bw.SJnrd0(decay.obs.hist)
      diff.decay.pattern = rev(diff.decay.factor^seq_len(max(0,model.time.of.obs-2)))
      decay.diff.obs.hist = colSums(diff.decay.pattern*diff.obs.mat[seq_len(max(0,model.time.of.obs-2)),,drop=FALSE])
      decay.diff.obs.bw = bw.SJnrd0(decay.diff.obs.hist)
      ## xxx versus recent differences?  check for similarity in slope?
      ## todo use some weighted kernel density bandwidth selection method
      ## instead of the unweighted bw.SJ/bw.nrd0 rules.
      for (sim.i in seq_len(n.sims)) {
        last.obs = ys[time.of.obs-1, sim.i]
        sum.obs = sum(ys[seq_len(model.time.of.obs-1)+time.of.obs-model.time.of.obs,sim.i])
        decay.obs = sum(decay.pattern*ys[seq_len(model.time.of.obs-1)+time.of.obs-model.time.of.obs,sim.i])
        diff.y = diff(ys[seq_len(model.time.of.obs-1)+time.of.obs-model.time.of.obs,sim.i])
        decay.diff.obs = sum(diff.decay.pattern*diff.y)
        last.obs.log.weights = dnorm(last.obs, last.obs.hist, last.obs.bw, log=TRUE)
        sum.obs.log.weights = dnorm(sum.obs, sum.obs.hist, sum.obs.bw, log=TRUE)
        decay.obs.log.weights = dnorm(decay.obs, decay.obs.hist, decay.obs.bw, log=TRUE)
        decay.diff.obs.log.weights = dnorm(decay.diff.obs, decay.diff.obs.hist, decay.diff.obs.bw, log=TRUE)
        log.weights =
          tradeoff.weights[1] * last.obs.log.weights +
          tradeoff.weights[2] * sum.obs.log.weights +
          tradeoff.weights[3] * decay.obs.log.weights +
          tradeoff.weights[4] * decay.diff.obs.log.weights
        log.weights <- shift.decay * log.weights
        weights = exp(log.weights-max(log.weights))
        ## print(max(weights)/sum(weights))
        weights <- weights/sum(weights)
        weights <- 0.9*weights+0.1/length(weights)
        index = sample.int(ncol(diff.obs.mat), 1, prob=weights, replace=TRUE)
        ## ys[time.of.obs,sim.i] <-
        ##     ys[time.of.obs-1,sim.i] + diff.obs.mat[model.time.of.obs-1,index] +
        ##     rnorm(1, 0, diff.obs.bws[model.time.of.obs-1])
        y.shrink = 0.9
        ## ys[time.of.obs,sim.i] <-
        ##     y.shrink*(ys[time.of.obs-1,sim.i] + diff.obs.mat[model.time.of.obs-1,index] +
        ##               rnorm(1, 0, diff.obs.bws[model.time.of.obs-1])
        ##     ) +
        ##     (1-y.shrink)*(obs.mat[model.time.of.obs,index] +
        ##                   rnorm(1, 0, obs.bws[model.time.of.obs])
        ##     )
        index2 = sample.int(ncol(obs.mat), 1)
        ys[time.of.obs,sim.i] <-
          y.shrink*(ys[time.of.obs-1,sim.i] +
                    diff.obs.mat[model.time.of.obs-1,index] +
                    rnorm(1, 0, diff.obs.bws[model.time.of.obs-1])
          ) +
          (1-y.shrink)*(obs.mat[model.time.of.obs,index2] +
                        rnorm(1, 0, obs.bws[model.time.of.obs])
          )
      }
    }
  }
  ## ys <- exp(ys) - 0.01
  ## weights = rep(1, n.sims)
  ## sim = list(ys=ys, weights=weights)

  ## Make a dummy control list, containing only model name
  control.list = list(model = "twkde")

  ## Return sim object
    sim = list(ys=ys,
               weights=sim.weights,
               control.list=control.list,
               old.dat = list(dat),
               new.dat = as.numeric(unlist(new.dat.sim[["ys"]])),
               old.season.labels = list(old.season.labels),
               new.season.label = list(new.season.label))
    class(sim) <- "sim"
    return (sim)
    }

## todo fix twkde =max.shifts= time smearing default to be generic, not flu-specific
## fixme do twkde methods actually fill in non_NA's with estimated noise resampling? don't think so, fix doc
