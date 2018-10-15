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


##' Creates a uniform distribution over discrete choices which can be used with \code{\link{get.eb.control.list}}.
##'
##' @param choices a vector of (discrete) choices
##'
##' @return a uniform discrete distribution over \code{choices}.
##'
##' @examples
##' 
##' uniform.seq = unifChoicePrior(letters[1:5])
##' 
##' ## The distributions used by EB can be broken down into buckets;
##' ## for the uniform discrete distribution, each bucket corresponds
##' ## (boringly) to a single choice from =choices=. However, it is
##' ## important to have a common interface.
##' random.bucket.indices = sample(seq_len(uniform.seq$n), 10000, replace=TRUE, prob=uniform.seq$probs)
##' random.elements = uniform.seq$sampler(random.bucket.indices)
##' random.elements.another.way = uniform.seq$choices[random.bucket.indices] # only works for =unifChoicePrior=
##' random.elements.a.third.way = letters[random.bucket.indices] # only works for this example
##'
##' @export
unifChoicePrior = function(choices) {
  probs = rep(1.0/length(choices), length(choices)) # uniform probabilities
  sampler = function(inds) choices[inds]
  return (list(n = length(choices), choices = choices, probs = probs, sampler = sampler))
}
## xxx refactor to discretePrior with probs (default to uniform) and possibly tapply-like groups (determining bucket grouping, default to seq_along(choices))

integerTriangleAroundZeroPrior = function(max.val) {
  val.ramp = Seq(1,max.val)
  choices = c(rev(-val.ramp), 0, val.ramp)
  weights = c(val.ramp, max.val+1, rev(val.ramp))
  probs = weights/sum(weights)
  sampler = function(inds) choices[inds]
  return (list(n = length(choices), choices = choices, probs = probs, sampler = sampler))
}

unifDistrUnifGridPrior = function(min.val, max.val, n.choices) {
  # We want to discretize the uniform distribution U[min.val, max.val].
  # Divide interval [min.val, max.val] into n.choices uniform-width buckets.
  # This code chooses the midpoints of these buckets as the representative elements (choices):
  choices = (1:n.choices - 0.5)/n.choices * (max.val-min.val) + min.val
  # Uniform probabilities for uniform grid:
  probs = rep(1.0/n.choices, n.choices)
  # Note for future: for other grids, actually compute the bucket boundaries and endpoints, then take differences of cdf
  mins = ((1:n.choices)-1)/n.choices * (max.val-min.val) + min.val
  maxes = (1:n.choices)/n.choices * (max.val-min.val) + min.val
  sampler = function(inds) stats::runif(length(inds), min=mins[inds], max=maxes[inds])
  return (list(n = n.choices, choices = choices, probs = probs, sampler=sampler))
}

logUnifDistrLogUnifGridPrior = function(min.val, max.val, n.choices) {
  # We want to discretize the uniform distribution U[min.val, max.val].
  # Divide interval [min.val, max.val] into n.choices uniform-width buckets.
  # This code chooses the midpoints of these buckets as the representative elements (choices):
    min.val <- log(min.val)
    max.val <- log(max.val)
  choices = (1:n.choices - 0.5)/n.choices * (max.val-min.val) + min.val
  # Uniform probabilities for uniform grid:
  probs = rep(1.0/n.choices, n.choices)
  # Note for future: for other grids, actually compute the bucket boundaries and endpoints, then take differences of cdf
  mins = ((1:n.choices)-1)/n.choices * (max.val-min.val) + min.val
  maxes = (1:n.choices)/n.choices * (max.val-min.val) + min.val
  sampler = function(inds) exp(stats::runif(length(inds), min=mins[inds], max=maxes[inds]))
  return (list(n = n.choices, choices = exp(choices), probs = probs, sampler=sampler))
}

unifDistrGaussianGridPrior = function(min.val, max.val, n.choices, mean=(max.val+min.val)/2, sd=(max.val-min.val)/4) {
  # Non-truncated normal CDF at window ends:
  left.cd = stats::pnorm(min.val, mean=mean, sd=sd)
  right.cd = stats::pnorm(max.val, mean=mean, sd=sd)
  
  endpoint.cds = (0:n.choices)/n.choices * (right.cd-left.cd) + left.cd
  endpoints = stats::qnorm(endpoint.cds,mean=mean,sd=sd)
  if(endpoints[1] + 1e-12 < min.val) stop("Bug in tnormcdf reasoning.")
  if(endpoints[length(endpoints)] - 1e-12 > max.val) stop("Bug in tnormcdf reasoning.")
  # Domain- and mass- midpoints of buckets are choosen as the choices (these are the same since we are dealing with a uniform distribution):
  choices = (endpoints[1:n.choices] + endpoints[2:(n.choices+1)])/2
  probs = diff(endpoints)/(max.val-min.val)

  # todo sampler

  return (list(n=n.choices,choices=choices,probs=probs))
}

unifLocGridPrior = function(fit) {
  locs = max.col(t(fit$f), ties.method="last")
  #locs = apply(fit$f,2,which.max)
  delta = mean(diff(sort(locs)))/2
  choices = Seq(max(round(min(locs)-delta),1),
                min(round(max(locs)+delta),nrow(fit$f)))
  return (unifChoicePrior(choices))
}

# Uniform grid spacing, uniform scale (underlying peak height) distribution.
unifScaleUnifGridPrior = function(fit, n.choices) {
  peaks = apply(fit$f,2,max)
  delta = mean(diff(sort(peaks)))/2
  
  min.val = min(peaks)-delta
  max.val = max(peaks)+delta
  
  return (unifDistrUnifGridPrior(min.val, max.val, n.choices))
}

# Uniform grid spacing, uniform scale (underlying peak height) distribution.
unifScaleGaussianGridPrior = function(fit, n.choices) {
  peaks = apply(fit$f,2,max)
  delta = mean(diff(sort(peaks)))/2
  
  min.val = min(peaks)-delta
  max.val = max(peaks)+delta
  
  return (unifDistrGaussianGridPrior(min.val, max.val, n.choices))
}
