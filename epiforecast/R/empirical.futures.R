## author_header begin
## Copyright (C) 2017 Logan C. Brooks, Sangwon Hyun
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

##' Simulate future in current trajectory with empirical (historical)
##' distribution
##'
##' @template sim.method_template
##'
##' @examples
##' ## National-level ILINet weighted %ILI data for recent seasons, excluding 2009 pandemic:
##' area.name = "nat"
##' full.dat = fetchEpidataFullDat("fluview", area.name, "wili",
##'                                min.points.in.season = 52L,
##'                                first.week.of.season = 31L,
##'                                cache.file.prefix=sprintf("fluview_%s_fetch", area.name))
##' full.dat <- full.dat[names(full.dat)!="S2009"]
##' ## Sample from conditional curve distribution estimate using the above data
##' ## and CDC's 2015 national %wILI onset threshold baseline of 2.1:
##' sim = empirical.futures.sim(full.dat, baseline=2.1, max.n.sims=100)
##' print(sim)
##' plot(sim, type="lineplot")
##'
##' @author Logan C. Brooks, David C. Farrow, Sangwon Hyun, Ryan J. Tibshirani, Roni Rosenfeld
##'
##' @export
empirical.futures.sim = function(full.dat, baseline=NA_real_, max.n.sims=2000L) {
  ## extract historical data and future data from full.dat
  dat = head(full.dat, -1L)
  dat <- match.dat(dat)
  new.dat.sim = tail(full.dat, 1L)[[1]]
  new.dat.sim <- match.new.dat.sim(new.dat.sim)
  old.season.labels = head(names(full.dat), -1L)
  new.season.label = tail(names(full.dat), 1L)
  baseline <- match.single.na.or.numeric(baseline) # (ignored by twkde though)
  max.n.sims <- match.single.nonna.integer.or.null(max.n.sims)

  n.out = nrow(new.dat.sim[["ys"]])
  ## sim object containing specified "pasts":
  past.sim = new.dat.sim
  ## sim object containing empirical "futures":
  future.sim = list(
    ys = sapply(dat, function(historical.trajectory) {
      approx(historical.trajectory, xout=seq_len(nrow(new.dat.sim[["ys"]])), rule=2L)[["y"]]
    }),
    weights = rep(1, length(dat))
  )
  ## We want to build a sim object for the random variable dplyr::coalesce(PAST,
  ## FUTURE), treating PAST and FUTURE as independent. If the n.sims for the
  ## past and future sim objects are small enough, we can calculate the
  ## resulting output sim exactly; otherwise, we will need to use sampling.
  n.past.sims = length(past.sim[["weights"]])
  n.future.sims = length(future.sim[["weights"]])
  if (n.past.sims * n.future.sims <= max.n.sims) {
    ## We can calculate the output sim exactly. (This case will typically be
    ## triggered when new.dat.sim represents a single possible trajectory rather
    ## than a distribution.)

    ## Pair each new.dat.sim partial trajectory with each dat trajectory:
    past.sim.is = rep(seq_len(n.past.sims), each=n.future.sims)
    future.sim.is = rep(seq_len(n.future.sims), times=n.past.sims)
    ## Construct wider sim objects that together correspond to a sim object over
    ## the cartesian product <PAST, FUTURE>:
    product.past.sim = list(
      ys=past.sim[["ys"]][,past.sim.is],
      weights=past.sim[["weights"]][past.sim.is]
    )
    product.future.sim = list(
      ys=future.sim[["ys"]][,future.sim.is],
      weights=future.sim[["weights"]][future.sim.is]
    )
  } else {
    ## Construct (potentially) wider sim objects that together correspond to a
    ## sim object over the cartesian product <PAST, FUTURE>:
    product.past.sim = upsample_sim(past.sim, max.n.sims, TRUE)
    product.future.sim = upsample_sim(future.sim, max.n.sims, TRUE)
  }

  ys = dplyr::coalesce(product.past.sim[["ys"]], product.future.sim[["ys"]])
  weights = product.past.sim[["weights"]] * product.future.sim[["weights"]]

  ## Make a dummy control list, containing only model name
  control.list = list(model = "empirical.futures")

  ## Return sim object
  sim = list(ys=ys,
             weights=weights,
             control.list=control.list,
             old.dat = list(dat),
             ## fake a vector new.dat if necessary:
             new.dat = rowMeans(new.dat.sim[["ys"]]),
             old.season.labels = list(old.season.labels),
             new.season.label = list(new.season.label))
  class(sim) <- "sim"
  return (sim)
}
