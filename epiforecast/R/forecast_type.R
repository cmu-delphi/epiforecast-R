## author_header begin
## Copyright (C) 2017 Logan C. Brooks
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
##'
##' @seealso weighted.bw.nrd0
##'
##' @export
weighted.bw.nrd0ish = function(x, w) {
  w <- match.nonnegative.numeric(w)
  ## fac =
  fac.a =
    min(diff(unname(Hmisc::wtd.quantile(x, w, c(0.25, 0.75))))/1.34,
        dplyr::coalesce(sqrt(Hmisc::wtd.var(x, w)), 0))
  ## if (fac==0) fac <- weighted.mean(abs(x), w)
  ## if (fac==0) fac <- 1
  fac.b = weighted.mean(abs(x), w)
  if (fac.b==0) fac.b <- 1
  fac <- weighted.mean(
    c(fac.a, fac.b),
    c(sum(w), 3)
  )
  result = 0.9 * fac * sum(w)^-0.2
  return (result)
}

get_na_value_or_empty_for_target = function(target.info, ...) {
  bin.info = target.info[["bin_info_for"]](...)
  breaks = bin.info[["breaks"]]
  if (bin.info[["include.na"]]) {
    ## get NA of appropriate type:
    breaks[NA_integer_][[1L]]
  } else {
    breaks[integer(0)]
  }
}

get_na_string_or_empty_for_target = function(target.info, ...) {
  target.info[["unit"]][["to_string"]](
    get_na_value_or_empty_for_target(target.info, ...),
    ...
  )
}

point.mae.forecast.type = list(
  Type = "Point",
  Value_from_weighted_univals = function(target.info,
                                         target.values, target.weights,
                                         ...) {
    matrixStats::weightedMedian(target.values, target.weights, na.rm=TRUE)
  },
  Bin_start_incl_for = function(target.info, ...) {
    return (NA_character_)
  },
  Bin_end_notincl_for = function(target.info, ...) {
    return (NA_character_)
  }
)

distr.logscore.forecast.type = list(
  Type = "Bin",
  Value_from_weighted_univals = function(target.info,
                                         target.values, target.weights,
                                         label.bins=TRUE,
                                         uniform.pseudoweight.total=3,
                                         smooth.sim.targets=TRUE,
                                         ...) {
    target.weights <- match.nonnegative.numeric(target.weights)
    ## get bin info and check:
    bin.info = target.info[["bin_info_for"]](...)
    stopifnot(bin.info[["include.na"]] || !any(is.na(target.values)))
    breaks = bin.info[["breaks"]]
    rightmost.closed = bin.info[["rightmost.closed"]]
    include.na = bin.info[["include.na"]]
    if (length(breaks) <= 1L || any(is.na(breaks))) {
      stop ("Must have at least 2 breaks, and all breaks must be non-NA.")
    }
    if (!all(diff(breaks) > 0)) {
      stop ("The break values must all be unique, and in strictly increasing order.")
    }
    if (!all(is.na(target.values) |
             target.values >= breaks[[1L]] &
             (if (rightmost.closed) target.values <= breaks[[length(breaks)]]
              else target.values < breaks[[length(breaks)]])
             )) {
      stop ("All non-NA target.values must lie between the first break and the last break.")
    }
    ## add uniform pseudoweight (if nonzero):
    uniform.pseudoweight.total <- match.single.nonna.numeric(uniform.pseudoweight.total)
    uniform.pseudoweight.total <- match.nonnegative.numeric(uniform.pseudoweight.total)
    break.bin.representatives = bin.info[["break.bin.representatives"]]
    if (is.null(break.bin.representatives)) {
      stop ("break.bin.representatives missing in returned bin.info")
    }
    if (uniform.pseudoweight.total > 0) {
      pseudo.target.values = c(
        ## break bin representatives:
        break.bin.representatives,
        ## NA bin representative if NA's are included; else 0-length vector:
        breaks[NA_integer_][include.na]
      )
      pseudo.target.weights = rep(
        uniform.pseudoweight.total/length(pseudo.target.values),
        length(pseudo.target.values)
      )
      target.values <- c(target.values, pseudo.target.values)
      target.weights <- c(target.weights, pseudo.target.weights)
    }
    ## split NA and non-NA target values (and corresponding weights):
    target.value.is.na = is.na(target.values)
    n.non.na.bins = length(breaks)-1L
    n.na.bins = as.integer(include.na)
    n.bins = n.non.na.bins + n.na.bins
    ## (tiny pre-optimization --- collapse all-FALSE vector to single FALSE):
    if (!include.na) target.value.is.na <- FALSE
    ## smooth or tabulate weights:
    if (smooth.sim.targets) {
      ## separate out non-NA target values and weights:
      na.bin.weight = sum(target.weights[target.value.is.na])
      non.na.target.values = target.values[!target.value.is.na]
      non.na.target.weights = target.weights[!target.value.is.na]
      ## bandwidth estimate for weighted data (density by default ignores the
      ## weights in the bandwidth estimate):
      bw = weighted.bw.nrd0ish(non.na.target.values, non.na.target.weights)
      ## kernel density estimate at function-determined points (=ks= package has
      ## a nice interface, but uses =KernSmooth::dpik= by default for bandwidth
      ## selection, which, like the default for =density=, ignores weights, and
      ## needs settings specified to avoid issues when a lot of target weight is
      ## highly concentrated on a single target value; and produces segfaults):
      pdf.fit = stats::density(non.na.target.values,
                               bw=bw,
                               weights=non.na.target.weights %>>%
                                 magrittr::divide_by(sum(.)))
      .GlobalEnv[["g.pdf.fit"]] <- pdf.fit
      cdf.fit = list(x=pdf.fit[["x"]], y=cumsum(pdf.fit[["y"]]))
      ## linear interpolation, constant extrapolation (so 0 mass assigned
      ## outside range(cdf.fit[["x"]])...):
      cdf.fit.at.breaks = cdf.fit %>>% stats::approx(xout=breaks, rule=2L)
      non.na.bin.weights =
        diff(cdf.fit.at.breaks[["y"]]) %>>%
        ## the input to the density estimate was based on a scale neglecting
        ## na.bin.weight, and the real line as the range; for now, just rescale
        ## all non-NA bin weights proportionally to leave room for na.bin.weight
        ## and incorporate any mass falling outside the range of the breaks:
        magrittr::multiply_by(sum(non.na.target.weights)/sum(.))
      bin.weights = c(
        non.na.bin.weights,
        if(include.na) na.bin.weight else numeric(0)
      )
    } else {
      bin = dplyr::coalesce(
                     findInterval(target.values, breaks,
                                  rightmost.closed=rightmost.closed,
                                  all.inside=TRUE),
                     length(breaks))
      bin.weights = weighted_tabulate(bin, n.bins, target.weights)
    }
    if (label.bins) {
      left.delimiter = "["
      left.break.string = target.info[["unit"]][["to_string"]](
        breaks[-length(breaks)], ...
      )
      right.break.string = target.info[["unit"]][["to_string"]](
        breaks[-1L], ...
      )
      right.delimiter =
        if (rightmost.closed) c(rep(")", length(breaks)-2L), "]")
        else ")"
      bin.names = c(
        paste0(left.delimiter, left.break.string, ", ",
               right.break.string, right.delimiter),
        get_na_string_or_empty_for_target(target.info, ...)
        )
      names(bin.weights) <- bin.names
    }
    ## normalize bin weights:
    bin.weights <- bin.weights/sum(bin.weights)
    return (bin.weights)
  },
  Bin_start_incl_for = function(target.info, ...) {
    bin.info = target.info[["bin_info_for"]](...)
    breaks = bin.info[["breaks"]]
    Bin_start_incl = c(
      target.info[["unit"]][["to_string"]](
        breaks[-length(breaks)], ...
      ),
      get_na_string_or_empty_for_target(target.info, ...)
    )
    return (Bin_start_incl)
  },
  Bin_end_notincl_for = function(target.info, ...) {
    bin.info = target.info[["bin_info_for"]](...)
    breaks = bin.info[["breaks"]]
    Bin_end_notincl = c(
      target.info[["unit"]][["to_string"]](
        breaks[-1L], ...
      ),
      get_na_string_or_empty_for_target(target.info, ...)
    )
    return (Bin_end_notincl)
  }
)

partial_spreadsheet_from_weighted_univals =
  function(forecast.type, target.info,
           target.values, target.weights,
           only.Value=FALSE,
           ...) {
    Value = forecast.type[["Value_from_weighted_univals"]](target.info, target.values, target.weights, label.bins=FALSE, ...)
    if (only.Value) {
      return (tibble::tibble(Value=Value))
    } else {
      partial.spreadsheet =
        tibble::tibble(
                  Target=target.info[["Target"]],
                  Type=forecast.type[["Type"]],
                  Unit=target.info[["unit"]][["Unit"]],
                  Bin_start_incl=forecast.type[["Bin_start_incl_for"]](target.info, ...),
                  Bin_end_notincl=forecast.type[["Bin_end_notincl_for"]](target.info, ...),
                  Value=Value
                )
      return (partial.spreadsheet)
    }
}

flusight2016.proxy.forecast.types = list(
  point.mae.forecast.type,
  distr.logscore.forecast.type
)
