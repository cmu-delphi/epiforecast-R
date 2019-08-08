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

##' @include ensemble.R
##' @import pipeR
NULL

##' @include match.R
##'
##' @seealso bw.nrd0
##'
##' @export
weighted.bw.nrd0ish = function(x, w=rep(1, length(x))) {
  w <- match.nonnegative.numeric(w)
  if (length(w) != length(x)) {
    stop ("w must have the same length as x")
  }
  fac.pos.pseudoweight = 3
  ## the core of the nrd type rules (taken along with the scaling at the end),
  ## potentially a problematic 0:
  fac.data = min(diff(unname(Hmisc::wtd.quantile(x, w, c(0.25, 0.75))))/1.34,
                 dplyr::coalesce(sqrt(Hmisc::wtd.var(x, w)), 0))
  ## a strictly positive scale fac based on the mean magnitude of x, unless it
  ## is too close to 0, in which case an arbitrary scale of 1 is assumed:
  fac.pos = weighted.mean(abs(x), w)
  if (fac.pos < .Machine[["double.eps"]]) fac.pos <- 1
  ## combine using w and pseudoweight:
  fac <- weighted.mean(
    c(fac.data, fac.pos),
    c(sum(w), fac.pos.pseudoweight)
  )
  ## nrd scaling rule, adjusted to incorporate the pseudoweight
  result = 0.9 * fac * (sum(w)+fac.pos.pseudoweight)^-0.2
  return (result)
}
## fixme this is slower than more complex (but unweighted) bandwidth selection rules, and seems to produce overly large bandwidth selections

get_na_value_or_empty_for_target = function(target.spec, ...) {
  bin.info = target.spec[["bin_info_for"]](...)
  breaks = bin.info[["breaks"]]
  if (bin.info[["include.na"]]) {
    ## get NA of appropriate type:
    breaks[NA_integer_][[1L]]
  } else {
    breaks[integer(0)]
  }
}

get_na_binlabelstart_or_empty_for_target = function(target.spec, ...) {
  target.spec[["unit"]][["to_binlabelstart"]](
    get_na_value_or_empty_for_target(target.spec, ...),
    ...
  )
}

point_prediction_error = function(forecast.value, observed.value) {
  ## calculate absolute error with some special treatment of NA's
  if (is.na(observed.value)) {
    NA_real_
  } else if (is.na(forecast.value)) {
    Inf
  } else {
    ## (note: all evaluate_forecast_value's should return numerics)
    return (as.numeric(forecast.value - observed.value))
  }
  ## xxx will probably also need the domain settings (e.g.,
  ## flusight2016.settings) to go with the target in some contexts, but only the
  ## target is passed in here
}

point.mae.forecast.type = list(
  Type = "Point",
  forecast_value_from_weighted_univals = function(target.spec,
                                                  target.values, target.weights,
                                                  ...) {
    matrixStats::weightedMedian(target.values, target.weights, na.rm=TRUE)
  },
  observed_value_from_observed_multival = function(target.spec,
                                                   observed.multival, ...) {
    ## xxx pretend that the median was observed instead; todo instead keep the
    ## multival and change the evaluation & maybe fitting functions
    median(observed.multival, na.rm=TRUE)
  },
  Value_from_forecast_value = function(forecast.value, target.spec, ...) {
    ## forecast.value is a point prediction in its natural representation; Value
    ## should be its numeric/character representation in a forecast spreadsheet;
    ## call the appropriate conversion method
    target.spec[["unit"]][["to_point"]](forecast.value, ...)
  },
  forecast_value_from_Value = function(Value, target.spec, ...) {
    target.spec[["unit"]][["from_point"]](Value, ...)
  },
  Bin_start_incl_for = function(target.spec, ...) {
    return (NA_character_)
  },
  Bin_end_notincl_for = function(target.spec, ...) {
    return (NA_character_)
  },
  evaluate_forecast_value = function(forecast.value, observed.value) {
    return (abs(point_prediction_error(forecast.value, observed.value)))
  },
  fit_ensemble_coefs = function(instance.method.forecast.values.listmat, instance.observed.values.list, total.instance.weight, fallback.method.index) {
    X = instance.method.forecast.values.listmat
    mode(X) <- "numeric"
    if (any(is.na(X[,fallback.method.index]))) {
      stop ("Fallback method must not contain any NA point predictions.")
    }
    y = instance.observed.values.list
    mode(y) <- "numeric"
    ## exclude instances where observation was NA:
    X <- X[!is.na(y),,drop=FALSE]
    y <- y[!is.na(y)]
    ## replace NA method point predictions with fallback values:
    X[is.na(X)] <- X[row(X)[is.na(X)], fallback.method.index]
    ## fit model and return coefficients:
    ## return (lasso_lad_coef(y, X, include.intercept=FALSE))
    return (simplex_lad_weights(y, X))
  }
)

unibin_log_score = function(forecast.value, observed.value) {
  ## unibin-like log score: if observed.value is not an indicator for a single
  ## value, treat it as a probability distribution over the "real" observed
  ## value, and calculate the expected unibin log score (KL divergence?):
  return (sum(observed.value[observed.value!=0]*
              log(forecast.value[observed.value!=0])))
}

multibin_log_score = function(forecast.value, observed.value) {
  ## (given different observed.value's than multibin_log_score: nearly always many-hot rather than nearly always one-hot)
  ##
  ## Sum all accepted mass and take log:
  log(sum(observed.value*forecast.value))
}

distr_from_weighted_univals =
  function(target.spec,
           target.values, target.weights,
           label.bins=TRUE,
           uniform.pseudoweight.total=3,
           smooth.sim.targets=TRUE,
           ...) {
    target.weights <- match.nonnegative.numeric(target.weights)
    ## get bin info and check:
    bin.info = target.spec[["bin_info_for"]](...)
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
      print(target.values)
      print(breaks)
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
      non.na.target.values.original = target.values[!target.value.is.na]
      non.na.target.values.to.smooth = target.spec[["unit"]][["shift_for_smoothing"]](non.na.target.values.original)
      non.na.target.weights = target.weights[!target.value.is.na]
      ## bandwidth estimate for weighted data (density by default ignores the
      ## weights in the bandwidth estimate):
      bw = weighted.bw.nrd0ish(non.na.target.values.to.smooth, non.na.target.weights)
      ## kernel density estimate at function-determined points (=ks= package has
      ## a nice interface, but uses =KernSmooth::dpik= by default for bandwidth
      ## selection, which, like the default for =density=, ignores weights, and
      ## needs settings specified to avoid issues when a lot of target weight is
      ## highly concentrated on a single target value; and produces segfaults):
      pdf.fit = stats::density(non.na.target.values.to.smooth,
                               bw=bw,
                               weights=non.na.target.weights %>>%
                                 magrittr::divide_by(sum(.)))
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
      left.break.string = target.spec[["unit"]][["to_binlabelstart"]](
        breaks[-length(breaks)], ...
      )
      right.break.string = target.spec[["unit"]][["to_binlabelend"]](
        breaks[-1L], ...
      )
      right.delimiter =
        if (rightmost.closed) c(rep(")", length(breaks)-2L), "]")
        else ")"
      bin.names = c(
        paste0(left.delimiter, left.break.string, ", ",
               right.break.string, right.delimiter),
        get_na_binlabelstart_or_empty_for_target(target.spec, ...)
      )
      names(bin.weights) <- bin.names
    }
    ## normalize bin weights:
    bin.weights <- bin.weights/sum(bin.weights)
    return (bin.weights)
  }

fit_distr_ensemble_coefs =
  function(instance.method.forecast.values.listmat, evaluate_forecast_value, instance.observed.values.list, total.instance.weight, fallback.method.index, excessively.low.fallback.log.score=-8) {
    ## fixme make consistent with updated metric --- weighting
    instance.method.log.scores.mat = Map(evaluate_forecast_value, instance.method.forecast.values.listmat, instance.observed.values.list) %>>%
      {
        dim(.) <- dim(instance.method.forecast.values.listmat)
        dimnames(.) <- dimnames(instance.method.forecast.values.listmat)
        mode(.) <- "numeric"
        .
      }
    if (any(instance.method.log.scores.mat[,fallback.method.index] <= excessively.low.fallback.log.score)) {
      stop ("Fallback method produced excessively low log score for at least one training instance.")
    }
    degenerate.em.weights = degenerate_em_weights(exp(instance.method.log.scores.mat))
    fallback.inflation.factor = min(1, 3/total.instance.weight)
    if (is.character(fallback.method.index)) {
      fallback.method.index <- which(colnames(instance.method.log.scores.mat)==fallback.method.index)
    }
    weights = (1-fallback.inflation.factor)*degenerate.em.weights +
      fallback.inflation.factor*as.numeric(seq_len(ncol(instance.method.log.scores.mat))==fallback.method.index)
    return (weights)
  }

distr.logscore.forecast.type = list(
  Type = "Bin",
  forecast_value_from_weighted_univals = distr_from_weighted_univals,
  observed_value_from_observed_multival = function(target.spec,
                                                   observed.multival,
                                                   label.bins=TRUE,
                                                   ...) {
    ## transform single observed val -> one-hot vector over bins
    ##
    ## xxx if there are multiple observed vals, produce a probability
    ## distribution over bins instead as if every val is equally likely
    later.args = list(...)
    ## set/override smoothing parameters:
    later.args[["uniform.pseudoweight.total"]] <- 0
    later.args[["smooth.sim.targets"]] <- FALSE
    do.call(distr_from_weighted_univals,
            c(list(
              target.spec,
              observed.multival,
              rep(1/length(observed.multival), length(observed.multival)),
              label.bins=label.bins
            ), later.args)
            )
  },
  Value_from_forecast_value = function(forecast.value, target.spec, ...) {
    forecast.value
  },
  forecast_value_from_Value = function(Value, target.spec, ...) {
    Value
  },
  Bin_start_incl_for = function(target.spec, ...) {
    bin.info = target.spec[["bin_info_for"]](...)
    breaks = bin.info[["breaks"]]
    Bin_start_incl = c(
      target.spec[["unit"]][["to_binlabelstart"]](
        breaks[-length(breaks)], ...
      ),
      get_na_binlabelstart_or_empty_for_target(target.spec, ...)
    )
    return (Bin_start_incl)
  },
  Bin_end_notincl_for = function(target.spec, ...) {
    bin.info = target.spec[["bin_info_for"]](...)
    breaks = bin.info[["breaks"]]
    Bin_end_notincl = c(
      target.spec[["unit"]][["to_binlabelend"]](
        breaks[-1L], ...
      ),
      get_na_binlabelstart_or_empty_for_target(target.spec, ...)
    )
    return (Bin_end_notincl)
  },
  evaluate_forecast_value = unibin_log_score,
  fit_ensemble_coefs = function(instance.method.forecast.values.listmat, instance.observed.values.list, total.instance.weight, fallback.method.index, excessively.low.fallback.log.score=-8) {
    fit_distr_ensemble_coefs(instance.method.forecast.values.listmat,
                             unibin_log_score,
                             instance.observed.values.list, total.instance.weight, fallback.method.index, excessively.low.fallback.log.score)
  }
)

multibin.logscore.forecast.type =
  distr.logscore.forecast.type[c(
    "Type",
    "Value_from_forecast_value",
    "forecast_value_from_Value",
    "Bin_start_incl_for",
    "Bin_end_notincl_for",
    "forecast_value_from_weighted_univals"
  )] %>>%
  c(list(
    observed_value_from_observed_multival =
      function(target.spec,
               observed.multival,
               label.bins=TRUE,
               ...) {
        unibin.observed.value = distr.logscore.forecast.type[["observed_value_from_observed_multival"]](
          target.spec, observed.multival, label.bins=label.bins, ...
        )
        ## flags for nonzero-weighted bins --- all are accepted centers for
        ## multibin score
        bin.is.valid.center = unibin.observed.value > 0
        ## expand around centers according to multibin neighbor matrix:
        bin.info = target.spec[["bin_info_for"]](...)
        multibin.neighbor.matrix = target.spec[["multibin_neighbor_matrix"]](bin.info, ...)
        bin.is.accepted = multibin.neighbor.matrix %*% bin.is.valid.center > 0
        ## flags to numeric:
        return (as.numeric(bin.is.accepted))
      },
    evaluate_forecast_value = multibin_log_score,
    ## forecast_value_from_weighted_univals = function(...) {
    ##   stop ("implementation decision todo: multibin forecast values --- reuse unibin forecast values vs. optimize")
    ## },
    fit_ensemble_coefs = function(instance.method.forecast.values.listmat, instance.observed.values.list, total.instance.weight, fallback.method.index, excessively.low.fallback.log.score=-8) {
      fit_distr_ensemble_coefs(instance.method.forecast.values.listmat,
                               multibin_log_score,
                               instance.observed.values.list, total.instance.weight, fallback.method.index, excessively.low.fallback.log.score)
    }
  ))

subspreadsheet_from_forecast_value = function(forecast.value,
                                              target.spec, forecast.type,
                                              only.Value=FALSE,
                                              ...) {
  Value = forecast.type[["Value_from_forecast_value"]](forecast.value, target.spec, ...)
  if (only.Value) {
    return (tibble::tibble(Value=Value))
  } else {
    subspreadsheet =
      ## tibble::tibble acts a bit slow in this case; forming and converting a
      ## data.frame instead is faster
      data.frame(
        Target=target.spec[["Target"]],
        Type=forecast.type[["Type"]],
        Unit=target.spec[["unit"]][["Unit"]],
        Bin_start_incl=forecast.type[["Bin_start_incl_for"]](target.spec, ...),
        Bin_end_notincl=forecast.type[["Bin_end_notincl_for"]](target.spec, ...),
        Value=Value,
        check.names=FALSE,
        stringsAsFactors=FALSE
      ) %>>%
      tibble::as_tibble()
    return (subspreadsheet)
  }
}

forecast_value_from_subspreadsheet =
  function(subspreadsheet, target.spec, forecast.type, ...) {
    ## tibble::tibble acts a bit slow in this case; forming and converting a
    ## data.frame instead is faster
    index.tbl =
      data.frame(
        Target=target.spec[["Target"]],
        Type=forecast.type[["Type"]],
        Unit=target.spec[["unit"]][["Unit"]],
        Bin_start_incl=forecast.type[["Bin_start_incl_for"]](target.spec, ...),
        Bin_end_notincl=forecast.type[["Bin_end_notincl_for"]](target.spec, ...),
        check.names=FALSE,
        stringsAsFactors=FALSE
      ) %>>%
      tibble::as_tibble() %>>%
      stats::setNames(tolower(names(.)))
    if (anyDuplicated(tolower(names(subspreadsheet))) != 0L) {
      stop ("Duplicate names in subspreadsheet when case-insensitive; need them to be unique after lower-casing")
    }
    subspreadsheet <- subspreadsheet %>>%
      stats::setNames(tolower(names(.)))
    if (!identical(index.tbl, subspreadsheet[names(index.tbl)])) {
      print("AAA")
      print(index.tbl)
      print("BBB")
      print(subspreadsheet[names(index.tbl)])
      print("CCC")
      missing.indices =
        dplyr::anti_join(index.tbl, subspreadsheet, names(index.tbl))
      if (nrow(missing.indices) != 0L) {
        stop (paste0("Missing indices:\n", paste(capture.output(print(missing.indices)),collapse="\n")))
      }
      extra.entries =
        dplyr::anti_join(subspreadsheet, index.tbl, names(index.tbl))
      if (nrow(extra.entries) != 0L) {
        stop (paste0("Extra/mislabeled entries:\n", paste(capture.output(print(extra.entries)),collapse="\n")))
      }
      if (nrow(index.tbl) != nrow(subspreadsheet)) {
        stop ("Number of rows in spreadsheet did not match the expected value; should be due to duplicate indices in subspreadsheet.")
      }
      stop (paste0("(After lowercasing colnames) index columns were not identical to expectations; feedback from all.equal (which may not cover all differences): ", capture.output(str(all.equal(subspreadsheet[names(index.tbl)], index.tbl, check.attributes=TRUE, ignore.col.order=FALSE, ignore.row.order=FALSE, tolerance=0)))))
    }
    forecast.type[["forecast_value_from_Value"]](subspreadsheet[["value"]], target.spec, ...)
  }

flusight2016.proxy.forecast.types = list(
  point.mae.forecast.type,
  distr.logscore.forecast.type
) %>>%
  setNames(sapply(., magrittr::extract2, "Type"))

flusight2016.evaluation.forecast.types = list(
  point.mae.forecast.type, # not really part of FluSight 2016 evaluations;
                           # including for extra evaluation information
  multibin.logscore.forecast.type
) %>>%
  setNames(sapply(., magrittr::extract2, "Type"))

forecast_value = function(target.spec, forecast.type, target.forecast, label.bins=TRUE) {
  if (length(target.forecast[["target.values"]]) != 1L) {
    stop ("Expected target.forecast$target.values to have length 1L.")
  } else {
    target.i = 1L
  }
  target.values = target.forecast[["target.values"]][[target.i]]
  target.weights = target.forecast[["target.weights"]]
  do.call(
    forecast.type[["forecast_value_from_weighted_univals"]],
    c(list(target.spec, target.values, target.weights,
           label.bins=label.bins),
      target.forecast[["target.settings"]],
      target.forecast[["method.settings"]])
  )
}

observed_value = function(target.spec, target.settings, forecast.type, observed.multival, label.bins=TRUE) {
  do.call(
    forecast.type[["observed_value_from_observed_multival"]],
    c(list(target.spec, observed.multival, label.bins=label.bins),
      target.settings))
}

## xxx Are break.bin.representatives used for labeling? If not, they can be shifted and the shift_for_smoothing functions can be removed.
