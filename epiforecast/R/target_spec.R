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

##' @import pipeR
NULL

##' Calculate the peak week(s) in a vector of weekly observations
##'
##' @param trajectory a vector of weekly observations
##' @param is.inseason logical vector length-compatible with \code{trajectory},
##' @param ... ignored
##' @return integer vector, typically of length 1; indices of elements that are
##'   \code{>=} all other elements
##'
##' @export
pwk = function(trajectory, is.inseason=rep_len(TRUE, length(trajectory)), ...) {
  if (any(is.inseason & is.na(trajectory))) {
    print(is.inseason)
    print(trajectory)
    print(length(is.inseason))
    print(length(trajectory))
    stop ("In-season NA elements are not allowed when calculating =pwk=.")
  }
  return (which(is.inseason & trajectory==max(trajectory[is.inseason])))
}

##' Calculate the peak height of a trajectory
##'
##' @param trajectory a vector
##' @param is.inseason logical vector length-compatible with \code{trajectory},
##' @param ... ignored
##' @return a scalar --- the maximum value of the vector
##'
##' @export
pht = function(trajectory, is.inseason=rep_len(TRUE, length(trajectory)), ...){
  return (max(trajectory[is.inseason]))
}

##' Calculate the onset of a trajectory (post-rounding)
##'
##' @param trajectory a vector
##' @param baseline the onset threshold
##' @param is.inseason logical vector length-compatible with \code{trajectory},
##'   acting as a mask on possible output values
##' @param ... ignored
##' @return the first index which is part of the in-season and is part of a run
##'   of consecutive observations above the onset threshold that lasts at least
##'   for two additional indices
##'
##' @export
ons = function(trajectory, baseline, is.inseason, ...) {
  above.baseline = trajectory >= baseline
  next.three.above.baseline =
      above.baseline &
      dplyr::lead(above.baseline, 1L) &
      dplyr::lead(above.baseline, 2L)
  return (which(is.inseason & next.three.above.baseline)[1L][[1L]])
}

##' Calculate the duration of a trajectory
##'
##' @param trajectory a vector
##' @param baseline the onset threshold
##' @param is.inseason logical vector length-compatible with \code{trajectory},
##'   acting as a mask on possible output values
##' @param ... ignored
##' @return a single non-\code{NA} integer: the number of indices which are part
##'   of the in-season and part of a run of at least three consecutive
##'   observations above the onset threshold
##'
##' @export
dur = function(trajectory, baseline, is.inseason, ...) {
  above.baseline = trajectory >= baseline
  next.three.above.baseline =
    above.baseline &
    dplyr::lead(above.baseline, 1L) &
    dplyr::lead(above.baseline, 2L)
  part.of.run =
    (next.three.above.baseline |
     dplyr::lag(next.three.above.baseline, 1L) |
     dplyr::lag(next.three.above.baseline, 2L))
  return (sum(is.inseason & part.of.run, na.rm=TRUE))
}

##' First epi week in weekly USA flu-related trajectories used by some methods
##'
##' @export
usa.flu.first.week.of.season = 31L

##' Get \code{logical} marking in-season times in a USA flu-related trajectory
##'
##' Given that a trajectory starts with the epi week
##' \code{\link{usa.flu.first.week.of.season}} and is \code{n.weeks.in.season}
##' long, marks parts of the trajectory that are part of the USA flu "in-season"
##' (epi week 40 of the first year of a season to epi week 20 of the following
##' year)
##'
##' @param n.weeks.in.season \code{52L} or \code{53L}; the number of weeks in
##'   the first year of a season
##'
##' @return logical vector of length equal to \code{n.weeks.in.season}; entries
##'   are \code{TRUE} at times corresponding to weeks that are part of the
##'   in-season, and \code{FALSE} otherwise (i.e., in the off-season)
##'
##' @export
usa_flu_inseason_flags = function(n.weeks.in.season) {
  seq_len(n.weeks.in.season) %>>%
    time_to_model_week(usa.flu.first.week.of.season) %>>%
    model_week_to_epi_week(usa.flu.first.week.of.season, n.weeks.in.season) %>>%
    dplyr::between(21L, 39L) %>>%
    magrittr::not() %>>%
    return()
}

flusurv_inseason_flags = function(n.weeks.in.season) {
  seq_len(n.weeks.in.season) %>>%
    time_to_model_week(usa.flu.first.week.of.season) %>>%
    model_week_to_epi_week(usa.flu.first.week.of.season, n.weeks.in.season) %>>%
    dplyr::between(18L, 39L) %>>%
    magrittr::not() %>>%
    return()
}

time_to_point = function(times, is.inseason, ...) {
  n.weeks.in.season = length(is.inseason)
  model.weeks <- time_to_model_week(times, usa.flu.first.week.of.season)
  epi.weeks <- model_week_to_epi_week(model.weeks, usa.flu.first.week.of.season, n.weeks.in.season)
  return (epi.weeks)
}
## fixme don't hardcode usa.flu.first.week.of.season

time_to_binlabelstart = function(times, is.inseason, ...) {
  epi.weeks = time_to_point(times, is.inseason, ...)
  epi.weeks.char = dplyr::coalesce(as.character(epi.weeks), "none")
}

time_to_binlabelend = function(times, is.inseason, ...) {
  epi.weeks.ish = 1L+time_to_point(times-1L, is.inseason, ...)
  epi.weeks.char = dplyr::coalesce(as.character(epi.weeks.ish), "none")
}

time_from_binlabel = function(epi.weeks.char, is.inseason, ...) {
  n.weeks.in.season = length(is.inseason)
  epi.weeks.numeric = as.numeric(dplyr::recode(epi.weeks.char, "none"=NA_character_))
  epi.weeks.integer = as.integer(epi.weeks.numeric)
  epi.weeks = if (all(epi.weeks.numeric==epi.weeks.integer)) {
                epi.weeks.integer
              } else {
                epi.weeks.numeric
              }
  model.weeks = epi_week_to_model_week(epi.weeks, usa.flu.first.week.of.season, n.weeks.in.season)
  times = model_week_to_time(model.weeks, usa.flu.first.week.of.season)
  return (times)
}

time_from_point = function(epi.weeks, is.inseason, ...) {
  n.weeks.in.season = length(is.inseason)
  model.weeks = epi_week_to_model_week(epi.weeks, usa.flu.first.week.of.season, n.weeks.in.season)
  times = model_week_to_time(model.weeks, usa.flu.first.week.of.season)
  return (times)
}

week.unit = list(
  Unit = "week",
  to_binlabelstart = function(target.values, is.inseason, ...) {
    time_to_binlabelstart(target.values, is.inseason, na.representation="none", ...)
  },
  to_binlabelend = function(target.values, is.inseason, ...) {
    time_to_binlabelend(target.values, is.inseason, na.representation="none", ...)
  },
  from_binlabelstart = function(binstrings, is.inseason, ...) {
    time_from_binlabel(binstrings, is.inseason, na.representation="none", ...)
  },
  to_point = function(target.values, is.inseason, ...) {
    time_to_point(target.values, is.inseason, na.representation="none", ...)
  },
  from_point = function(points, is.inseason, ...) {
    time_from_point(points, is.inseason, na.representation="none", ...)
  },
  shift_for_smoothing = function(target.values, ...) {
    target.values + 0.5
  }
)

percentage.unit = list(
  Unit = "percent",
  to_binlabelstart = function(x, ...) as.character(x, ...),
  to_binlabelend = function(x, ...) as.character(x, ...),
  from_binlabel = function(x, ...) as.numeric(x, ...),
  to_point = function(x, ...) x,
  from_point = function(x, ...) as.numeric(x),
  shift_for_smoothing = function(x, ...) {
    x + 0.05
  }
)

flusight2016.onset.target.spec = list(
  Target = "Season onset",
  unit = week.unit,
  for_processed_trajectory = function(processed.trajectory, baseline, is.inseason, ...) {
    if (length(is.inseason) != length(processed.trajectory)) {
      print("traj")
      print(processed.trajectory)
      print(length(processed.trajectory))
      print("seas")
      print(is.inseason)
      print(length(is.inseason))
      stop ("length(is.inseason) != length(processed.trajectory)")
    }
    if (!all(0 <= processed.trajectory & processed.trajectory <= 100 &
             round(processed.trajectory, 1L) == processed.trajectory)) {
      stop ("Observations in processed.trajectory must lie between 0 and 100, inclusive, and be rounded to the first decimal place.")
    }
    return (ons(processed.trajectory, baseline, is.inseason))
  },
  bin_info_for = function(is.inseason, ...) {
    ## 1 bin per in-season week (need to add a break at the end of the in-season
    ## week indices for the last week index), plus an NA(=none) bin:
    breaks = which(is.inseason) %>>% c(.[[length(.)]]+1L)
    return (list(
      breaks=breaks,
      break.bin.representatives=head(breaks, -1L),
      rightmost.closed=FALSE, include.na=TRUE
    ))
  }
)

flusight2016.peak.week.target.spec = list(
  Target = "Season peak week",
  unit = week.unit,
  for_processed_trajectory = function(processed.trajectory, is.inseason, ...) {
    return (pwk(processed.trajectory, is.inseason, ...))
  },
  bin_info_for = function(is.inseason, ...) {
    ## 1 bin per in-season week (need to add a break at the end of the in-season
    ## week indices for the last week index), with no NA(=none) bin (unlike
    ## Season onset):
    breaks = which(is.inseason) %>>% c(.[[length(.)]]+1L)
    return (list(
      breaks=breaks,
      break.bin.representatives=head(breaks, -1L),
      rightmost.closed=FALSE, include.na=FALSE
    ))
  }
)

flusight2016.percentage.bin.info = list(
  breaks=c(0:130/10, 100),
  break.bin.representatives = 0:130/10,
  rightmost.closed=TRUE, include.na=FALSE
)

flusight2016.peak.percentage.target.spec = list(
  Target = "Season peak percentage",
  unit = percentage.unit,
  for_processed_trajectory = function(processed.trajectory, is.inseason, ...) {
    return (pht(processed.trajectory, is.inseason, ...))
  },
  bin_info_for = function(...) {
    return (flusight2016.percentage.bin.info)
  }
)

flusight2016_percentage_target_spec_for_lookahead = function(lookahead) {
  list(
    Target = paste0(lookahead, " wk ahead"),
    unit = percentage.unit,
    for_processed_trajectory = function(processed.trajectory, forecast.time, ...) {
      return (processed.trajectory[[forecast.time+lookahead]])
    },
    bin_info_for = function(...) {
      return (flusight2016.percentage.bin.info)
    }
  )
}

flusight2016.target.specs = list(
  flusight2016.onset.target.spec,
  flusight2016.peak.week.target.spec,
  flusight2016.peak.percentage.target.spec,
  flusight2016_percentage_target_spec_for_lookahead(1L),
  flusight2016_percentage_target_spec_for_lookahead(2L),
  flusight2016_percentage_target_spec_for_lookahead(3L),
  flusight2016_percentage_target_spec_for_lookahead(4L)
) %>>%
  setNames(sapply(., magrittr::extract2, "Target"))

flusight2016_target_trajectory_preprocessor = function(trajectory) {
  round(pmin(pmax(trajectory, 0), 100), 1L)
}

flusurv2017.age.group.percentage.bin.infos = c(
  "Overall"=130L,
  "0-4 yr"=130L,
  "5-17 yr"=130L,
  "18-49 yr"=130L,
  "50-64 yr"=130L,
  "65+ yr"=600L
) %>>%
  lapply(function(max.bin.start.times.ten) {
    list(
      breaks=c(0:max.bin.start.times.ten/10, 100),
      break.bin.representatives = 0:max.bin.start.times.ten/10,
      rightmost.closed=TRUE, include.na=FALSE
    )
  })

flusurv2017.peak.week.target.spec = flusight2016.peak.week.target.spec

flusurv2017.peak.percentage.target.spec = list(
  Target = "Season peak percentage",
  unit = percentage.unit,
  for_processed_trajectory = function(processed.trajectory, is.inseason, ...) {
    return (pht(processed.trajectory, is.inseason, ...))
  },
  bin_info_for = function(age.group, ...) {
    return (flusurv2017.age.group.percentage.bin.infos[[age.group]])
  }
)

flusurv2017_percentage_target_spec_for_lookahead = function(lookahead) {
  list(
    Target = paste0(lookahead, " wk ahead"),
    unit = percentage.unit,
    for_processed_trajectory = function(processed.trajectory, forecast.time, ...) {
      return (processed.trajectory[[forecast.time+lookahead]])
    },
    bin_info_for = function(age.group, ...) {
      return (flusurv2017.age.group.percentage.bin.infos[[age.group]])
    }
  )
}

flusurv2017.target.specs = list(
  flusurv2017.peak.week.target.spec,
  flusurv2017.peak.percentage.target.spec,
  flusurv2017_percentage_target_spec_for_lookahead(1L),
  flusurv2017_percentage_target_spec_for_lookahead(2L),
  flusurv2017_percentage_target_spec_for_lookahead(3L),
  flusurv2017_percentage_target_spec_for_lookahead(4L)
) %>>%
  setNames(sapply(., magrittr::extract2, "Target"))

flusurv2017_target_trajectory_preprocessor = function(trajectory) {
  result = c(rep(NA_real_, 9L), round(pmin(pmax(trajectory, 0), 100), 1L), rep(NA_real_, 13))
  if (length(result) != 52L && length(result) != 53L) {
    stop ("Trajectory is the wrong length, or there was an error in appending NA's to the ends.")
  }
  return (result)
}

Unit_for_target.spec = function(target.spec) {
  target.spec[["unit"]][["Unit"]]
}

bin_representatives_for_target.spec = function(target.spec, voxel.data) {
  bin.info = do.call(target.spec[["bin_info_for"]], voxel.data[["target.settings"]])
  bin.representatives = c(
    bin.info[["break.bin.representatives"]],
    do.call(get_na_value_or_empty_for_target,
            c(list(target.spec),
              voxel.data[["target.settings"]]))
  )
  return (bin.representatives)
}

## todo names in target list vs. Target in target...
## todo unit as sublist or concatenated?
