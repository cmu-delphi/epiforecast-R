
##' @include utils.R
##' @include match.R
NULL

##' Returns a possibly-named integer-class vector of weekday numbers in 0:6 (\%w
##' format) given an integer-valued numeric vector with values in 0:7 (\%w or \%u
##' format), or throws an error if the input seems inappropriate.
##'
##' @param wday vector of weekday numbers (each \code{\%in\% 0:7})
##'
##' @return integer-class \%w format weekday numbers
match.wday.w = function(wday) {
  wday <- match.integer(wday)
  if (!all(wday %in% 0:7)) {
    stop (sprintf("All values of =%s= must be =%%in%% 0:7=.",
                  paste(deparse(match.call()$wday), collapse="")))
  }
  wday[wday==7] <- 0
  return (wday)
}

##' Returns a possibly-named length-1 integer-class vector of weekday numbers in 0:6 (\%w
##' format) given an integer-valued length-1 numeric vector with values in 0:7 (\%w or \%u
##' format), or throws an error if the input seems inappropriate.
##'
##' @param wday vector of weekday numbers (each \code{\%in\% 0:7})
##'
##' @return integer-class \%w format weekday numbers
match.single.wday.w = function(wday) {
  wday <- match.single.nonna.integer(wday)
  if (!all(wday %in% 0:7)) {
    stop (sprintf("All values of =%s= must be =%%in%% 0:7=.",
                  paste(deparse(match.call()$wday), collapse="")))
  }
  wday[wday==7] <- 0
  return (wday)
}

##' Convert dates to a data frame of year-week-wday according to a specified
##' convention.
##'
##' @param date object compatible with \code{\link{as.Date}}: the dates to
##'   convert
##' @param first.wday weekday number(s) (0--7, Sunday can be 0 or 7): the
##'   weekday that is considered the beginning of the week; typically Sunday or
##'   Monday
##' @param owning.wday weekday number(s) (0--7, Sunday can be 0 or 7): a week is
##'   assigned to a given year if the owning weekday of that week falls in that
##'   year; typically \code{first.wday} or \code{first.wday+3}
##'
##' @return a data frame with three columns, \code{$year}, \code{$week}, and
##'   \code{$wday}, corresponding to the given dates using the convention
##'   specified by \code{first.wday} and \code{owning.wday}; each wday entry
##'   will be \code{\%in\% 0:6}.
##'
##' @export
DateToYearWeekWdayDF = function(date, first.wday, owning.wday) {
  ## Do some checks and make wday's lie in 0:6 (%w format).
  date <- as.Date(date)
  first.wday <- match.single.wday.w(first.wday)
  owning.wday <- match.single.wday.w(owning.wday)
  n = names(date)

  if (any(is.na(date))) {
    stop("NA dates are not currently supported (they at least interfere with some bug-detection checks).")
  }

  ## Extract %w wday's from input |date|:
  input.wday = lubridate::wday(date) - 1

  ## Determine the offsets between the weekdays of the input, first
  ## weekday, and owning weekday.  R's =%%= on testing system appears to
  ## take sign of RHS, which is what we want; however, just in case,
  ## make sure:
  first.to.input = (input.wday - first.wday + 7) %% 7
  first.to.owning = (owning.wday - first.wday + 7) %% 7
  if (!all(first.to.input %in% 0:6) || !all(first.to.owning %in% 0:6)) {
    stop ("Bug detected: logic error when determining weekday shifts (1st check).")
  }

  owning.date = date - first.to.input + first.to.owning
  if (!all(as.integer(format(owning.date,"%w"))==owning.wday))
    stop("Bug detected: logic error when determining weekday shifts (2nd check).")

  owning.year = as.integer(lubridate::year(owning.date))
  owning.week = (as.integer(lubridate::yday(owning.date))-1L) %/% 7L + 1L

  result = data.frame(year=owning.year, week=owning.week, wday=input.wday)
  rownames(result) <- n
  return (result)
}

##' Convert a structure containing year, week, and weekday numbers into the
##' corresponding dates under a specified week numbering convention. Reverses
##' \code{\link{DateToYearWeekWdayDF}}.
##'
##' @param ywwd data frame with columns \code{$year}, \code{$week}, and
##'   \code{$wday} containing the year, week, and weekday numbers respectively
##'   (weekday numbers are 0--7; Sunday can be inputted either as 0 or as 7)
##' @param first.wday wday number corresponding to the first weekday of any week
##' @param owning.wday wday number; a week is assigned to a given year if the
##'   owning weekday of that week falls in the given year
##' @param error.on.wrap TRUE or FALSE: if TRUE, an error is generated if a
##'   given week falls outside the range of possible week numbers for the
##'   corresponding year (i.e., if it is nonpositive or greater than the number
##'   of weeks assigned to the corresponding year in \code{year}); if FALSE, the
##'   week number will wrap around to future or previous years (e.g., week 0 of
##'   1997 will be considered the last week of 1996)
##'
##' @return a \code{Date} vector: the dates corresponding to the year-week-wdays
##'   and week numbering convention
##'
##' @export
yearWeekWdayDFToDate = function(ywwd, first.wday, owning.wday, error.on.wrap=TRUE) {
  n = rownames(ywwd)
  if (identical(n, as.character(seq_len(nrow(ywwd))))) {
    n <- NULL
  }
  return (yearWeekWdayVecsToDate(structure(ywwd$year, names=n), ywwd$week, ywwd$wday,
                                 first.wday, owning.wday,
                                 error.on.wrap=error.on.wrap))
}

##' Like \code{\link{yearWeekWdayDFToDate}}, but with each column of \code{ywwd}
##' provided as separate parameters.
##'
##' @param year integer-valued vector: the year numbers
##' @param week integer-valued vector: the week numbers
##' @param wday integer-valued vector: the wday numbers (0--7, either 0 or 7 can
##'   be used for Sunday)
##' @param first.wday wday number of the first day of each week
##' @param owning.wday if the owning weekday of a week falls in a particular
##'   year, the entire week is assigned to that year
##' @param error.on.wrap if TRUE, throws errors when week numbers do not fall in
##'   the specified years; when FALSE, wraps them around into other years
##'
##' @return the corresponding dates
##'
##' @export
yearWeekWdayVecsToDate = function(year, week, wday, first.wday, owning.wday, error.on.wrap=TRUE) {
  year <- match.integer(year)
  week <- match.integer(week)
  wday <- match.wday.w(wday)
  first.wday <- match.single.wday.w(first.wday)
  owning.wday <- match.single.wday.w(owning.wday)
  len = max(sapply(list(year,week,wday,first.wday,owning.wday), length))
  if (!is.null(names(year))) {
    n = rep_len(names(year), len)
  } else if (!is.null(names(week))) {
    n = rep_len(names(week), len)
  } else if (!is.null(names(first.wday))) {
    n = rep_len(names(first.wday), len)
  } else if (!is.null(names(owning.wday))) {
    n = rep_len(names(owning.wday), len)
  } else {
    n = NULL
  }

  if (error.on.wrap) {
    if (any(week <= 0)) {
      stop ("Week numbers must be positive (in particular, week 0 is not supported) when |error.on.wrap| is TRUE.")
    } else if (any(week > lastWeekNumber(year, owning.wday))) {
      stop ("Week numbers must not exceed the last week number in the corresponding year unless |error.on.wrap| is FALSE.")
    }
  }

  ## Determine the offsets between the weekdays of the input, first
  ## weekday, and owning weekday.  R's |%%| on testing system appears to
  ## take sign of RHS, which is what we want; however, just in case,
  ## make sure:
  first.to.input = (wday - first.wday + 7) %% 7
  first.to.owning = (owning.wday - first.wday + 7) %% 7
  if (!all(first.to.input %in% 0:6) || !all(first.to.owning %in% 0:6)) {
    stop ("Bug detected: logic error when determining weekday shifts.")
  }

  jan1 = as.Date(sprintf("%d-01-01",year))
  ## lubridate assigns wday 1 to Sunday; adjust:
  jan1.wday = lubridate::wday(jan1) - 1
  jan1.to.owning = (owning.wday - jan1.wday + 7) %% 7
  if (!all(jan1.to.owning %in% 0:6)) {
    stop ("Bug detected: logic error when determining weekday shifts.")
  }
  week1.owning = jan1 + jan1.to.owning

  date = week1.owning + (week-1)*7 - first.to.owning + first.to.input
  if (!all(as.integer(format(date,"%w"))==wday)) {
    stop ("Bug detected: logic error when determining weekday shifts.")
  }

  names(date) <- n
  return (date)
}

##' Like \code{\link{yearWeekWdayDFToDate}}, but allows parameters to be
##' provided in several (is.)list structures and vectors.
##'
##' @param ... vectors and (is.)lists to feed to
##'   \code{\link{yearWeekWdayVecsToDate}}; each entry in a list (e.g., column
##'   in a data.frame) is treated as a vector. Names can be provided by the
##'   names attribute of list arguments or optional parameter names for the
##'   vector arguments. (With no names, parameters of
##'   \code{\link{yearWeekWdayVecsToDate}} are matched sequentially.):
##'
##' year integer-valued vector: the year numbers
##'
##' week integer-valued vector: the week numbers
##'
##' wday integer-valued vector: the wday numbers (0--7, either 0 or 7 can be
##'   used for Sunday)
##'
##' first.wday wday number of the first day of each week
##'
##' owning.wday if the owning weekday of a week falls in a particular year, the
##'   entire week is assigned to that year
##'
##' error.on.wrap if TRUE, throws errors when week numbers do not fall in the
##'   specified years; when FALSE, wraps them around into other years
##'
##' @return the corresponding dates
##'
##' @export
yearWeekWdayListsToDate = function(...) {
  args = list(...)
  args <- lapply(args, function(arg) {
    listified.arg = as.list(arg)
    if (is.data.frame(arg)) {
      for (col.i in seq_along(listified.arg)) {
        names(listified.arg[[col.i]]) <- rownames(arg)
      }
    }
    listified.arg
  })
  args <- do.call(c,args)
  return (do.call(yearWeekWdayVecsToDate, args))
}

##' The number of weeks assigned to a given year or years.
##'
##' The length of one input should be evenly divisible by the length of the
##' other.
##'
##' @param year integer/NA-valued vector of years
##' @param owning.wday integer vector of wday numbers (0--7, 0 and 7 are both
##'   Sunday): a week is assigned to a given year based on whether this weekday
##'   is contained in that year
##'
##' @return the number of weeks assigned to each year in \code{years}
##'
##' @examples
##' ## The number of epi weeks in each year from 1990 to 2020:
##' lastWeekNumber(1990:2020, 3)
##'
##' @export
lastWeekNumber = function(year, owning.wday) {
  year <- match.integer(year)
  owning.wday <- match.wday.w(owning.wday)
  return (structure(52L+(year==lubridate::year(yearWeekWdayVecsToDate(year, 53L, owning.wday, 0L, owning.wday, error.on.wrap=FALSE)))))
}

namedSeason = function(season) {
  season <- match.integer(season)
  names(season) <- paste0("S", season)
  return (season)
}

##' Get the season number associated with a particular year and week
##'
##' Seasons are 52/53-week-long time spans that start with a particular week
##' number (from 1 to 52); years are a special case that start with week 1.
##' Seasons contain weeks from two consecutive years except in the case that the
##' season starts with week 1, in which case it coincides exactly with a single
##' year. Each season ends the week before the next season (with the same
##' starting week) begins. Seasons are numbered by the first distinct year from
##' which they take weeks, and labeled with an "S" prefix; for example, a season
##' containing weeks from 2015 and 2016 would be numbered 2015 and labeled
##' "S2015". This function gives the season number associated with the inputted
##' weeks specified as year-week combinations.
##'
##' @param year integer-valued vector: year in which the weeks fall
##' @param week integer-valued vector: associated week numbers (each \code{\%in\% 1:53})
##' @param first.week the week on which each season should start
##'
##' @return named integer-valued vector: season numbers associated with the
##'   inputted weeks with names giving the season labels
##'
##' @export
seasonOfYearWeek = function(year, week, first.week) {
  year <- match.integer(year)
  week <- match.integer(week)
  first.week <- match.integer(first.week)
  if (!all(1L <= first.week & first.week <= 52L)) {
    violators = unique(first.week)
    violators <- violators[!(1L <= violators & violators <= 52L)]
    stop(paste("Season starting weeks must be %in% 1:52; violators: ", violators, collapse=","))
  } else {
    first.year = year - as.integer(week < first.week)
    return (namedSeason(first.year))
  }
}

##' Get the season number associated with a particular date
##'
##' Seasons are 52/53-week-long time spans that start with a particular week
##' number (from 1 to 52); years are a special case that start with week 1.
##' Seasons contain weeks from two consecutive years except in the case that the
##' season starts with week 1, in which case it coincides exactly with a single
##' year. Seasons are numbered by the first distinct year from which they take
##' weeks, and labeled with an "S" prefix; for example, a season containing
##' weeks from 2015 and 2016 would be numbered 2015 and labeled "S2015". This
##' function gives the season number associated with the dates.
##'
##' @param date object convertible to \code{Date}
##' @param first.week the week on which each season should start
##' @param first.wday weekday number(s) (0--7, Sunday can be 0 or 7): the
##'   weekday that is considered the beginning of the week; typically Sunday or
##'   Monday
##' @param owning.wday weekday number(s) (0--7, Sunday can be 0 or 7): a week is
##'   assigned to a given year if the owning weekday of that week falls in that
##'   year; typically \code{first.wday} or \code{first.wday+3}
##'
##' @return named integer-valued vector: season numbers associated with the
##'   inputted weeks with names giving the season labels
##'
##' @export
seasonOfDate = function(date, first.week, first.wday, owning.wday) {
  date <- as.Date(date)
  first.week <- match.integer(first.week)
  ywwd = DateToYearWeekWdayDF(date, first.wday, owning.wday)
  return (seasonOfYearWeek(ywwd$year, ywwd$week, first.week))
}

##' Get the first weekday of every week in the given seasons
##'
##' @param season integer-valued vector: season numbers
##' @param first.week integer-valued vector: first week number of each season
##' @param first.wday integer-valued vector of weekday numbers: first weekday
##'   number in each week
##' @param owning.wday integer-valued vector of weekday numbers: weekday that
##'   determines the year to which a week will be assigned
##'
##' @return $\code{length(season)}$-length named list of 52/53-length
##'   \code{Date} vectors
##'
##' @export
DatesOfSeason = function(season, first.week, first.wday, owning.wday) {
  season <- match.integer(season)
  first.week <- match.integer(first.week)
  if (is.null(names(season))) {
    season <- namedSeason(season)
  }
  if (length(first.week) < length(season)) {
    if (length(season) %% length(first.week) != 0) {
      stop("The greater of =length(season)= and =length(first.week)= must be a multiple of the lesser.")
    }
    first.week <- rep_len(first.week, length(season))
  }
  first.year = season
  last.year = season+1L
  first.date = yearWeekWdayVecsToDate(first.year, first.week, first.wday, first.wday, owning.wday)
  last.date = yearWeekWdayVecsToDate(last.year, first.week, first.wday, first.wday, owning.wday)-1
  season.dates = mapply(seq.Date, first.date, last.date, by="1 week", SIMPLIFY=FALSE)
  names(season.dates) <- names(season)
  ## Validate output: all dates generated for a given season should indeed by part of that season:
  if (!isTRUE(all.equal(as.list(season),
                        structure(lapply(seq_along(season.dates), function(season.i) {
                          unique(seasonOfDate(season.dates[[season.i]], first.week[season.i], first.wday, owning.wday))
                        }), names=names(season))))) {
    stop("Bug detected: =seasonDates= season.datess inconsistent with =seasonOf=.")
  }
  return (season.dates)
}

##' Convert year-week to season-model.week
##'
##' Seasons starting with a week number \code{n} other than 1 contain week
##' numbers \code{n} to 52/53 from the season's first year and 1 to \code{n-1}
##' of the season's second year. Sometimes we would like a numbering of weeks
##' within a season that coincides with the week number when possible, but does
##' not have a jump down from 52/53 to 1. Model weeks fulfill this purpose: they
##' begin with the starting week of a season and increase by 1 for each
##' subsequent week; they coincide with week numbers in the starting year, and
##' are 52/53 plus the week number in the second year.
##'
##' @param year integer-valued vector: year from the year-week numbering
##' @param week integer-valued vector: week from the year-week numbering
##' @param first.week integer-valued vector: week number that the seasons start
##'   on
##' @param owning.wday integer-valued vector of weekday numbers: weekday that
##'   determines the year to which a week will be assigned
##'
##' @return \code{data.frame} with two columns giving the corresponding
##'   season-model.week breakdown:
##'
##' \code{$season}: integer-class vector: season numbers
##'
##' \code{$model.week}: integer-class vector: model week numbers
yearWeekToSeasonModelWeekDF = function(year, week, first.week, owning.wday) {
    year <- as.integer(year)
    week <- as.integer(week)
    first.week <- as.integer(first.week)
    owning.wday <- as.integer(owning.wday)
    season = seasonOfYearWeek(year, week, first.week)
    model.week = ifelse(week >= first.week, week, lastWeekNumber(year-1L, owning.wday)+week)
    return (data.frame(season=season, model.week=model.week))
}

##' Convert year-week in a \code{data.frame} to season-model.week
##'
##' Delegates to \code{\link{yearWeekToSeasonModelWeekDF}}.
##'
##' @param yearWeek \code{data.frame} (or other list) with columns \code{$year}
##'   and \code{$week}:
##'
##' \code{$year} integer-valued vector: year from the year-week numbering
##'
##' \code{$week} integer-valued vector: week from the year-week numbering
##' @param first.week integer-valued vector: week number that the seasons start
##'   on
##' @param owning.wday integer-valued vector of weekday numbers: weekday that
##'   determines the year to which a week will be assigned
##'
##' @return \code{data.frame} with two columns giving the corresponding
##'   season-model.week breakdown:
##'
##' \code{$season}: integer-class vector: season numbers
##'
##' \code{$model.week}: integer-class vector: model week numbers
##'
##' @export
yearWeekDFToSeasonModelWeekDF = function(yearWeek, first.week, owning.wday) {
  return (yearWeekToSeasonModelWeekDF(yearWeek$year, yearWeek$week, first.week, owning.wday))
}

##' Convert season-model.week to year-week
##'
##' Like \code{\link{yearWeekToSeasonModelWeekDF}}, but in opposite direction.
##'
##' @param season integer-valued vector: season numbers
##' @param model.week integer-valued vector: model week numbers
##' @param first.week integer-valued vector: week number that the seasons start
##'   on
##' @param owning.wday integer-valued vector of weekday numbers: weekday that
##'   determines the year to which a week will be assigned
##'
##' @return \code{data.frame} with two columns giving the corresponding
##'   year-week breakdown:
##'
##' \code{$year}: integer-class vector: years
##'
##' \code{$week}: integer-class vector: week numbers
##'
##' @export
seasonModelWeekToYearWeekDF = function(season, model.week, first.week, owning.wday) {
  season <- match.integer(season)
  model.week <- match.integer(model.week)
  first.week <- match.integer(first.week)
  owning.wday <- match.wday.w(owning.wday)
  if (length(season) < length(model.week)) {
    season <- rep_len(season, length(model.week))
  } else if (length(model.week) < length(season)) {
    model.week <- rep_len(model.week, length(season))
  }
  first.week <- as.integer(first.week)
  owning.wday <- as.integer(owning.wday)
  n.weeks = lastWeekNumber(season, owning.wday)
  year = ifelse(model.week > n.weeks, season+1L, season)
  week = ifelse(model.week > n.weeks, model.week - n.weeks, model.week)
  return (data.frame(year=year, week=week))
}

##' Convert season-model.week in \code{data.frame} to year-week
##'
##' Like \code{\link{yearWeekDFToSeasonModelWeekDF}}, but in opposite direction.
##'
##' @param seasonModelWeek \code{data.frame} (or other list) with columns
##'
##' \code{$season}: integer-valued vector: season numbers
##'
##' \code{$model.week}: integer-valued vector: model week numbers
##' @param first.week integer-valued vector: week number that the seasons start
##'   on
##' @param owning.wday integer-valued vector of weekday numbers: weekday that
##'   determines the year to which a week will be assigned
##'
##' @return \code{data.frame} with two columns giving the corresponding
##'   year-week breakdown:
##'
##' \code{$year}: integer-class vector: years
##'
##' \code{$week}: integer-class vector: week numbers
##'
##' @examples
##' dates = as.Date("2015-01-01")+seq.int(0L, 1000L, 7L)
##' ywwd = DateToYearWeekWdayDF(dates, 0L, 3L) # epi week convention
##' dates.duplicate1 = yearWeekWdayDFToDate(ywwd, 0L, 3L)
##' identical(dates, dates.duplicate1)
##' smw = yearWeekDFToSeasonModelWeekDF(ywwd, 21L, 3L) # seasons starting on week number 21
##' yw = seasonModelWeekDFToYearWeekDF(smw, 21L, 3L)
##' identical(ywwd[,c("year","week")], yw)
##' dates.duplicate2 = seasonModelWeekWdayDFToDate(cbind(smw, wday=ywwd$wday), 21L, 0L, 3L)
##' identical(dates, dates.duplicate2)
##'
##' @export
seasonModelWeekDFToYearWeekDF = function(seasonModelWeek, first.week, owning.wday) {
  return (seasonModelWeekToYearWeekDF(seasonModelWeek$season, seasonModelWeek$model.week, first.week, owning.wday))
}

##' Convert season-model.week-wday to \code{Date}
##'
##' @param season season number
##' @param model.week model week number
##' @param wday weekday number
##' @param first.week first week number of season
##' @param first.wday first weekday number of week
##' @param owning.wday integer-valued vector of weekday numbers: weekday that
##'   determines the year to which a week will be assigned
##' @param error.on.wrap TRUE or FALSE: if TRUE, an error is generated if a
##'   given week falls outside the range of possible week numbers for the
##'   corresponding year (i.e., if it is nonpositive or greater than the number
##'   of weeks assigned to the corresponding year in \code{year}); if FALSE, the
##'   week number will wrap around to future or previous years (e.g., week 0 of
##'   1997 will be considered the last week of 1996)
##'
##' @return \code{Date} vector with the corresponding dates
##'
##' @export
seasonModelWeekWdayToDate = function(season, model.week, wday, first.week, first.wday, owning.wday, error.on.wrap=TRUE) {
  season <- match.integer(season)
  model.week <- match.integer(model.week)
  wday <- match.wday.w(wday)
  return (yearWeekWdayDFToDate(cbind(seasonModelWeekToYearWeekDF(season, model.week, first.week, owning.wday), wday=wday), first.wday, owning.wday, error.on.wrap=error.on.wrap))
}

##' Convert season-model.week-wday \code{data.frame} to \code{Date}
##'
##' @param seasonModelWeekWday \code{data.frame} with \code{$season},
##'   \code{model.week}, and \code{wday} columns
##' @param first.week first week number of season
##' @param first.wday first weekday number of week
##' @param owning.wday integer-valued vector of weekday numbers: weekday that
##'   determines the year to which a week will be assigned
##' @param error.on.wrap TRUE or FALSE: if TRUE, an error is generated if a
##'   given week falls outside the range of possible week numbers for the
##'   corresponding year (i.e., if it is nonpositive or greater than the number
##'   of weeks assigned to the corresponding year in \code{year}); if FALSE, the
##'   week number will wrap around to future or previous years (e.g., week 0 of
##'   1997 will be considered the last week of 1996)
##'
##' @return \code{Date} vector with the corresponding dates
##'
##' @export
seasonModelWeekWdayDFToDate = function(seasonModelWeekWday, first.week, first.wday, owning.wday, error.on.wrap=TRUE) {
  return (seasonModelWeekWdayToDate(seasonModelWeekWday$season, seasonModelWeekWday$model.week, seasonModelWeekWday$wday, first.week, first.wday, owning.wday, error.on.wrap=error.on.wrap))
}

##' \code{first.wday} and \code{owning.wday} for some week numbering conventions
##'
##' Covers four common week numbering conventions:
##' \itemize{
##'  \item{"epi": }{Epidemiological weeks or "epi weeks": weeks begin on Sunday,
##'  and are assigned to years based on what year the majority of days fall in
##'  (i.e., what year Wednesday falls in)}
##'  \item{"iso": }{ISO 8601 weeks: weeks begin on Monday, and are assigned to
##'  years based on what year the majority of days fall in (i.e., what year
##'  Thursday falls in)}
##'  \item{"usa": }{USA convention: weeks begin on Sunday, and are assigned to
##'  years based on what year Sunday falls in}
##'  \item{"uk": }{UK convention: weeks begin on Monday, and are assigned to years
##'  based on what year Monday falls in}
##' }
##'
##' There are two rows, named \code{"first.wday"} and \code{"owning.wday"}.
##' There are four columns, corresponding to the four conventions above.
##'
##' @export
weekConventions = cbind(epi=c(0L,3L), iso=c(1L,4L), usa=c(0L,0L), uk=c(1L,1L))
rownames(weekConventions) <- c("first.wday","owning.wday")

## todo combine manual entries for time conversion functions?
## todo copy or provide link to examples in each method?
## todo S3-ify time conversion?  improve interface
## todo week convention strings as alternatives to specifying first&owning wday
## todo options interface for week convention, season first.week
## todo function to ensure rep_len compatibility / array conformability between args
