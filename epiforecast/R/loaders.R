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



##' @import tibble
##' @import httr
##' @import pipeR
##' @include utils.R
##' @include match.R
##' @include weeks.R
##' @include delphi_epidata.R
NULL

##' Epiweek before which no data should exist
##'
##' Assume that the universe was created in 1234 AD EW01, and that no data will
##' be from times before this.
firstEpiweekOfUniverse = 123401L

##' Augment df with epiweek/year/week/season/model.week, make weekly, full seasons
##'
##' Given a df with with either (a) $epiweek, (b) $year and $week, or (c) $date,
##' fills in (if missing) $epiweek, $year, $week, $date, $season, and
##' $model.week. Fills in missing weekly data from all seasons so that each
##' season in \code{df$season} has all of its model weeks in
##' \code{df$model.week}. Assumes epi week convention.
##'
##' Entries in \code{data.frame} are assumed without any checks to be sorted and
##' weekly (potentially with some skipped weeks).
##'
##' @param df \code{data.frame} with week numbers and other data
##' @param first.week.of.season the first week number in each season or
##'   \code{NULL} (the default); if \code{NULL}, then the first week of the
##'   season is assumed to be the week of the first data point.
##'
##' @export
augmentWeeklyDF = function(df, first.week.of.season=NULL) {
  if ("date" %in% names(df)) {
    date = df$date
    ywwd = DateToYearWeekWdayDF(date, 0L, 3L)
    year = ywwd$year
    week = ywwd$week
    if (length(unique(ywwd$wday)) > 1L) {
      stop("More than one day of the week is contained in =df$date=.")
    }
    wday = unique(ywwd$wday)
    if (length(wday)==0)
      wday <- 0L
    if ("year" %in% names(df) && !isTRUE(all.equal(df$year, year))) {
      stop("=df$date= appears not to match with =df$year=")
    }
    if ("week" %in% names(df) && !isTRUE(all.equal(df$week, week))) {
      stop("=df$date= appears not to match with =df$week=")
    }
    df$year <- year
    df$week <- week
  } else if ("epiweek" %in% names(df) || all(c("year","week") %in% names(df))) {
    if ("epiweek" %in% names(df)) {
      df$year <- as.integer(substr(df$epiweek, 1L, 4L))
      df$week <- as.integer(substr(df$epiweek, 5L, 6L))
    }
    year = as.integer(df$year)
    week = as.integer(df$week)
    df$epiweek <- 100L*year + week
    df$year <- year
    df$week <- week
    wday = 0L
    date = yearWeekWdayVecsToDate(year, week, wday, 0L,3L)
    df$date <- date
  } else {
    stop("=df= must contain either (a) $year and $week, or (b) $date.")
  }
  if (is.null(first.week.of.season))
    first.week.of.season <- utils::head(week, n=1L)
  ## season.model.week = yearWeekToSeasonModelWeekDF(year, week, first.week.of.season, 3L)
  season.model.week = yearWeekDFToSeasonModelWeekDF(df, first.week.of.season, 3L)
  ## if ("season" %in% names(df) && !isTRUE(all.equal(df$season, season.model.week$season))) {
  ##     stop("=df$season= appears not to match with other timing info")
  ## }
  ## if ("model.week" %in% names(df) && !isTRUE(all.equal(df$model.week, season.model.week$model.week))) {
  ##     stop("=df$model.week= appears not to match with other timing info")
  ## }
  df$season <- season.model.week$season
  df$model.week <- season.model.week$model.week

  ## fill in missing weekly data points:
  last.season = utils::tail(df$season, 1L)
  ## season.dates.but.last = seasonDates(utils::head(df$season, 1L), last.season-1, first.week.of.season, 0L,3L)
  ## season.dates.last = structure(
  ##     list(df$date[df$season==last.season]),
  ##     names=names(namedSeason(last.season, first.week.of.season)))
  ## season.dates = c(season.dates.but.last, season.dates.last)
  season.dates = DatesOfSeason(Seq(utils::head(df$season, 1L), last.season), first.week.of.season, 0L, 3L)
  new.date = do.call(base::c, season.dates) # unlist(season.dates) w/o class change
  inds = match(new.date, date)
  df <- do.call(tibble::data_frame, lapply(df, `[`, inds)) # df[inds,] w/ desired NA behavior a/f
  new.ywwd = DateToYearWeekWdayDF(new.date, 0L, 3L)
  new.epiweek = 100L*new.ywwd$year + new.ywwd$week
  new.year = new.ywwd$year
  new.week = new.ywwd$week
  new.season.model.week = yearWeekToSeasonModelWeekDF(new.year, new.week, first.week.of.season, 3L)
  df$epiweek <- new.epiweek
  df$year <- new.year
  df$week <- new.week
  df$season <- new.season.model.week$season
  df$model.week <- new.season.model.week$model.week
  df$date <- yearWeekWdayVecsToDate(new.year, new.week, wday, 0L, 3L)
  return (df)
}

##' Trims incomplete past seasons from a \code{data.frame}
##'
##' Removes rows from \code{df} corresponding to "past" seasons (i.e., all but
##' the last season in \code{df}) for which \code{df} has less than
##' \code{min.points.in.season} non-missing entries in \code{df[[signal.ind]]}.
##'
##' @param df data frame with columns \code{df$season}, and
##'   \code{df[[signal.ind]]}
##' @param signal.ind single non-NA character/integer-valued index for column of
##'   \code{df}
##' @param min.points.in.season the minimum number of non-\code{NA} values for
##'   \code{signal.ind} that a season must have in order to be retained; all
##'   rows corresponding to seasons containing less observations will be removed
##'   from \code{df}
##'
##' @export
trimPartialPastSeasons = function(df, signal.ind, min.points.in.season) {
  if (!is.character(signal.ind) && !is.numeric(signal.ind)) {
    stop("signal.ind must satisfy is.character or is.numeric")
  } else if (length(signal.ind)!=1) {
    stop("signal.ind must have length 1")
  } else if (is.na(signal.ind)) {
    stop("signal.ind must not be NA")
  } else if (is.numeric(signal.ind)) {
    signal.ind <- match.single.nonna.integer(signal.ind)
  }
  min.points.in.season <- match.integer(min.points.in.season)
  separate = split(df, df$season) # list of df's, one per season
  separate.prev = utils::head(separate, n=-1L)
  record.counts = sapply(separate.prev, function(season.df) {
    sum(!is.na(season.df[[signal.ind]]))
  })
  included = c(separate.prev[record.counts >= min.points.in.season], utils::tail(separate, n=1L))
  result = dplyr::bind_rows(included)
  return (result)
}

##' Fetch & cache \href{https://github.com/undefx/delphi-epidata}{delphi-epidata} data, convert it to a \code{data.frame}
##'
##' @param source length-1 character vector; name of data source; one of
##'   \itemize{
##'
##'   \item{\code{"fluview"}: }{U.S. Outpatient Influenza-like Illness Surveillance
##'  Network (ILINet) data from CDC FluView reports; weekly, national or by HHS
##'  region; updated weekly}
##'
##'   \item{\code{"ilinet"}: }{Estimates of ILINet data at state level; weekly; static}
##'
##'   \item{\code{"gft"}: }{Google Flu Trends data}
##'
##'   \item{\code{"ght"}: }{Google Health Trends data; currently restricted, not
##'   supported by this function}
##'
##'   \item{\code{"twitter"}: }{Twitter data; currently restricted, not supported by this
##'   function}
##'
##'   \item{\code{"wiki"}: }{Wikipedia access log data; currently not supported by this
##'   function}
##'
##'   \item{\code{"nidss.flu"}: }{Taiwan National Infectious Disease Statistics System
##'   (NIDSS) outpatient ILI data}
##'
##'   \item{\code{"nidss.dengue"}: }{Taiwan National Infectious Disease
##'   Statistics System (NIDSS) dengue incidence data}
##'
##' }
##'
##' @param area length-1 character vector; name of area; possibilities by
##'   \code{source} are: \itemize{
##'
##'   \item{\code{"fluview"}: }{\code{"hhs1"}, \code{"hhs2"}, ...,
##'   \code{"hhs10"}, \code{"nat"}}
##'
##'   \item{\code{"ilinet"}: }{2-character state abbreviation or \code{"DC"}}
##'
##'   \item{\code{"gft"}: }{listed at
##'   \url{https://github.com/undefx/delphi-epidata#gft-parameters}}
##'
##'   \item{\code{"nidss.flu"}, \code{"nidss.dengue"}: }{one of
##'   \code{"nationwide"}, \code{"central"}, \code{"eastern"}, \code{"kaoping"},
##'   \code{"northern"}, \code{"southern"}, \code{"taipei"}} }
##'
##' @param lag single integer value or \code{NULL}; for supported data sources
##'   for which observations for a particular time are revised at later times as
##'   more data is received, gives access to some of the older versions of the
##'   data. \code{NULL} gives the latest revision of the observation. An integer
##'   value gives the value for each observation \code{lag} weeks after its
##'   initial report; the value \code{0} corresponds to using the initial report
##'   for each data point. (This function does not reconstruct what a the data
##'   looked like at any particular point in the past, but can be used to do
##'   so.)
##' @param first.week.of.season single integer value giving the week number that
##'   seasons start with, or \code{NULL} to have seasons start with the week
##'   number of the first data point
##' @param first.epiweek year-epiweek as string in form \code{"YYYYww"}
##'   specifying a lower limit on times for which data will be fetched, or
##'   \code{NULL} to not impose a lower limit
##' @param last.epiweek year-epiweek as string in form \code{"YYYYww"}
##'   specifying an upper limit on times for which data will be fetched, or
##'   \code{NULL} to not impose an upper limit
##' @param cache.prefix length-1 \code{character}; file path of \code{.rds} file to
##'   create to cache the data from the delphi-epidata server, or \code{NULL} to
##'   not cache
##' @param cache.invalidation.period single \code{difftime}: time duration that
##'   must pass from last fetch for a new fetch to be performed instead of
##'   reading from the cache; default is one day
##' @param force.cache.invalidation single non-\code{NA} logical; if
##'   \code{TRUE}, then the \code{cache.invalidation.period} and a fetch will
##'   always be performed (the cache will not be read)
##'
##' @return data frame (\code{tbl_df}) with \code{$date}, \code{$year},
##'   \code{$week}, and corresponding data (fields differ based on
##'   \code{source})
##'
##' @examples
##' ## All fluview data at the national level:
##' fluview.nat.all.df =
##'    trimPartialPastSeasons(fetchEpidataDF("fluview", "nat",
##'                           first.week.of.season=21L,
##'                           cache.file="fluview_nat_allfetch.rds"),
##'            "wili", min.points.in.season=33L)
##' ## Fluview data at the national level for more recent seasons for which data
##' ## was reported for all weeks (not just those in the influenza season), plus
##' ## the current season:
##' fluview.nat.recent.df =
##'    trimPartialPastSeasons(fetchEpidataDF("fluview", "nat",
##'                           first.week.of.season=21L,
##'                           cache.file="fluview_nat_allfetch.rds"),
##'            "wili", min.points.in.season=52L)
##' @export
fetchEpidataDF = function(source, area, lag=NULL,
                          first.week.of.season=NULL,
                          first.epiweek=NULL, last.epiweek=NULL,
                          cache.file=NULL, cache.invalidation.period=as.difftime(1L, units="days"), force.cache.invalidation=FALSE) {
  if (!is.null(lag) && length(lag) > 1) stop("Fetching epidata for multiple lags simultaneously is not supported.")
  ## todo drop first season if not enough data for specified first week of season

  if (is.null(first.epiweek)) first.epiweek <- firstEpiweekOfUniverse
  if (is.null(last.epiweek)) last.epiweek <- as.integer(format(Sys.Date(),"%Y53"))

  should.fetch.now =
    if (is.null(cache.file)) { TRUE
    } else if (force.cache.invalidation) { TRUE
    } else if (!file.exists(cache.file)) { TRUE
    } else {
      fetch.info = readRDS(cache.file) # ==> fetch.info
      should.refetch =
        difftime(Sys.time(), fetch.info$fetch.time) > cache.invalidation.period ||
        ## !(fetch.info$fetch.response$result==1 && fetch.info$fetch.response$message=="success") ||
        first.epiweek != fetch.info$first.epiweek || last.epiweek != fetch.info$last.epiweek
      ## print(difftime(Sys.time(), file.info(cache.file)$mtime))
      ## should.refetch = difftime(Sys.time(), file.info(cache.file)$mtime) > cache.invalidation.period
      should.refetch
    }

  if (should.fetch.now) {
    message("No cache, empty cache, expired cache, or forced cache invalidation; fetching data from server.")
    ## todo prompt for confirmation on fetch/refetch
    ## todo option to read from cache without considering refetch
    if (is.null(lag)) {
      fetch.response = Epidata[[source]](area, Epidata$range(first.epiweek, last.epiweek))
    } else {
      fetch.response = Epidata[[source]](area, Epidata$range(first.epiweek, last.epiweek), lag=lag)
    }
    fetch.time = Sys.time()
  } else {
    message("Cached version of data used.")
    fetch.response = fetch.info$fetch.response
    fetch.time = fetch.info$fetch.response
  }

  if (!(fetch.response$result==1 && fetch.response$message=="success"))
    stop(sprintf('Failed to read data; result=%d, message="%s".', fetch.response$result, fetch.response$message))

  if (!is.null(cache.file) && should.fetch.now) {
    ## Update cache:
    fetch.info = list(fetch.time=fetch.time, fetch.response=fetch.response,
                      first.epiweek=first.epiweek, last.epiweek=last.epiweek)
    saveRDS(fetch.info, cache.file)
  }

  ## make list of lists a data.frame:
  ## df = do.call(rbind, lapply(fetch.response$epidata, as.data.frame)) # slower, doesn't handle NULL's
  df = tibble::as_data_frame(t(sapply(fetch.response$epidata, identity)))
  ## each column remains a list, even though it's a data.frame; fix this:
  df <- tibble::as_data_frame(lapply(df, function(col.as.list) {
    ## NULL's will be excluded when unlisted / do.call(c,.)'d; make them NA's first:
    col.as.list[sapply(col.as.list, is.null)] <- NA
    ## unlist(col.as.list) # slower
    do.call(c, col.as.list)
  }))
  if ("release_date" %in% names(df))
    df$release_date <- as.Date(df$release_date)

  df <- augmentWeeklyDF(df, first.week.of.season)

  attr(df, "source") <- source
  attr(df, "area") <- area

  return (df)
}

##' Fetch data from delphi-epidata, trim partial seasons, and convert to list of trajectories
##' 
##' @param source length-1 character vector; name of data source; one of
##'   \itemize{
##'
##'   \item{\code{"fluview"}: }{U.S. Outpatient Influenza-like Illness Surveillance
##'  Network (ILINet) data from CDC FluView reports; weekly, national or by HHS
##'  region; updated weekly}
##'
##'   \item{\code{"ilinet"}: }{Estimates of ILINet data at state level; weekly; static}
##'
##'   \item{\code{"gft"}: }{Google Flu Trends data}
##'
##'   \item{\code{"ght"}: }{Google Health Trends data; currently restricted, not
##'   supported by this function}
##'
##'   \item{\code{"twitter"}: }{Twitter data; currently restricted, not supported by this
##'   function}
##'
##'   \item{\code{"wiki"}: }{Wikipedia access log data; currently not supported by this
##'   function}
##'
##'   \item{\code{"nidss.flu"}: }{Taiwan National Infectious Disease Statistics System
##'   (NIDSS) outpatient ILI data}
##'
##'   \item{\code{"nidss.dengue"}: }{Taiwan National Infectious Disease
##'   Statistics System (NIDSS) dengue incidence data}
##'
##' }
##'
##' @param area length-1 character vector; name of area; possibilities by
##'   \code{source} are: \itemize{
##'
##'   \item{\code{"fluview"}: }{\code{"hhs1"}, \code{"hhs2"}, ...,
##'   \code{"hhs10"}, \code{"nat"}}
##'
##'   \item{\code{"ilinet"}: }{2-character state abbreviation or \code{"DC"}}
##'
##'   \item{\code{"gft"}: }{listed at
##'   \url{https://github.com/undefx/delphi-epidata#gft-parameters}}
##'
##'   \item{\code{"nidss.flu"}, \code{"nidss.dengue"}: }{one of
##'   \code{"nationwide"}, \code{"central"}, \code{"eastern"}, \code{"kaoping"},
##'   \code{"northern"}, \code{"southern"}, \code{"taipei"}}}
##'
##' @param signal.ind single non-NA character/integer-valued index for column of
##'   \code{df}
##' @param min.points.in.season the minimum number of non-\code{NA} values for
##'   \code{signal.ind} that a season must have in order to be retained; all
##'   rows corresponding to seasons containing less observations will be removed
##'   from \code{df}
##' @param lag single integer value or \code{NULL}; for supported data sources
##'   for which observations for a particular time are revised at later times as
##'   more data is received, gives access to some of the older versions of the
##'   data. \code{NULL} gives the latest revision of the observation. An integer
##'   value gives the value for each observation \code{lag} weeks after its
##'   initial report; the value \code{0} corresponds to using the initial report
##'   for each data point. (This function does not reconstruct what a the data
##'   looked like at any particular point in the past, but can be used to do
##'   so.)
##' @param first.week.of.season single integer value giving the week number that
##'   seasons start with, or \code{NULL} to have seasons start with the week
##'   number of the first data point
##' @param first.epiweek year-epiweek as string in form \code{"YYYYww"}
##'   specifying a lower limit on times for which data will be fetched, or
##'   \code{NULL} to not impose a lower limit
##' @param last.epiweek year-epiweek as string in form \code{"YYYYww"}
##'   specifying an upper limit on times for which data will be fetched, or
##'   \code{NULL} to not impose an upper limit
##' @param cache.file single string; file path of \code{.rds} file to create to
##'   cache the data from the delphi-epidata server, or \code{NULL} to not cache
##' @param cache.invalidation.period single \code{difftime}: time duration that
##'   must pass from last fetch for a new fetch to be performed instead of
##'   reading from the cache; default is one day
##' @param force.cache.invalidation single non-\code{NA} logical; if
##'   \code{TRUE}, then the \code{cache.invalidation.period} and a fetch will
##'   always be performed (the cache will not be read)
##'
##' @return named list of 52/53-length is.numeric vectors; each vector is the
##'   trajectory of the given signal for a single season, with NA's used to fill
##'   in for missing and future data; the names are of the form "SYYYY", where
##'   YYYY is the first year of the season, e.g., "S2003" corresponds to the
##'   2003--2004 season.
##'
##' @examples
##' ## All fluview data at the national level:
##' fluview.nat.all.full.dat =
##'   fetchEpidataFullDat("fluview", "nat", "wili", 33L, first.week.of.season=21L,
##'                       cache.file="fluview_nat_allfetch.rds")
##' ## Fluview data at the national level for more recent seasons for which data was
##' ## reported for all weeks (not just those in the influenza season), plus the
##' ## current season:
##' fluview.nat.recent.full.dat =
##'   fetchEpidataFullDat("fluview", "nat", "wili", 52L, first.week.of.season=21L,
##'                       cache.file="fluview_nat_allfetch.rds")
##'
##' @export
fetchEpidataFullDat = function(source,
                               area,
                               signal.ind,
                               min.points.in.season=52L,
                               lag=NULL,
                               first.week.of.season=NULL,
                               first.epiweek=NULL,
                               last.epiweek=NULL,
                               cache.file=NULL,
                               cache.invalidation.period=as.difftime(1L, units="days"),
                               force.cache.invalidation=FALSE) {
  df = trimPartialPastSeasons(
    fetchEpidataDF(source, area, lag=lag,
                   first.week.of.season=first.week.of.season,
                   first.epiweek=first.epiweek, last.epiweek=last.epiweek,
                   cache.file=cache.file, cache.invalidation.period=cache.invalidation.period),
    signal.ind, min.points.in.season=min.points.in.season)
  full.dat = split(df[[signal.ind]], df$season)
  names(full.dat) <- sprintf("S%s", names(full.dat))
  return (full.dat)
}

##' Combine current and lagged epidata into a data frame
##'
##' Take the epidata data frames for the specified lags plus the current
##' (unlagged) epidata data frame and combine them together into a single data
##' frame (\code{tbl_df}).
##'
##' @seealso fetchEpidataHistoryDT
##'
##' @export
fetchEpidataHistoryDF = function(source, area, lags,
                                 first.week.of.season=NULL,
                                 first.epiweek=NULL, last.epiweek=NULL,
                                 cache.file.prefix, cache.invalidation.period=as.difftime(1L, units="days"), force.cache.invalidation=FALSE) {
  lag.dfs = lapply(lags, function(lag) {
    fetchEpidataDF(source, area, lag,
                   first.week.of.season=first.week.of.season,
                   first.epiweek=first.epiweek, last.epiweek=last.epiweek,
                   cache.file=paste0(cache.file.prefix,"_lag",lag,".rds"), cache.invalidation.period=cache.invalidation.period, force.cache.invalidation=force.cache.invalidation)
  })
  current.df = fetchEpidataDF(source, area, NULL,
                   first.week.of.season=first.week.of.season,
                   first.epiweek=first.epiweek, last.epiweek=last.epiweek,
                   cache.file=paste0(cache.file.prefix,"_current.rds"), cache.invalidation.period=cache.invalidation.period, force.cache.invalidation=force.cache.invalidation)
  history.df = dplyr::bind_rows(dplyr::bind_rows(lag.dfs), current.df) %>>%
    dplyr::distinct()
  return (history.df)
}

##' Combine current and lagged epidata into a \code{data.table}
##'
##' Data table form of \code{\link{fetchEpidataHistoryDF}}; sets the keys as
##' \code{epiweek} and \code{issue} and adds a column \code{actual.issue}
##' duplicating \code{issue}, to be used in rolling joins used in one possible
##' implementation of \code{\link{mimicPastEpidataDF}}.
##'
##' @export
fetchEpidataHistoryDT = function(source, area, lags,
                                 first.week.of.season=NULL,
                                 first.epiweek=NULL, last.epiweek=NULL,
                                 cache.file.prefix, cache.invalidation.period=as.difftime(1L, units="days"), force.cache.invalidation=FALSE) {
  history.df = fetchEpidataHistoryDF(source, area, lags,
                                     first.week.of.season=first.week.of.season,
                                     first.epiweek=first.epiweek, last.epiweek=last.epiweek,
                                     cache.file.prefix=cache.file.prefix, cache.invalidation.period=cache.invalidation.period, force.cache.invalidation=force.cache.invalidation)
  history.dt = data.table::as.data.table(history.df)
  data.table::setkey(history.dt, epiweek, issue)
  history.dt[,actual.issue:=issue]
  return (history.dt)
}

mimicPastEpidataDF1 = function(history.dt.or.df, forecast.epiweek) {
  history.df = tibble::as_data_frame(history.dt.or.df)
  available.and.recorded.issues =
    history.df %>>%
    dplyr::select(epiweek, issue) %>>%
    dplyr::filter(issue <= forecast.epiweek) %>>%
    dplyr::group_by(epiweek) %>>%
    dplyr::summarize(issue=max(issue, na.rm=TRUE)) %>>%
    dplyr::ungroup()
  future.fillin.for.missing.issues =
    history.df %>>%
    dplyr::select(epiweek, issue) %>>%
    dplyr::filter(epiweek <= forecast.epiweek &
                  issue > forecast.epiweek) %>>%
    dplyr::group_by(epiweek) %>>%
    dplyr::summarize(issue=min(issue, na.rm=TRUE)) %>>%
    dplyr::ungroup()
  available.and.recorded.issues %>>%
    dplyr::bind_rows(future.fillin.for.missing.issues) %>>%
    dplyr::group_by(epiweek) %>>%
    dplyr::filter(issue==min(issue, na.rm=TRUE)) %>>%
    dplyr::ungroup() %>>%
    dplyr::arrange(epiweek) %>>%
    dplyr::left_join(history.df, c("epiweek","issue")) %>>%
    ## remove near-duplicate records (lag=NULL vs not)
    dplyr::group_by(epiweek) %>>%
    dplyr::filter(seq_along(epiweek)==1L) %>>%
    dplyr::ungroup() %>>%
    dplyr::mutate(forecast.epiweek=forecast.epiweek) %>>%
    augmentWeeklyDF(history.df[["week"]][1L]) %>>%
    {.}
}

mimicPastEpidataDF2 = function(history.dt, forecast.epiweek) {
  timing.dt = unique(history.dt[epiweek <= forecast.epiweek, list(epiweek)])
  timing.dt[,issue:=forecast.epiweek]
  available.and.recorded = history.dt[timing.dt,,roll=Inf]
  future.fillin.for.missing = history.dt[timing.dt,,roll=-Inf]
  dplyr::bind_rows(tibble::as_data_frame(available.and.recorded),
                   tibble::as_data_frame(future.fillin.for.missing)) %>>%
    dplyr::rename(forecast.epiweek=issue) %>>%
    dplyr::rename(issue=actual.issue) %>>%
    dplyr::arrange(epiweek) %>>%
    dplyr::group_by(epiweek) %>>%
    dplyr::filter(issue==min(issue, na.rm=TRUE)) %>>%
    ## remove near-duplicate records (lag=NULL vs not)
    dplyr::filter(seq_along(issue)==1L) %>>%
    dplyr::ungroup() %>>%
    augmentWeeklyDF(history.dt[1L,week]) %>>%
    {.}
}

##' Mimic what fetchEpidataDF would have given in the past
##'
##' Attempt to reproduce what fetchEpidataDF would have produced given data
##' through issue \code{forecast.epiweek} (i.e., data through the time when the
##' initial report for \code{forecast.epiweek} would have been available).
##' Historical data is not always complete. Try the following in order for each
##' observation week: (1) look for data for the report corresponding to the
##' forecast week, (2) look for data in earlier reports (taking the latest
##' report with data), (3) look for data in later reports (taking the earliest
##' report with data). For example, the ILINet (\code{source="fluview"}) US
##' National (\code{area="nat"}) historical data skips issue 200352; when
##' \code{forecast.epiweek} is chosen as 200352, observations for \code{epiweek}
##' 200340 through 200351 are provided by the previous issue, 200351, while the
##' observation for \code{epiweek} 200352 uses data from the \emph{next} issue,
##' 200353. (Observations for earlier seasons also use "future" data from issue
##' 201352, because issue 201352 is special, filling in data for epiweeks
##' missing from all available past issues.) Isolated skipped weeks are
##' uncommon; usually when backfill data is missing, it is for earlier seasons,
##' during the off-season, or at the beginning of a season; for example, in HHS
##' Region 1 (\code{area="hhs1"}), the first recorded issue was 200949, so
##' mimicing any reports from the 2005/2006 season will use finalized data from
##' issue 201352 instead, and mimicing report 200940 will fill in 200940 with
##' data from issue 200949. Similarly, data is missing from the following
##' off-season through issue 201040, inclusive, so mimicing report 201040 will
##' fill in off-season data with finalized data from issue 201352, and the
##' observation for 201040 with data from issue 201041.
##'
##' @examples
##'
##' ## set up a cache directory:
##' epidata.cache.dir = "~/.epiforecast-cache"
##' if (!dir.exists(epidata.cache.dir)) {
##'   dir.create(epidata.cache.dir)
##' }
##' ## fetch HHS Region 1 ILINet version history, assuming no significant
##' ## revisions past lag 51:
##' fluview.area = "hhs1"
##' history.df = fetchEpidataHistoryDF(
##'   "fluview", fluview.area, 0:51,
##'   first.week.of.season = 31L,
##'   cache.file.prefix=file.path(epidata.cache.dir,paste0("fluview_",fluview.area))
##' )
##' ## Simulate fetchEpidataHistoryDF when issue 201540 (the first ILINet report
##' ## containing data for epiweek 201540) just came out:
##' mimicPastEpidataDF(history.df, 201540L)
##'
##' @export
mimicPastEpidataDF = mimicPastEpidataDF1

## ## todo turn into test
## history.dt = fetchEpidataHistoryDT("fluview", "hhs1", 0:51,
##                            first.week.of.season = 31L,
##                            cache.file.prefix="~/.epiforecast-cache/fluview_hhs1")
## list(mimicPastEpidataDF1, mimicPastEpidataDF2) %>>%
##   lapply(function(mimicPastEpidataDFn) {
##     ## mimicPastEpidataDFn(history.dt, 201540L) %>>%
##     mimicPastEpidataDFn(history.dt, 201040L) %>>%
##       dplyr::arrange(-epiweek) %>>%
##       dplyr::select(epiweek, issue, forecast.epiweek, wili)
##   }) %>>%
##   do.call(what=identical) %>>%
##   {.}

fit.to.oldfit = function(fit) {
  f = lapply(fit, `[[`, "f")
  f <- sapply(f, `[`, seq_len(min(lengths(f))))
  tau = sapply(fit, `[[`, "tau")
  return (list(f=f, tau=tau))
}

oldfit.to.fit = function(oldfit, type="Gaussian") {
  return (lapply(seq_along(oldfit$tau), function(fit.s.i) {
    list(f=oldfit$f[,fit.s.i], tau=oldfit$tau[fit.s.i], type=type)
  }))
}

smooth.curves.to.fit = function(smooth.curves, type="Gaussian") {
  return (lapply(seq_along(smooth.curves$sigma.hat), function(fit.s.i) {
    list(f=smooth.curves$smooth.obj[[fit.s.i]], tau=smooth.curves$sigma.hat[fit.s.i], type=type)
  }))
}

## todo fail gracefully with curl errors in Epidata
## todo check augmentWeeklyDF input
## todo trimPartialPastSeasons setting to trim all incomplete (like min # being # of weeks in season)
## todo check fetching input, add default caching
