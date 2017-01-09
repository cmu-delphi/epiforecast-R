##
## You should have received a copy of the GNU General Public License
## along with epiforecast.  If not, see <http://www.gnu.org/licenses/>.
## license_header end

##' @include utils.R
##' @include match.R
##' @include weeks.R
##' @include delphi_epidata.R
##' @import httr
NULL

##' Epiweek before which no data should exist
##'
##' Assume that the universe was created in 1234 AD EW01, and that no data will
##' be from times before this.
firstEpiweekOfUniverse = 123401L

##' Augment df with year/week/season/model.week, make weekly, full seasons
##'
##' Given a df with with either (a) $year and $week, or (b) $date, fills in (if
##' missing) $year, $week, $date, $season, and $model.week. Fills in missing
##' weekly data from all seasons so that each season in \code{df$season} has all
##' of its model weeks in \code{df$model.week}.  Assumes epi week convention.
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
  } else if (all(c("year","week") %in% names(df))) {
    year = as.integer(df$year)
    week = as.integer(df$week)
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
  df <- do.call(base::data.frame, lapply(df, `[`, inds)) # df[inds,] w/ desired NA behavior a/f
  new.ywwd = DateToYearWeekWdayDF(new.date, 0L, 3L)
  new.year = new.ywwd$year
  new.week = new.ywwd$week
  new.season.model.week = yearWeekToSeasonModelWeekDF(new.year, new.week, first.week.of.season, 3L)
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
  has.default.rownames = identical(rownames(df), as.character(seq_len(nrow(df))))
  separate = split(df, df$season) # list of df's, one per season
  separate.prev = utils::head(separate, n=-1L)
  record.counts = sapply(separate.prev, function(season.df) {
    sum(!is.na(season.df[[signal.ind]]))
  })
  included = c(separate.prev[record.counts >= min.points.in.season], utils::tail(separate, n=1L))
  result = do.call(rbind, included)
  if (has.default.rownames)
    rownames(result) <- NULL # reset row names
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
##' @param cache.file single string; file path of Rdata file to create to cache
##'   the data from the delphi-epidata server, or \code{NULL} to not cache
##' @param cache.invalidation.period single \code{difftime}: time duration that
##'   must pass from last fetch for a new fetch to be performed instead of
##'   reading from the cache; default is one day
##' @param force.cache.invalidation single non-\code{NA} logical; if
##'   \code{TRUE}, then the \code{cache.invalidation.period} and a fetch will
##'   always be performed (the cache will not be read)
##'
##' @return \code{data.frame} with \code{$date}, \code{$year}, \code{$week}, and
##'   corresponding data (fields differ based on \code{source})
##'
##' @examples
##' ## All fluview data at the national level:
##' fluview.nat.all.df =
##'    trimPartialPastSeasons(fetchEpidataDF("fluview", "nat",
##'                           first.week.of.season=21L,
##'                           cache.file="fluview_nat_allfetch.Rdata"),
##'            "wili", min.points.in.season=33L)
##' ## Fluview data at the national level for more recent seasons for which data
##' ## was reported for all weeks (not just those in the influenza season), plus
##' ## the current season:
##' fluview.nat.recent.df =
##'    trimPartialPastSeasons(fetchEpidataDF("fluview", "nat",
##'                           first.week.of.season=21L,
##'                           cache.file="fluview_nat_allfetch.Rdata"),
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
      load(cache.file) # ==> fetch.info
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
    ## load(cache.file) # ==> fetch.info
    fetch.response = fetch.info$fetch.response
    fetch.time = fetch.info$fetch.response
  }

  if (!(fetch.response$result==1 && fetch.response$message=="success"))
    stop(sprintf('Failed to load data; result=%d, message="%s".', fetch.response$result, fetch.response$message))

  if (!is.null(cache.file) && should.fetch.now) {
    ## Update cache:
    fetch.info = list(fetch.time=fetch.time, fetch.response=fetch.response,
                      first.epiweek=first.epiweek, last.epiweek=last.epiweek)
    save(fetch.info, file=cache.file)
  }

  ## make list of lists a data.frame:
  ## df = do.call(rbind, lapply(fetch.response$epidata, as.data.frame)) # slower, doesn't handle NULL's
  df = as.data.frame(t(sapply(fetch.response$epidata, identity)))
  ## each column remains a list, even though it's a data.frame; fix this:
  df <- as.data.frame(lapply(df, function(col.as.list) {
    ## NULL's will be excluded when unlisted / do.call(c,.)'d; make them NA's first:
    col.as.list[sapply(col.as.list, is.null)] <- NA
    ## unlist(col.as.list) # slower
    do.call(c, col.as.list)
  }))
  if ("release_date" %in% names(df))
    df$release_date <- as.Date(df$release_date)

  df$year = as.integer(substr(df$epiweek, 1L, 4L))
  df$week = as.integer(substr(df$epiweek, 5L, 6L))

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
##' @param cache.file single string; file path of Rdata file to create to cache
##'   the data from the delphi-epidata server, or \code{NULL} to not cache
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
##'                       cache.file="fluview_nat_allfetch.Rdata")
##' ## Fluview data at the national level for more recent seasons for which data was
##' ## reported for all weeks (not just those in the influenza season), plus the
##' ## current season:
##' fluview.nat.recent.full.dat =
##'   fetchEpidataFullDat("fluview", "nat", "wili", 52L, first.week.of.season=21L,
##'                       cache.file="fluview_nat_allfetch.Rdata")
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
## fixme!!! restrict access to Twitter, GHT data?
## todo check fetching input, add default caching
