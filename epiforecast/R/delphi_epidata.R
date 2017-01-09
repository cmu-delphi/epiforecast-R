## This file is redistributed from https://github.com/undefx/delphi-epidata with modification; below is the original license for https://github.com/undefx/delphi-epidata.
## The only modification is the addition of namespace qualifiers to some function calls.
## 
## Redistribution license information (note disclaimer of warranty):
## The MIT License (MIT)
## 
## Copyright (c) 2015
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
## FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS,
## COPYRIGHT HOLDERS, OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
## ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
## 
## Original license information:
## The MIT License (MIT)
## 
## Copyright (c) 2015 
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.
## 
# A module for DELPHI's Epidata API.
#
# https://github.com/undefx/delphi-epidata
#
# Notes:
#  - Requires the `httr` library.

# External libraries
library(httr)


# Because the API is stateless, the Epidata class only contains static methods
Epidata <- (function() {

  # API base url
  BASE_URL <- 'http://delphi.midas.cs.cmu.edu/epidata/api.php'

  # Helper function to cast values and/or ranges to strings
  .listitem <- function(value) {
    if(is.list(value) && 'from' %in% names(value) && 'to' %in% names(value)) {
      return(paste0(toString(value$from), '-', toString(value$to)))
    } else {
      return(toString(value))
    }
  }

  # Helper function to build a list of values and/or ranges
  .list <- function(values) {
    if(!is.list(values) || ('from' %in% names(values) && 'to' %in% names(values))) {
      values <- list(values)
    }
    return(paste(sapply(values, .listitem), collapse=','))
  }

  # Helper function to request and parse epidata
  .request <- function(params) {
    # API call
    return(httr::content(httr::GET(BASE_URL, query=params), 'parsed'))
  }

  # Build a `range` object (ex: dates/epiweeks)
  range <- function(from, to) {
    if(to <= from) {
        temp <- to
        to <- from
        from <- temp
    }
    return(list(from=from, to=to))
  }

  # Fetch FluView data
  fluview <- function(regions, epiweeks, issues, lag) {
    # Check parameters
    if(missing(regions) || missing(epiweeks)) {
      stop('`regions` and `epiweeks` are both required')
    }
    if(!missing(issues) && !missing(lag)) {
      stop('`issues` and `lag` are mutually exclusive')
    }
    # Set up request
    params <- list(
      source = 'fluview',
      regions = .list(regions),
      epiweeks = .list(epiweeks)
    )
    if(!missing(issues)) {
      params$issues <- .list(issues)
    }
    if(!missing(lag)) {
      params$lag <- lag
    }
    # Make the API call
    return(.request(params))
  }

  # Fetch ILINet data
  ilinet <- function(locations, epiweeks) {
    # Check parameters
    if(missing(locations) || missing(epiweeks)) {
      stop('`locations` and `epiweeks` are both required')
    }
    # Set up request
    params <- list(
      source = 'ilinet',
      locations = .list(locations),
      epiweeks = .list(epiweeks)
    )
    # Make the API call
    return(.request(params))
  }

  # Fetch Google Flu Trends data
  gft <- function(locations, epiweeks) {
    # Check parameters
    if(missing(locations) || missing(epiweeks)) {
      stop('`locations` and `epiweeks` are both required')
    }
    # Set up request
    params <- list(
      source = 'gft',
      locations = .list(locations),
      epiweeks = .list(epiweeks)
    )
    # Make the API call
    return(.request(params))
  }

  # Fetch Google Health Trends data
  ght <- function(auth, locations, epiweeks, query) {
    # Check parameters
    if(missing(auth) || missing(locations) || missing(epiweeks) || missing(query)) {
      stop('`auth`, `locations`, `epiweeks`, and `query` are all required')
    }
    # Set up request
    params <- list(
      source = 'ght',
      auth = auth,
      locations = .list(locations),
      epiweeks = .list(epiweeks),
      query = query
    )
    # Make the API call
    return(.request(params))
  }

  # Fetch HealthTweets data
  twitter <- function(auth, locations, dates, epiweeks) {
    # Check parameters
    if(missing(auth) || missing(locations)) {
      stop('`auth` and `locations` are both required')
    }
    if(!xor(missing(dates), missing(epiweeks))) {
      stop('exactly one of `dates` and `epiweeks` is required')
    }
    # Set up request
    params <- list(
      source = 'twitter',
      auth = auth,
      locations = .list(locations)
    )
    if(!missing(dates)) {
      params$dates <- .list(dates)
    }
    if(!missing(epiweeks)) {
      params$epiweeks <- .list(epiweeks)
    }
    # Make the API call
    return(.request(params))
  }

  # Fetch Wikipedia access data
  wiki <- function(articles, dates, epiweeks, hours) {
    # Check parameters
    if(missing(articles)) {
      stop('`articles` is required')
    }
    if(!xor(missing(dates), missing(epiweeks))) {
      stop('exactly one of `dates` and `epiweeks` is required')
    }
    # Set up request
    params <- list(
      source = 'wiki',
      articles = .list(articles)
    )
    if(!missing(dates)) {
      params$dates <- .list(dates)
    }
    if(!missing(epiweeks)) {
      params$epiweeks <- .list(epiweeks)
    }
    if(!missing(hours)) {
      params$hours <- .list(hours)
    }
    # Make the API call
    return(.request(params))
  }

  # Fetch NIDSS flu data
  nidss.flu <- function(regions, epiweeks, issues, lag) {
    # Check parameters
    if(missing(regions) || missing(epiweeks)) {
      stop('`regions` and `epiweeks` are both required')
    }
    if(!missing(issues) && !missing(lag)) {
      stop('`issues` and `lag` are mutually exclusive')
    }
    # Set up request
    params <- list(
      source = 'nidss_flu',
      regions = .list(regions),
      epiweeks = .list(epiweeks)
    )
    if(!missing(issues)) {
      params$issues <- .list(issues)
    }
    if(!missing(lag)) {
      params$lag <- lag
    }
    # Make the API call
    return(.request(params))
  }

  # Fetch NIDSS dengue data
  nidss.dengue <- function(locations, epiweeks) {
    # Check parameters
    if(missing(locations) || missing(epiweeks)) {
      stop('`locations` and `epiweeks` are both required')
    }
    # Set up request
    params <- list(
      source = 'nidss_dengue',
      locations = .list(locations),
      epiweeks = .list(epiweeks)
    )
    # Make the API call
    return(.request(params))
  }

  # Fetch Delphi's forecast
  delphi <- function(system, epiweek) {
    # Check parameters
    if(missing(system) || missing(epiweek)) {
      stop('`system` and `epiweek` are both required')
    }
    # Set up request
    params <- list(
      source = 'delphi',
      system = system,
      epiweek = epiweek
    )
    # Make the API call
    return(.request(params))
  }

  # Fetch Delphi's digital surveillance signals
  signals <- function(auth, names, locations, epiweeks) {
    # Check parameters
    if(missing(auth) || missing(names) || missing(locations) || missing(epiweeks)) {
      stop('`auth`, `names`, `locations`, and `epiweeks` are all required')
    }
    # Set up request
    params <- list(
      source = 'signals',
      auth = auth,
      names = .list(names),
      locations = .list(locations),
      epiweeks = .list(epiweeks)
    )
    # Make the API call
    return(.request(params))
  }

  # Fetch Delphi's wILI nowcast
  nowcast <- function(locations, epiweeks) {
    # Check parameters
    if(missing(locations) || missing(epiweeks)) {
      stop('`locations` and `epiweeks` are both required')
    }
    # Set up request
    params <- list(
      source = 'nowcast',
      locations = .list(locations),
      epiweeks = .list(epiweeks)
    )
    # Make the API call
    return(.request(params))
  }

  # Export the public methods
  return(list(
    range = range,
    fluview = fluview,
    ilinet = ilinet,
    gft = gft,
    ght = ght,
    twitter = twitter,
    wiki = wiki,
    nidss.flu = nidss.flu,
    nidss.dengue = nidss.dengue,
    delphi = delphi,
    signals = signals,
    nowcast = nowcast
  ))
})()
