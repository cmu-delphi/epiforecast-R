
##' Test if Date is Christmas Day (vectorized)
##'
##' @param Date \code{Date} vector to test
##' @return logical vector
##' 
##' @export
is_christmas = function(Date) {
  posixlt = as.POSIXlt(Date)
  return (posixlt$mon + 1L == 12L & posixlt$mday == 25L)
}

##' Test if Date is (Gregorian) New Year's Day (vectorized)
##'
##' @param Date \code{Date} vector to test
##' @return logical vector
##'
##' @export
is_newyear = function(Date) {
  posixlt = as.POSIXlt(Date)
  return (posixlt$mon + 1L == 1L & posixlt$mday == 1L)
}

##' Test if Date is Thanksgiving Day (vectorized)
##'
##' @param Date \code{Date} vector to test
##' @return logical vector
##'
##' @export
is_thanksgiving = function(Date) {
  posixlt = as.POSIXlt(Date)
  return (posixlt$mon + 1L == 11L & (posixlt$mday - 1L) %/% 7L + 1L == 4L & posixlt$wday == 4L)
}

## todo: Chinese New Year; cannot use seasonal::cny since it is GPL3.0;
## investigate seasonal-cited source,
## http://www.chinesenewyears.info/chinese-new-year-calendar.php (license does
## not mention reuse) or ConvCalendar package (projectpluto.com calendars)

## alternative: timeDate::holiday / chron::is.holiday
