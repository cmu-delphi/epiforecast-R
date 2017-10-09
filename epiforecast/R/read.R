## author_header begin
## Copyright (C) 2016 Sangwon Hyun
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

##' Get index (time) of last non-\code{NA} in a vector or other new.dat.sim
##' @export
get.latest.time = function(new.dat.sim){
  new.dat.sim <- match.new.dat.sim(new.dat.sim)
    ys = new.dat.sim[["ys"]]

    is.na.ys = is.na(ys)

    if (!any(is.na.ys)) {
        stop("There are no missing values in new.dat.sim, i.e., the new season has been fully observed!")
    }

    ## check that NA occurrences match between all columns
    if (ncol(is.na.ys) > 1L && # preoptimization...: no need to check if only 1 trajectory
        any(xor(is.na.ys[,1], is.na.ys))) {
        stop ("All trajectories in new.dat.sim must have NA at exactly the same indices.")
    }

    ## Obtain index of non-NA
    time.of.forecast = max(0L, which(!is.na.ys[,1]))

    if (!all(time.of.forecast < which(is.na.ys[,1]))) {
        stop("new.dat.sim should be formatted so that non-NA's are followed by all NA's")
    }

    return(time.of.forecast)
}


##' Function that reads from a csv file
##' @param filename name of data file with each column equal to each season,
##'   with the first (n-1) columns to be, and the n'th column with the current
##'   season.
##' @export
read.from.file = function(filename){

    ## Sanity check of data file formatting
    check.file.contents(filename)

    ## Reformat table into list
    full.dat = read.csv(filename)
    return (table.to.list(full.dat))
}



##' Function to check if a data file is properly formatted (i.e. contains the
##' headers, is filled with numeric values, etc.) Current implementation is very
##' memory-inefficient.
##'
##' @param filename name of data filewith each column equal to each season, with
##'   the first (n-1) columns to be, and the n'th column with the current
##'   season.
check.file.contents = function(filename) {
  my.full.dat = read.csv(filename)
  check.table.format(my.full.dat)
}




##' Performs various checks on full.dat as a table (matrix/data.frame).
check.table.format = function(full.dat){

    ## Check if data frame
    stopifnot(class(full.dat)%in%c("matrix","data.frame"))

    ## Check if any column names are missing or NA
    default.colnames = Map(paste0, rep("V",ncol(full.dat)), 1:ncol(full.dat))
    if(any(colnames(full.dat) == default.colnames)) stop("Must supply all column names!")
    if(any(is.na(colnames(full.dat)))) stop("Column names contains NA values!")
    if(any(colnames(full.dat)=="NA")) stop("Column names contains NA values!")
    
    ## Check if it contains all numeric values.
    all.is.numeric = all(apply(full.dat,2, function(mycol){all(is.numeric(mycol))}))
    if(!all.is.numeric){
        stop("Table must contain all numeric values!")
    }

    ## Helper function
    check.column.to.see.if.partially.observed = function(mycolumn){
        ## Continue scan until you find the first three consecutive missing values
        scan.boolean.vector.for.three.consecutive.trues = function(ind, boolean.vector.to.scan){
            if(ind+3 > length(boolean.vector.to.scan)){stop("ind is bigger than n-3")}
            return(all(boolean.vector.to.scan[ind:(ind+3)]))}
        ## Make function to deal with |is.na(mycolumn)|
        myscan = function(ind){scan.boolean.vector.for.three.consecutive.trues(ind, is.na(mycolumn))}
        max.ind.three.NA.starts = max(which(!sapply(1:(length(mycolumn)-3), myscan)))
        return(max.ind.three.NA.starts)
    }

    ## Check if last column is partially observed
    last.column = full.dat[,ncol(full.dat)]
    max.ind.three.NA.starts = check.column.to.see.if.partially.observed(last.column)
    if(max.ind.three.NA.starts >= length(last.column) - 3){
        stop("I roughly scanned the last column of your table, but it doesn't seem to be partially observed!")
    }
}



##' Performs various checks on full.dat as a list of numeric vectors.
check.list.format = function(full.dat){
    ## Check if list is formatted correctly.
    if(!class(full.dat)%in%c("list")) stop("Type of input is not list!")

    ## Check if all numeric values
    all.is.numeric = all(sapply(full.dat, function(mycol){all(is.numeric(mycol))}))
    if(!all.is.numeric){
        stop("All vectors must only contain numeric values!")
    }

    ## Check if any column names are missing or NA
    default.colnames = Map(paste0, rep("V",length(full.dat)), 1:length(full.dat))
    if(any(is.na(names(full.dat)))) stop("Column names contains NA values!")
    if(any(names(full.dat) == default.colnames)) stop("Must supply all column names!")
    if(any(nchar(names(full.dat))==0)) stop("Some column names are missing!")

    ## ## Check partially observed last season
    ## if(!any(is.na(full.dat[[length(full.dat)]]))) stop("Last element (vector) of dataframe is not partially observed! (i.e. there should be some NA values!)")


    ## Check if last column is partially observed.
    ## a = (which(!is.na(last.column.weird.full.dat[,ncol(last.column.weird.full.dat)])))
    ## a.diff = c(a[2:length(a)],NA) - c(a[1:(length(a)-1)],NA) 
    ## if(any(a.diff[!is.na(a.diff)]!=1)) stop("Last column does not seem to be partially observed! i.e. no missing values!")
}



##' Function to change full.dat from table to list.
table.to.list = function(full.dat){
    mylist = lapply(seq_len(ncol(full.dat)), function(i) full.dat[,i])
    names(mylist) = colnames(full.dat)
    return(mylist)
}

  


## ##' Produces a |sim| object.
## ##' @param full.dat List of numeric vectors, each with proper names.
## ## ' @param area.name One of "nat" or "hhs*" where * is 1~10 #
## ## ' @param full.dat Either a properly formatted list (like the output of
## ## '     fetchepidataDF(); a list of numeric vectors, with the last vector being
## ## '     the 'new' season for which forecasts are being made), or a matrix/data
## ## '     frame with appropriate column names (and of which the last column
## ## '     represents the 'new' seaosn for which forecasts are being made.)
## ## ' @param filename The file name of the csv file that contains the table of
## ## '     data, with the first row being the names of the seasons, and each
## ## '     subsequent i'th row being the dataset for, say, week i. The last column
## ## '     should contain the 'new' season.
## ## ##' @examples 
## ## ##' sim = make.forecast(full.dat = NULL,
## ## ##'                     area.name = "hhs1",
## ## ##'                     first.model.week = 21L,
## ## ##'                     method = "br")
## ## ##' matplot(sim$ys[,1:100], type = 'l', col = 'cyan', lty=1)
## make.eb.forecast = function(full.dat,                
##                             first.model.week,
##                             min.points.in.season=52L,
##                             n.sim,
##                             eb.control.list = function(...){ eb.control.list(n.sim, ... ))){

##     ## Split into old dat (list) and new dat (vector)
##     old.dat = head(full.dat, -1L)
##     new.dat = tail(full.dat, 1L)[[1]]

##     time.of.forecast = get.latest.time(new.dat)

##     ## Unbundle control list that contains options for simulation
##     ## Then make forecast

##     if(is.null(eb.option)){
##         eb.control.list = get.eb.control.list()
##     } else {
##         option.names = names(eb.options)
##         eb.control.list = get.eb.control.list()
##     }

##     ## Make sim = eb.sim(dat, new.dat, control.list = eb.control.list)

##     ## Return the simulated curves
##     return(sim)
## }

