
library("pipeR")
source("covidplot.R")

## year = 2020L
## epi.week = 10L
## crowdcast.tgz.dirpath = "~/Dropbox/private/crowdcast_spreadsheets_to_filter"
## submissions.root.dirpath = "~/Dropbox/private/test_submissions"
## do.plot = TRUE
## plots.root.dirpath = "~/Dropbox/private/test_plots"
## create.destination.dirs = FALSE

## --- Provide usage message / Read in command-line arguments: ---
command.args = commandArgs(trailingOnly=TRUE)
if (length(command.args) != 7L) {
    print(command.args)
    stop (paste(collapse="\n",capture.output({
        cat(fill=TRUE,
            '
Usage: Rscript prepare-covid19-submissions.R <year> <epi.week> <crowdcast.tgz.dirpath> <submissions.root.dirpath> <do.plot> <plots.root.dirpath> <create.destination.dirs>; where:
 - year is probably 2020
 - epi.week is a two-digit epi-week-of-year number
 - crowdcast.tgz.dirpath is the path of the directory containing dfarrow epicast tarballs
 - submissions.root.dirpath is the draft submission directory or the submission repository directory
 - do.plot is a single Boolean (TRUE/FALSE) indicating whether or not to prepare plots
 - plots.root.dirpath is a directory indicating where to place plots if requested (and can be an arbitrary string otherwise)
 - create.destination.dirs is a single Boolean (TRUE/FALSE) indicating whether or not to create spreadsheet and plot destination directories --- use FALSE to double-check that pre-existing directories have the right structure.
')
    })))
}
year = as.integer(command.args[[1L]])
epi.week = as.integer(command.args[[2L]])
crowdcast.tgz.dirpath = as.character(command.args[[3L]])
submissions.root.dirpath = as.character(command.args[[4L]])
do.plot = as.logical(command.args[[5L]])
plots.root.dirpath = as.character(command.args[[6L]])
create.destination.dirs = as.logical(command.args[[7L]])


## --- Utils: ---
ensure_dirpath = function(dirpath) {
    if (length(dirpath) != 1L) {
        stop ('Only tested for one dirpath at a time.')
    }
    if (!dir.exists(dirpath)) {
        if (create.destination.dirs) {
            dir.create(dirpath)
        } else {
            stop (sprintf('Destination directory "%s" does not already exist.  Maybe one or both of {submissions.root.dirpath, spreadsheets.root.dirpath} have not been used before and needs a single-time set-up by using create.destination.dirs=TRUE, or one or both of these dirpaths were misspecified.'))
        }
    }
}



## --- Set up destinations: ---
cmu.delphi.entry.names = c("CMU_Delphi-Crowdcast","CMU_Delphi-Crowdcast_MTurk","CMU_Delphi-Stat_Nowcast")
ensure_dirpath(submissions.root.dirpath)
natreg.submissions.dirpath = file.path(submissions.root.dirpath, "nation-region-forecast-data")
state.submissions.dirpath = file.path(submissions.root.dirpath, "state-forecast-data")
for (geography.dirpath in c(natreg.submissions.dirpath, state.submissions.dirpath)) {
    ensure_dirpath(geography.dirpath)
    for (cmu.delphi.entry.name in cmu.delphi.entry.names) {
        ensure_dirpath(file.path(geography.dirpath, cmu.delphi.entry.name))
    }
}
if (do.plot) {
    ensure_dirpath(plots.root.dirpath)
    natreg.plots.dirpath = file.path(plots.root.dirpath, "nation-region-forecast-plots")
    state.plots.dirpath = file.path(plots.root.dirpath, "state-forecast-plots")
    for (geography.dirpath in c(natreg.plots.dirpath, state.plots.dirpath)) {
        ensure_dirpath(geography.dirpath)
        for (cmu.delphi.entry.name in cmu.delphi.entry.names) {
            ensure_dirpath(file.path(geography.dirpath, cmu.delphi.entry.name))
        }
    }
}



## --- Prepare Crowdcast submission: ---
## Select a tarball:
crowdcast.tgz.pattern = sprintf("^dfarrow_covid_epicast_%d%02d_([0-9]+).tgz$", year, epi.week)
candidate.crowdcast.tgz.filenames = list.files(crowdcast.tgz.dirpath) %>>%
    {.[grep(crowdcast.tgz.pattern, .)]}
selected.crowdcast.tgz.filename =
    if (length(candidate.crowdcast.tgz.filenames) == 0L) {
        stop ('Cannot find tgz for this week.')
    } else if (length(candidate.crowdcast.tgz.filenames) == 1L) {
        candidate.crowdcast.tgz.filenames[[1L]]
    } else {
        version.numbers = as.integer(sub(crowdcast.tgz.pattern, "\\1", candidate.crowdcast.tgz.filenames))
        selected.filename = candidate.crowdcast.tgz.filenames[[which.max(version.numbers)]]
        warning(sprintf('There appeared to be multiple tarballs for the given week.  Selected "%s"', selected.filename))
        selected.filename
    }
selected.crowdcast.tgz.filepath = file.path(crowdcast.tgz.dirpath, selected.crowdcast.tgz.filename)
## Extract tarball to temporary directory:
extraction.dirpath = tempdir()
untar.code = untar(selected.crowdcast.tgz.filepath, c("epicast-regional.csv","epicast-state.csv"),
                   exdir=extraction.dirpath,
                   compressed="gzip")
if (untar.code != 0L) {
    stop (sprintf('Error extracting tarball; error code: %d',untar.code))
}
cat('Extracted Crowdcast tarball.', fill=TRUE)
## Filter targets within the csv's, write to csv's in destinations:
crowdcast.natreg.spreadsheet.filepath = file.path(natreg.submissions.dirpath, "CMU_Delphi-Crowdcast", sprintf("%d-ew%02d-CMU_Delphi-Crowdcast.csv", year, epi.week))
readr::read_csv(file.path(extraction.dirpath, "epicast-regional.csv"),
                col_types="ccccc") %>>%
    dplyr::filter(target %in% c("1 wk ahead","2 wk ahead")) %>>%
    readr::write_csv(crowdcast.natreg.spreadsheet.filepath)
crowdcast.state.spreadsheet.filepath = file.path(state.submissions.dirpath, "CMU_Delphi-Crowdcast", sprintf("%d-ew%02d-CMU_Delphi-Crowdcast.csv", year, epi.week))
readr::read_csv(file.path(extraction.dirpath, "epicast-state.csv"),
                col_types="ccccc") %>>%
    dplyr::filter(target %in% c("1 wk ahead","2 wk ahead")) %>>%
    readr::write_csv(crowdcast.state.spreadsheet.filepath)
cat('Wrote out spreadsheets.', fill=TRUE)
## Prepare plots:
if (do.plot) {
    plot.covid.forecast(crowdcast.natreg.spreadsheet.filepath, file.path(natreg.plots.dirpath, "CMU_Delphi-Crowdcast"))
    plot.covid.forecast(crowdcast.state.spreadsheet.filepath, file.path(state.plots.dirpath, "CMU_Delphi-Crowdcast"))
}
cat('Wrote out plots.', fill=TRUE)
## Clean up temporary directory immediately:
file.remove(file.path(extraction.dirpath, "epicast-regional.csv"))
file.remove(file.path(extraction.dirpath, "epicast-state.csv"))
unlink(extraction.dirpath)
cat('Cleaned up temporary directory.', fill=TRUE)
