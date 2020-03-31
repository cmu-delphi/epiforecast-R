
## xxx Approximate conversion.  Not appropriate if significant mass is above 13.0 threshold.

devtools::load_all("../epiforecast")

epidata.cache.dir = "~/.epiforecast-cache"

covid19ilinet.natreg.spreadsheet.template =
    fetchUpdatingResource(
        function() {
            read.csv(textConnection(RCurl::getURL("https://raw.githubusercontent.com/cdcepi/COVID-19-ILI-forecasting/master/templates-and-data/covid19-ili-forecast-national-regional-template.csv")), check.names=FALSE, stringsAsFactors=FALSE)
        },
        function(fetch.response) {
            return ()
        },
        cache.file.prefix=file.path(epidata.cache.dir,"covid19ilinet_natreg_spreadsheet_template"),
        cache.invalidation.period=as.difftime(2L, units="days"),
        force.cache.invalidation=TRUE
    )
covid19ilinet.state.spreadsheet.template =
    fetchUpdatingResource(
        function() {
            read.csv(textConnection(RCurl::getURL("https://raw.githubusercontent.com/cdcepi/COVID-19-ILI-forecasting/master/templates-and-data/covid19-ili-forecast-state-template.csv")), check.names=FALSE, stringsAsFactors=FALSE)
        },
        function(fetch.response) {
            return ()
        },
        cache.file.prefix=file.path(epidata.cache.dir,"covid19ilinet_state_spreadsheet_template"),
        cache.invalidation.period=as.difftime(2L, units="days"),
        force.cache.invalidation=TRUE
    )

overlap_portions = function(start, end, of.start, of.end) {
    overlap.start = pmax(start, of.start)
    overlap.end = pmin(end, of.end)
    overlap.start <- pmin(overlap.start, overlap.end)
    (overlap.end-overlap.start)/(of.end-of.start)
}

covid19ilinet.effective.breaks =
    covid19ilinet.percentage.bin.info[["breaks"]]
flusight2016.effective.breaks =
    flusight2016.percentage.bin.info[["breaks"]] %>>%
    {.[.!=0 & .!=100] <- .[.!=0 & .!=100]-0.05; .} %>>%
    {.}

covid19ilinet.percentage.bin.weight.from.flusight2016.mat =
    Matrix::Matrix(sparse=TRUE,
                   outer(seq_len(length(covid19ilinet.effective.breaks)-1L),
                         seq_len(length(flusight2016.effective.breaks)-1L),
                         function(is, js) {
                             overlap_portions(
                                 covid19ilinet.effective.breaks[is],
                                 covid19ilinet.effective.breaks[is+1L],
                                 flusight2016.effective.breaks[js],
                                 flusight2016.effective.breaks[js+1L]
                             )
                         })
                   ) %>>%
    `rownames<-`(levels(cut(0, covid19ilinet.effective.breaks, include.lowest=covid19ilinet.percentage.bin.info$rightmost.closed, right=FALSE))) %>>%
    `colnames<-`(levels(cut(0, flusight2016.effective.breaks, include.lowest=flusight2016.percentage.bin.info$rightmost.closed, right=FALSE))) %>>%
    {.}

covid19ilinet.percentage.bin.weight.from.flusight2016.mat %>>%
    (mat ~ {
        summary(mat) %>>%
            dplyr::mutate(rowname=rownames(mat)[i],
                          colname=colnames(mat)[j])
    })

stopifnot (all(abs(Matrix::colSums(covid19ilinet.percentage.bin.weight.from.flusight2016.mat)-1)<1e-8))

covid19ilinet_nearcast_spreadsheet_from_flusight = function(flusight.ilinet.spreadsheet, covid19ilinet.spreadsheet.template) {
    targets.to.translate = c("1 wk ahead","2 wk ahead")
    mass.on.13ups =
        flusight.ilinet.spreadsheet %>>%
        dplyr::filter(Target %in% targets.to.translate,
                      Type == "Bin",
                      Bin_start_incl=="13") %>>%
        (~ {
            if (nrow(.)==0L) {
                stop ('Check is checking nothing.')
            }
        }) %>>%
        {.}
    max.mass.on.13up = max(mass.on.13ups[["Value"]])
    if (max.mass.on.13up > 0.05) {
        print(dplyr::arrange(mass.on.13ups, -Value))
        print(max.mass.on.13up)
        warning('Greater than 5% of mass was placed in a [13,100] Bin for a target that is supposed to be translated over to the new bin system.')
    }
    covid19ilinet.spreadsheet =
        flusight.ilinet.spreadsheet %>>%
        dplyr::filter(Target %in% targets.to.translate) %>>%
        dplyr::group_by(Location, Target, Type, Unit) %>>%
        dplyr::do(
                   tibble::tibble(
                               location = .[["Location"]][[1L]],
                               target = .[["Target"]][[1L]],
                               type = tolower(.[["Type"]][[1L]]),
                               bin =
                                   if (.[["Type"]][[1L]] == "Point") {
                                       NA_character_
                                   } else if (.[["Type"]][[1L]] == "Bin") {
                                       sprintf("%.1f", covid19ilinet.percentage.bin.info[["breaks"]] %>>% {.[-length(.)]})
                                   } else stop('Unexpected Type encountered.'),
                               value =
                                   if (.[["Type"]][[1L]] == "Point") {
                                       .[["Value"]]
                                   } else if (.[["Type"]][[1L]] == "Bin"){
                                       as.vector(covid19ilinet.percentage.bin.weight.from.flusight2016.mat %*% .[["Value"]])
                                   } else stop('Unexpected Type encountered.')
                           )
               ) %>>%
        dplyr::ungroup() %>>%
        dplyr::select(-Location, -Target, -Type, -Unit) %>>%
        {.}
    bad.rows = dplyr::anti_join(covid19ilinet.spreadsheet, covid19ilinet.spreadsheet.template, by=c("location","target", "type", "bin"))
    if (nrow(bad.rows)!=0L) {
        print(bad.rows)
        stop ('Formed rows with index column entries that are not present in the template.')
        ## fixme todo check for duplicates
    }
    ## It's okay to omit rows the other way around when submitting forecasts for only a subset of targets.
    covid19ilinet.spreadsheet %>>%
        dplyr::mutate(
                   location=ordered(location, unique(covid19ilinet.spreadsheet.template[["location"]])),
                   target=ordered(target, unique(covid19ilinet.spreadsheet.template[["target"]])),
                   type=ordered(type, unique(covid19ilinet.spreadsheet.template[["type"]]))
               ) %>>%
        dplyr::arrange(location, target, type) %>>%
        {.}
}

ew = 11L

readr::read_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/cdc_challenge_2019/flusight-natreg-run/stat-spreadsheets/EW%02d-Delphi-Stat-2020-03-23.csv",ew), col_types="ccccccd") %>>%
    covid19ilinet_nearcast_spreadsheet_from_flusight(covid19ilinet.natreg.spreadsheet.template) %>>%
    readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/nation-region-forecast-data/CMU_Delphi-Stat_Nowcast/2020-ew%02d-CMU_Delphi-Stat_Nowcast.csv",ew))

readr::read_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/cdc_challenge_2019/flusight-state-run/stat-spreadsheets/EW%02d-Delphi-Stat-StateILI-2020-03-23.csv",ew), col_types="ccccccd") %>>%
    covid19ilinet_nearcast_spreadsheet_from_flusight(covid19ilinet.state.spreadsheet.template) %>>%
    readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/state-forecast-data/CMU_Delphi-Stat_Nowcast/2020-ew%02d-CMU_Delphi-Stat_Nowcast.csv",ew))



readr::read_csv(sprintf("~/Dropbox/private/EW%02d-delphi-epicast-regional-2020-03-23.csv",ew), col_types="ccccccd") %>>%
    covid19ilinet_nearcast_spreadsheet_from_flusight(covid19ilinet.natreg.spreadsheet.template) %>>%
    readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/nation-region-forecast-data/CMU_Delphi-Crowdcast/2020-ew%02d-CMU_Delphi-Crowdcast.csv",ew))

readr::read_csv(sprintf("~/Dropbox/private/EW%02d-delphi-epicast-state-2020-03-23.csv",ew), col_types="ccccccd") %>>%
    covid19ilinet_nearcast_spreadsheet_from_flusight(covid19ilinet.state.spreadsheet.template) %>>%
    readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/state-forecast-data/CMU_Delphi-Crowdcast/2020-ew%02d-CMU_Delphi-Crowdcast.csv",ew))



## readr::read_csv(sprintf("~/Dropbox/private/EW%02d-delphi-epicast-mturk-regional-2020-03-23.csv",ew), col_types="ccccccd") %>>%
##     covid19ilinet_nearcast_spreadsheet_from_flusight(covid19ilinet.natreg.spreadsheet.template) %>>%
##     readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/nation-region-forecast-data/CMU_Delphi-Crowdcast_MTurk/2020-ew%02d-CMU_Delphi-Crowdcast_MTurk.csv",ew))

## readr::read_csv(sprintf("~/Dropbox/private/EW%02d-delphi-epicast-mturk-states-2020-03-23.csv",ew), col_types="ccccccd") %>>%
##     covid19ilinet_nearcast_spreadsheet_from_flusight(covid19ilinet.state.spreadsheet.template) %>>%
##     readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/state-forecast-data/CMU_Delphi-Crowdcast_MTurk/2020-ew%02d-CMU_Delphi-Crowdcast_MTurk.csv",ew))
