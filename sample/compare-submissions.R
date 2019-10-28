

year.week.df = seasonModelWeekToYearWeekDF(s.prospective.seasons, w.prospective.model.weeks, 31L,3L)
stopifnot(nrow(year.week.df)==1L)
week.label = sprintf("EW%02d-%d",year.week.df[["week"]],year.week.df[["year"]])

tbl1 =
    "~/files/nosync/cdc-flusight-ensemble/model-forecasts/real-time-component-models/Delphi_ExtendedDeltaDensity/%s-Delphi_ExtendedDeltaDensity.csv" %>>%
    ## "~/files/nosync/cdc-flusight-ensemble/model-forecasts/real-time-component-models/Delphi_MarkovianDeltaDensity/%s-Delphi_MarkovianDeltaDensity.csv" %>>%
    sprintf(week.label) %>>%
    readr::read_csv()

tbl2 =
    ## "~/Dropbox/DELPHI/forecasting/flu/cdc_challenge_2018/collab-submissions/%1$s/%1$s-Delphi_DeltaDensity_PackageDefaults.csv" %>>%
    ## sprintf(week.label) %>>%
    "~/Dropbox/DELPHI/forecasting/flu/cdc_challenge_2019/flusight-natreg-run/stat-spreadsheets/EW42-Delphi-Stat-2019-10-25.csv" %>>%
    readr::read_csv()

list(tbl1, tbl2) %>>%
    lapply(function(tbl) {
        tbl %>>%
            dplyr::filter(Type=="Point") %>>%
            {.}
    }) %>>%
    {dplyr::full_join(.[[1L]],.[[2L]],by=c("Location","Target","Type","Unit","Bin_start_incl","Bin_end_notincl"))} %>>%
    dplyr::mutate(absdiff = abs(Value.x-Value.y)) %>>%
    dplyr::mutate(reldiff = 2*abs(Value.x-Value.y)/(Value.x+Value.y)) %>>%
    dplyr::select(Location,Target,Bin_start_incl,Value.x,Value.y,absdiff,reldiff,dplyr::everything()) %>>%
    ## dplyr::arrange(-reldiff)
    dplyr::arrange(-absdiff)

list(tbl1, tbl2) %>>%
    lapply(function(tbl) {
        tbl %>>%
            dplyr::filter(Type=="Bin") %>>%
            {.}
    }) %>>%
    {dplyr::full_join(.[[1L]],.[[2L]],by=c("Location","Target","Type","Unit","Bin_start_incl","Bin_end_notincl"))} %>>%
    dplyr::mutate(absdiff = abs(Value.x-Value.y)) %>>%
    dplyr::mutate(reldiff = 2*abs(Value.x-Value.y)/(Value.x+Value.y)) %>>%
    dplyr::select(Location,Target,Bin_start_incl,Value.x,Value.y,absdiff,reldiff,dplyr::everything()) %>>%
    ## dplyr::arrange(-reldiff)
    dplyr::arrange(-absdiff) %>>%
    dplyr::select(Location, Target, Bin_start_incl, absdiff, Value.x, Value.y) %>>%
    {.}
