library("pipeR")

natreg.stat.spreadsheet.dir = "~/files/nosync/epiforecast-epiproject/flusight-natreg-run/stat-spreadsheets"
state.stat.spreadsheet.dir = "~/files/nosync/epiforecast-epiproject/flusight-state-run/stat-spreadsheets"
combined.stat.spreadsheet.dir = "~/files/nosync/epiforecast-epiproject/flusight-combined-run/stat-spreadsheets"
if (!dir.exists(combined.stat.spreadsheet.dir)) {
  dir.create(combined.stat.spreadsheet.dir, recursive=TRUE)
}

spreadsheet_name_to_epi_week = function(spreadsheet.name) {
  as.integer(substring(spreadsheet.name, 3L,4L))
}

spreadsheet_info_tbl = function(spreadsheet.names, tag) {
  duplicate.info.tbl =
    tibble::tibble(spreadsheet.name=spreadsheet.names) %>>%
    dplyr::mutate(epi.week=spreadsheet_name_to_epi_week(spreadsheet.name)) %>>%
    dplyr::group_by(epi.week) %>>%
    dplyr::mutate(epi.week.is.duplicated=n()>1L) %>>%
    dplyr::ungroup()
  duplicated.epi.weeks = duplicate.info.tbl %>>%
    dplyr::filter(epi.week.is.duplicated) %>>%
    dplyr::distinct(epi.week) %>>%
    (epi.week)
  print(paste0("Don't know what to do in case of multiple spreadsheets for the same EW in the same project ",tag,"; skipping spreadsheets for the following EW's (if any): ", paste(duplicated.epi.weeks, collapse=", ")))
  return (
    duplicate.info.tbl %>>%
    dplyr::filter(!epi.week.is.duplicated) %>>%
    dplyr::select(-epi.week.is.duplicated) %>>%
    stats::setNames(c(paste0(tag,".spreadsheet.name"), names(.)[-1L])) %>>%
    {.}
  )
}

natreg.stat.spreadsheet.info.tbl =
  spreadsheet_info_tbl(list.files(natreg.stat.spreadsheet.dir), "natreg.stat")
state.stat.spreadsheet.info.tbl =
  spreadsheet_info_tbl(list.files(state.stat.spreadsheet.dir), "state.stat")

print("Skipping EW in natreg but not state (if any):")
print(setdiff(natreg.stat.spreadsheet.info.tbl[["epi.week"]],
              state.stat.spreadsheet.info.tbl[["epi.week"]]))
print("Skipping EW in state but not natreg (if any):")
print(setdiff(state.stat.spreadsheet.info.tbl[["epi.week"]],
              natreg.stat.spreadsheet.info.tbl[["epi.week"]]))

combined.stat.spreadsheet.info.tbl =
  dplyr::inner_join(natreg.stat.spreadsheet.info.tbl, state.stat.spreadsheet.info.tbl, by="epi.week")

combined.stat.spreadsheet.info.tbl %>>%
  dplyr::rowwise() %>>%
  dplyr::mutate(success={
    tryCatch({
      natreg.spreadsheet = readr::read_csv(file.path(natreg.stat.spreadsheet.dir,
                                                     natreg.stat.spreadsheet.name),
                                           col_types=readr::cols())
      state.spreadsheet = readr::read_csv(file.path(state.stat.spreadsheet.dir,
                                                    state.stat.spreadsheet.name),
                                          col_types=attr(natreg.spreadsheet,"spec"))
      combined.spreadsheet = dplyr::bind_rows(natreg.spreadsheet, state.spreadsheet)
      natreg.date = as.Date(stringr::str_sub(natreg.stat.spreadsheet.name, -14L, -5L))
      state.date = as.Date(stringr::str_sub(natreg.stat.spreadsheet.name, -14L, -5L))
      combined.date = max(natreg.date, state.date)
      combined.stat.spreadsheet.name =
        natreg.stat.spreadsheet.name %>>%
        stringr::str_replace_all("Delphi-Stat", "Delphi-Stat-combined") %>>%
        stringr::str_replace_all(as.character(natreg.date), as.character(combined.date)) %>>%
        {.}
      readr::write_csv(combined.spreadsheet,
                       file.path(combined.stat.spreadsheet.dir, combined.stat.spreadsheet.name))
      TRUE
    },
    error=function(e) {
      FALSE
    })
  }) %>>%
  dplyr::select(success, dplyr::everything())
