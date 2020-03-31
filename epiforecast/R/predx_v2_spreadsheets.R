
##' Reformat a predx v1 spreadsheet to a predx v2 spreadsheet
##'
##' @param old.spreadsheet predx v1 spreadsheet
##' @param old.spreadsheet predx v2 spreadsheet template
##' @return predx v2 spreadsheet
##'
##' @export
reformat_to_predx_v2_spreadsheet = function(old.spreadsheet, new.template) {
    new.spreadsheet.minus.ordering = old.spreadsheet %>>%
        dplyr::group_by(Location, Target, Type, Unit) %>>%
        dplyr::do(
                   tibble::tibble(
                               location = .[["Location"]][[1L]],
                               target = .[["Target"]][[1L]],
                               type = tolower(.[["Type"]][[1L]]),
                               bin =
                                   if (.[["Type"]][[1L]] == "point") {
                                       NA_character_
                                   } else if (.[["Type"]][[1L]] == "bin") {
                                       sprintf("%.1f", covid19ilinet.percentage.bin.info[["breaks"]] %>>% {.[-length(.)]})
                                   } else stop(sprintf('Unexpected Type encountered: %s.', .[["Type"]][[1L]])),
                               value = .[["Value"]]
                           )
               ) %>>%
        dplyr::ungroup() %>>%
        dplyr::select(-Location, -Target, -Type, -Unit) %>>%
        {.}
    bad.rows = dplyr::anti_join(new.spreadsheet.minus.ordering, new.template, by=c("location","target", "type", "bin"))
    if (nrow(bad.rows)!=0L) {
        stop (paste(collapse="\n", capture.output({
            cat('Formed rows with index column entries that are not present in the template:', fill=getOption('width')-nchar('Error: '))
            print(bad.rows)
        })))
        ## fixme todo check for duplicates
    }
    overlarge.groups = new.spreadsheet.minus.ordering %>>%
        dplyr::group_by(location, target, type, bin) %>>%
        dplyr::summarize(count=dplyr::n()) %>>%
        dplyr::ungroup() %>>%
        dplyr::filter(count != 1L) %>>%
        {.}
    if (nrow(overlarge.groups) != 0L) {
        stop (paste(collapse="\n", capture.output({
            cat('Formed overlarge groups:', fill=getOption('width')-nchar('Error: '))
            print(overlarge.groups)
        })))
    }
    ## It's okay to omit rows the other way around when submitting forecasts for only a subset of targets.
    new.spreadsheet = new.spreadsheet.minus.ordering %>>%
        dplyr::mutate(
                   location=ordered(location, unique(new.template[["location"]])),
                   target=ordered(target, unique(new.template[["target"]])),
                   type=ordered(type, unique(new.template[["type"]]))
               ) %>>%
        dplyr::arrange(location, target, type) %>>%
        {.}
    return (new.spreadsheet)
}
