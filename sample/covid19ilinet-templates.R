
covid19ilinet.natreg.spreadsheet.template =
    fetchUpdatingResource(
        function() {
            read.csv(textConnection(RCurl::getURL("https://raw.githubusercontent.com/cdcepi/COVID-19-ILI-forecasting/master/templates-and-data/covid19-ili-forecast-national-regional-template.csv")), check.names=FALSE, stringsAsFactors=FALSE)
        },
        function(fetch.response) {
            return ()
        },
        cache.file.prefix=file.path(epidata.cache.dir,"covid19ilinet_natreg_spreadsheet_template"),
        cache.invalidation.period=as.difftime(8L, units="hours"),
        force.cache.invalidation=FALSE
    ) %>>%
    tibble::as_tibble()

covid19ilinet.state.spreadsheet.template =
    fetchUpdatingResource(
        function() {
            read.csv(textConnection(RCurl::getURL("https://raw.githubusercontent.com/cdcepi/COVID-19-ILI-forecasting/master/templates-and-data/covid19-ili-forecast-state-template.csv")), check.names=FALSE, stringsAsFactors=FALSE)
        },
        function(fetch.response) {
            return ()
        },
        cache.file.prefix=file.path(epidata.cache.dir,"covid19ilinet_state_spreadsheet_template"),
        cache.invalidation.period=as.difftime(8L, units="hours"),
        force.cache.invalidation=FALSE
    ) %>>%
    tibble::as_tibble()
