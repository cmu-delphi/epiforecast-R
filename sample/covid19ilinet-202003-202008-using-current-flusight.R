
library("pipeR")

devtools::load_all("../epiforecast")

config_env = function(filename) {
    source(filename, chdir=TRUE, local=TRUE)
    environment()
}

flusight.natreg.config = config_env("natreg-config.R")
flusight.high.state.config = config_env("high-state-config.R")
flusight.low.state.config = config_env("low-state-config.R")

covid19ilinet.202003.202008.natreg.config = config_env("covid19ilinet-202003-202008-natreg-config.R")

derived_config_env = function(config, covid19config) {
    result = new.env(parent=config)
    ## todo voxel data target settings...
    ## todo observed/valuation trajectories...
    obj.names = c(
        "last.extended.epi.week",
        "epidata_df_to_chopped_trajectory_df",
        "target_trajectory_preprocessor",
        "full_dat_fixup"
    )
    for (obj.name in obj.names) {
        result[[obj.name]] <- covid19config[[obj.name]]
    }
    ## todo enable?:
    ## result[["b.backcasters"]][c(
    ##           "ignorant",
    ##           "quantile_arx_backnowcast",
    ##           "quantile_arx_pancast"
    ##       )] <-
    ##     covid19config[["b.backcasters"]][c(
    ##                      "RevisionIgnorant",
    ##                      "QARXBacknowcastAllShifts",
    ##                      "QARXPancastAllShifts"
    ##                  )]
    ## todo f.forecasters ???
    result[["f.forecasters"]] <- get0("f.forecasters", config) %>>%
        {
            ## Having issues with `[[<-` with functions; use `[<-` instead,
            ## which has dimnames names dropping issue, and fix up:
            .[grep("(Extended)|(_Delta)",namesp(.))] <- list(
                function(full.dat, baseline=0, max.n.sims=1000L) {
                    twkde.sim(full.dat, baseline=baseline, max.n.sims=max.n.sims
                            , max.shifts=c(rep(10L,10L),10:1,rep(0L,3L),1:10,rep(10L,20L),rep(10L,covid19ilinet.202003.202008.natreg.config[["last.extended.epi.week"]]-usa.flu.first.week.of.season+1L))
                              )
                }
            )
            . %>>% with_dimnamesnames("Forecaster")
        }
    result[["t.target.specs"]] <- get0("t.target.specs", config) %>>%
        {
            .[c("1 wk ahead","2 wk ahead")] <-
                covid19config[["t.target.specs"]][c("1 wk ahead","2 wk ahead")];
            ## This set-by-name operation destroys the dimnamesnames; restore:
            .[c("1 wk ahead","2 wk ahead"),drop=FALSE] %>>% with_dimnamesnames("Target")
        }
    result[["m.forecast.types"]] <- get0("m.forecast.types", config) %>>%
        {
            .[c("Point","Bin")] <-
                covid19config[["m.forecast.types"]][c("point","bin")]
            ## This set-by-name operation destroys the dimnamesnames; restore:
            . %>>% with_dimnamesnames("Type")
        }
    result[["e.ensemble.partial.weighting.scheme.wgt.indexer.lists"]] <-
        get0("e.ensemble.partial.weighting.scheme.wgt.indexer.lists", config)["target-based",drop=FALSE]
    result[["epiproject.cache.dir"]] <- file.path(covid19config[["epiproject.cache.dir"]], basename(config[["epiproject.cache.dir"]]))
    return (result)
}

flusight.configs = list(
    natreg=flusight.natreg.config,
    low.state=flusight.low.state.config,
    high.state=flusight.high.state.config
)

derived.configs = lapply(flusight.configs, derived_config_env, covid19ilinet.202003.202008.natreg.config)

selected.ensemble.weighting.scheme = "target-based"

## xxx Copy in flusight-weightsets beforehand.
for (derived.config.i in seq_along(derived.configs)) {
    derived.config = derived.configs[[derived.config.i]]
    derived.config.name = names(derived.configs)[[derived.config.i]]
    weightsets.dirpath = file.path(
        derived.config[["epiproject.cache.dir"]],
        "..", "flusight-weightsets",
        sprintf(
            "e.prospective.ensemble.weightsets.%s",
            gsub("\\.","-",derived.config.name)
        )
    )
    destpath = file.path(
        derived.config[["epiproject.cache.dir"]],
        "e.prospective.ensemble.weightsets/"
    )
    if (!file.exists(destpath)) {
        R.utils::copyDirectory(weightsets.dirpath, destpath, recursive=TRUE)
    }
}

for (derived.config.i in seq_along(derived.configs)) {
    derived.config = derived.configs[[derived.config.i]]
    (function() {
        source("gen-prospective-component-forecasts.R", local=TRUE)
        swgtmbf.retro.component.forecast.values = no_join(NULL)
        source("gen-prospective-ensemble-forecasts.R", local=TRUE)
        save_spreadsheets(
            swge.prospective.ensemble.target.multicasts[,,,selected.ensemble.weighting.scheme,drop=FALSE],
            swg.prospective.voxel.data,
            t.target.specs, m.forecast.types,
            epigroup.colname,
            file.path(epiproject.cache.dir,"stat-spreadsheets")
        )
    }) %>>%
        `environment<-`(derived.config) %>>%
        {.()}
}

epidata.cache.dir = covid19ilinet.202003.202008.natreg.config[["epidata.cache.dir"]]

source("covid19ilinet-templates.R")

dir.create("~/Dropbox/private/covid19natreg/")
dir.create("~/Dropbox/private/covid19state/")



## readr::read_csv(
##            file.path(derived.configs[["natreg"]][["epiproject.cache.dir"]],"stat-spreadsheets","2019-2020.MW63.target-based.csv")
##          , col_types="ccccccd"
##        ) %>>%
##     tibble::as_tibble() %>>%
##     reformat_to_predx_v2_spreadsheet(covid19ilinet.natreg.spreadsheet.template) %>>%
##     (? .) %>>%
##     (? nrow(.)) %>>%
##     ## readr::write_csv("~/Dropbox/private/covid19natreg/2020-ew11-CMU_Delphi-Stat_Nowcast.csv")
##     readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/nation-region-forecast-data/CMU_Delphi-Stat_Nowcast/2020-ew%02d-CMU_Delphi-Stat_Nowcast.csv",11L))



## dplyr::bind_rows(
##            readr::read_csv(file.path(derived.configs[["low.state"]][["epiproject.cache.dir"]],"stat-spreadsheets","2019-2020.MW63.target-based.csv")
##                          , col_types="ccccccd"
##                            ),
##            readr::read_csv(
##                       file.path(derived.configs[["high.state"]][["epiproject.cache.dir"]],"stat-spreadsheets","2019-2020.MW63.target-based.csv")
##                     , col_types="ccccccd"
##                   )
##        ) %>>%
##     reformat_to_predx_v2_spreadsheet(covid19ilinet.state.spreadsheet.template) %>>%
##     (? .) %>>%
##     (? nrow(.)) %>>%
##     ## readr::write_csv("~/Dropbox/private/covid19state/2020-ew11-CMU_Delphi-Stat_Nowcast.csv")
##     readr::write_csv(sprintf("~/Dropbox/DELPHI/forecasting/flu/covid19ilinet_202003_202008/draft-submissions/state-forecast-data/CMU_Delphi-Stat_Nowcast/2020-ew%02d-CMU_Delphi-Stat_Nowcast.csv",11L))



## todo check for Crowdcast sensor availability
