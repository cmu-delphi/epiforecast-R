
source("gen-retro-component-forecasts.R")
source("gen-retro-ensemble-forecasts.R")



## apply(swgtmbf.retro.component.evaluations[(2003:2009)%>>%paste0("/",.+1L),,,,,,,drop=FALSE],
##       c(5L,7L), mean, na.rm=TRUE)
## apply(swgtme.retro.ensemble.evaluations[(2003:2009)%>>%paste0("/",.+1L),,,,,,drop=FALSE],
##       5:6, mean, na.rm=TRUE)

## apply(swgtmbf.retro.component.evaluations, c(7L,4:5), mean, na.rm=TRUE)
## apply(swgtme.retro.ensemble.evaluations, c(6L,4:5), mean, na.rm=TRUE)

## apply(swgtme.retro.ensemble.evaluations, c(1,4:5), mean, na.rm=TRUE)

## todo another level of stacking?  ensemble of ensembles?

## swbf.retro.forecast.spreadsheets = map_join(
##   target_multicast_epigroup_forecast_table,
##   swgbf.retro.component.target.multicasts,
##   swg.retro.voxel.data,
##   no_join(t.target.specs), no_join(m.forecast.types),
##   epigroup.colname
## ) %>>%
##   apply(c(1:2,4:5), dplyr::bind_rows)

## spreadsheet.template = readr::read_csv("~/long_flu_submission_template_1718.csv", col_types=readr::cols())
## attr(spreadsheet.template, "spec") <- NULL
## spreadsheet.to.check = swbf.retro.forecast.spreadsheets[["2010/2011",1L,1L,1L]]

## spreadsheet.to.check %>>% dplyr::mutate(A="a") %>>% dplyr::full_join(qwer %>>% dplyr::mutate(B="b"), by=c("Location","Target","Type","Unit","Bin_start_incl","Bin_end_notincl")) %>>% dplyr::filter(is.na(A) | is.na(B))
## class(spreadsheet.to.check[["Value"]])==class(spreadsheet.template[["Value"]])
## identical(sapply(spreadsheet.to.check, class), sapply(spreadsheet.template, class))
## identical(sapply(spreadsheet.to.check, typeof), sapply(spreadsheet.template, typeof))
## identical(sapply(spreadsheet.to.check, mode), sapply(spreadsheet.template, mode))
## identical(sapply(spreadsheet.to.check, storage.mode), sapply(spreadsheet.template, storage.mode))
## Map(all.equal,
##     spreadsheet.to.check %>>% dplyr::select(-Value),
##     spreadsheet.template %>>% dplyr::select(-Value))
## Map(identical,
##     spreadsheet.to.check %>>% dplyr::select(-Value),
##     spreadsheet.template %>>% dplyr::select(-Value))
## identical(class(spreadsheet.to.check), class(spreadsheet.template))
## identical(typeof(spreadsheet.to.check), typeof(spreadsheet.template))
## identical(mode(spreadsheet.to.check), mode(spreadsheet.template))
## identical(storage.mode(spreadsheet.to.check), storage.mode(spreadsheet.template))
## ## all.equal overload does not check attributes; if "spec" attribute set by
## ## readr is included, the all.equal check passes but the identical check fails;
## ## removing this attribute makes the identical check pass
## all.equal(spreadsheet.to.check %>>% dplyr::select(-Value),
##           spreadsheet.template %>>% dplyr::select(-Value))
## identical(spreadsheet.to.check %>>% dplyr::select(-Value),
##           spreadsheet.template %>>% dplyr::select(-Value))

source("gen-prospective-component-forecasts.R")
source("gen-prospective-ensemble-forecasts.R")

## todo automatic appropriate parallel_dim_i setting
## fixme finish documenting map_join
## todo use simplified sim's to strip out unnecessary elements
## todo sw.epidata --- if will be borrowing across epigroups...
## todo add prospective forecasts...
## todo allow for subsets in ensemble... need to use cv_apply to get the forecasts again
## xxx target.settings -> voxel.settings?
## todo also need to strip out information from target forecasts if remains after using simplified sim's
## todo dots to lists in target spec functions, maybe elsewhere as well
## todo could also have save season-evaluation epidata to hopefully speed up calculating observed values
## todo reshape so that can include labels without exploding in size
## todo don't store sim's
## reverse order on the indexing?
## fixme no uniform fallback for Point; maybe make fallback the empirical distribution since it is smoothed now and has uniform pseudoweight?  or just remove since everything has uniform pseudoweight?
## fixme no Bin smoothing by default? warn about plotting target forecasts?
## todo joint epidata? voxel settings?
## joint backcast? joint forecast?
## fixme Virgin Islands --- little history
## todo document no_join
## fixme check on Season onset Point predictions of "none" vs. NA vs. ...
## fixme EB weight sum is too large
## fixme smooth sim targets with a Laplace kernel? or a spike + slab type --- can inflate bw to make up for mass on spike?
## fixme should adjust dimension ordering... for col major feel
## todo deal with parallel-related memory issues --- duplicate swg.voxel.data...
## todo work on speed getting voxel data --- possible to avoid storage?
## xxx consider just basing everything on filesystem contracts... no need to hold everything in memory
## fixme try to solve memory issues with mclapply env's? require interaction with disk?
## todo weighted cv_apply smearing schemes (boxcar kernel -> other kernels)
## todo tiny subset run to do some testing and development on
## fixme table verification
## fixme empirical distribution is using target's baseline with empirical curves again instead of pairing baselines and curves
## todo better fill-in in 2010 for states
## todo remove Season onset target for states
## fixme modularize into functions, document, ...
