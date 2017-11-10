## author_header begin
## Copyright (C) 2017 Logan C. Brooks
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

##' @include cv_apply.R
##' @include map_join.R
##' @include namesp.R
NULL

get_backcast = function(voxel.data, signal.name, backcaster) {
  set.seed(42L)
  backcaster(voxel.data, signal.name)
}

get_forecast = function(voxel.data, full.dat, forecaster) {
  set.seed(42L)
  forecaster(full.dat, baseline=voxel.data[["baseline"]])
}

target_forecast2 = function(voxel.data, target_trajectory_preprocessor, target.spec, simlike) {
  target.forecast = do.call(
    target_forecast,
    c(list(
      simlike,
      target.name=target.spec[["Target"]],
      target.fun=target.spec[["for_processed_trajectory"]],
      target_trajectory_preprocessor=target_trajectory_preprocessor,
      target.spec=target.spec,
      target_value_formatter=target.spec[["to_point"]],
      compute.estimates=FALSE
    ), voxel.data[["target.settings"]])
  )
  ## for use with forecast_value2; remove target.settings from this object and
  ## rely on voxel.data
  target.forecast[["target.settings"]] <- NULL
  return (target.forecast)
}

forecast_value2 = function(voxel.data, target.spec, forecast.type, target.forecast, label.bins=FALSE) {
  if ("target.settings" %in% names(target.forecast)) {
    stop ("Expected target.settings to be missing in target.forecast (as it will be filled in using voxel.data).")
  } else {
    ## todo write c and assignment overloads for target_forecast
    target.forecast[["target.settings"]] <- voxel.data[["target.settings"]]
  }
  forecast_value(target.spec, forecast.type, target.forecast, label.bins=label.bins)
}

target_multicast = function(voxel.data, full.dat, forecaster, target_trajectory_preprocessor, target.specs, forecast.types, method.settings.overrides=list()) {
  simlike = forecaster(full.dat, baseline=voxel.data[["baseline"]])
  target.forecasts = map_join(
    target_forecast2,
    no_join(voxel.data),
    target_trajectory_preprocessor, target.specs,
    no_join(simlike),
    lapply_variant=lapply,
    progress.output=FALSE
  )
  target.forecasts <- map_join(
    function(target.forecast) {
      if ("method.settings" %in% names(target.forecast)) {
        target.forecast[["method.settings"]][names(method.settings.overrides)] <- method.settings.overrides
      } else {
        target.forecast[["method.settings"]] <- method.settings.overrides
      }
      target.forecast
    },
    target.forecasts,
    lapply_variant=lapply
  )
  forecast.values = map_join(
    forecast_value2,
    no_join(voxel.data),
    target.specs, forecast.types, target.forecasts,
    lapply_variant=lapply,
    progress.output=FALSE
  )
  structure(
    c(#voxel.data[c("season","model.week","epigroup","issue")],
      list(#baseline=voxel.data[["baseline"]],
           simplified.simlike=simplified_simlike(simlike),
           forecast.values=forecast.values
           )
      ),
    class = "target_multicast"
  )
}

##' @method plot target_multicast
##' @export
##' @export plot.target_multicast
plot.target_multicast = function(x, voxel.data, target.specs, forecast.values, ...) {
  target.multicast = x
  stop ("fixme port over fanplus plots")
}

EW_NA_to_MWplus1 = function(EW_NA, n.weeks.in.season) {
  MW_NA = dplyr::if_else(EW_NA < 40L, EW_NA+n.weeks.in.season, EW_NA)
  plus1 = n.weeks.in.season+21L
  storage.mode(plus1) <- storage.mode(MW_NA)
  dplyr::if_else(is.na(MW_NA), plus1, MW_NA)
}
## fixme don't hardcode week numbers

EWnone_to_MWplus1 = function(EWnone, n.weeks.in.season) {
  EW_NA = as.integer(dplyr::recode(EWnone, none=NA_character_))
  EW_NA_to_MWplus1(EW_NA, n.weeks.in.season)
}

target_multicast_epigroup_forecast_table = function(target.multicast, voxel.data, t.target.specs, m.forecast.types) {
  epigroup.forecast.table = map_join(
    function(forecast.type, target.spec, voxel.data, forecast.value) {
      subspreadsheet =
        do.call(subspreadsheet_from_forecast_value,
                c(list(forecast.value, target.spec, forecast.type),
                  voxel.data[["target.settings"]]))
      subspreadsheet.previous.colnames = colnames(subspreadsheet)
      subspreadsheet[[epigroup.colname]] <- voxel.data[["epigroup"]]
      subspreadsheet <- subspreadsheet[c(epigroup.colname,subspreadsheet.previous.colnames)]
      subspreadsheet
    },
    m.forecast.types, t.target.specs, no_join(voxel.data),
    target.multicast[["forecast.values"]],
    lapply_variant=lapply,
    progress.output=FALSE
  ) %>>%
    dplyr::bind_rows()
  return (epigroup.forecast.table)
}

multi_spreadsheet_linlog_plot = function(multi.spreadsheet, binlabel_to_x, point_to_x, group.colname, weight.colname=NULL) {
  multi.spreadsheet %>>%
    dplyr::filter(Type == "Bin") %>>%
    dplyr::mutate_at(dplyr::vars(Bin_start_incl, Bin_end_notincl), binlabel_to_x) %>>%
    dplyr::mutate(Scale = "Linear") %>>%
    dplyr::bind_rows(dplyr::mutate(., Value=log(Value), Scale="Log")) %>>%
    ggplot2::ggplot(ggplot2::aes_string(
                               x="Bin_start_incl", y="Value",
                               size=weight.colname,
                               colour=group.colname, fill=group.colname, group=group.colname)) +
    ggplot2::facet_grid(Scale ~ ., scales="free_y") +
    ggplot2::geom_line() +
    ggplot2::geom_point(
               ggplot2::aes(x=Value, y=0),
               multi.spreadsheet %>>%
               dplyr::filter(Type == "Point") %>>%
               dplyr::mutate(Value = point_to_x(Value)) %>>%
               {.}
             ) +
    ggplot2::scale_size(range=c(1,3)) +
    ggplot2::theme(legend.position=c(1,1), legend.justification=c(1,1))
}

target_multicast_week_plot = function(target.multicast, voxel.data, t.target.specs, m.forecast.types) {
  n.weeks.in.season = length(voxel.data[["target.settings"]][["is.inseason"]])
  ##
  epigroup.forecast.table = target_multicast_epigroup_forecast_table(target.multicast, voxel.data, t.target.specs, m.forecast.types)
  ##
  epigroup.forecast.table %>>%
    dplyr::filter(Unit == "week") %>>%
    multi_spreadsheet_linlog_plot(
      function(binlabel) {
        binlabel %>>% EWnone_to_MWplus1(n.weeks.in.season)
      },
      function(point) {
        point %>>% EW_NA_to_MWplus1(n.weeks.in.season)
      },
      "Target")
}

target_multicast_percent_plot = function(target.multicast, voxel.data, t.target.specs, m.forecast.types) {
  epigroup.forecast.table = target_multicast_epigroup_forecast_table(target.multicast, voxel.data, t.target.specs, m.forecast.types)
  ##
  epigroup.forecast.table %>>%
    dplyr::filter(Unit == "percent") %>>%
    multi_spreadsheet_linlog_plot(as.numeric,
                                 identity,
                                 "Target")
}

target_multicast_combined_target_plot = function(target.multicasts, voxel.data, t.target.specs, m.forecast.types, plot.Target, Mtm.weights=NULL) {
  n.weeks.in.season = length(voxel.data[["target.settings"]][["is.inseason"]])
  Unit = t.target.specs[[plot.Target]][["unit"]][["Unit"]]
  if (Unit == "week") {
    binlabel_to_x = function(binlabel) {
      binlabel %>>% EWnone_to_MWplus1(n.weeks.in.season)
    }
    point_to_x = function(point) {
      point %>>% EW_NA_to_MWplus1(n.weeks.in.season)
    }
  } else if (Unit == "percent") {
    binlabel_to_x = as.numeric
    point_to_x = identity
  } else {
    stop ("Missing or unrecognized Unit.")
  }
  target.spreadsheet =
    lapply(target.multicasts,
           target_multicast_epigroup_forecast_table,
           voxel.data, t.target.specs, m.forecast.types) %>>%
    dplyr::bind_rows(.id="Model") %>>%
    {
      if (is.null(Mtm.weights)) {
        .
      } else {
        . %>>%
          dplyr::left_join(
                   reshape2::melt(Mtm.weights, value.name="weight") %>>%
                   dplyr::mutate_if(is.factor, as.character)
                 , by=c("Model","Target","Type")) %>>%
          dplyr::bind_rows(
                   . %>>%
                   dplyr::group_by_(.dots=colnames(.)[!colnames(.) %in% c("Model","Value","weight")]) %>>%
                   dplyr::summarize(Value=weighted.mean(Value, weight)) %>>%
                   dplyr::ungroup() %>>%
                   dplyr::mutate(Model="Selected Ensemble", weight=1)
                 )
      }
    } %>>%
    dplyr::filter(Target == plot.Target)
  multi_spreadsheet_linlog_plot(target.spreadsheet, binlabel_to_x, point_to_x, "Model",
                                if (is.null(Mtm.weights)) NULL else "weight")
}

get_ensemble_weightset = function(swgtmbf.forecast.values, swgtm.observation.values, forecast.types, weighting.scheme.indexer.list) {
  ## todo generalize to allow specification of which dimensions are "indices"
  ## (here: season, model week, epigroup, target spec, forecast type) vs. which
  ## are "methods" (here: backcaster & forecaster)
  cv_apply(
    swgtmbf.forecast.values,
    weighting.scheme.indexer.list,
    function(train, test) {
      fallback.method.index = "ignorant.Delphi_Uniform"
      instance.method.forecast.values.listmat =
        R.utils::wrap(train, list(1:5,6:7))
      instance.observation.values.list =
        do.call(`[`,
                c(list(swgtm.observation.values),
                  dimnames(train)[1:5])) %>>%
        {dim(.) <- NULL; .}
      if (dimp(train)[[5L]] != 1L) {
        stop ("Must break down by forecast type (each=ALL).")
      }
      forecast.type = forecast.types[[dimnamesp(train)[[5L]]]]
      fold.coefs = forecast.type[["fit_ensemble_coefs"]](
        instance.method.forecast.values.listmat ,
        instance.observation.values.list,
        ## only count separate season-location pairs as separate observations:
        prod(dim(train)[c(1L,3L)]),
        fallback.method.index
      )
      dim(fold.coefs) <- dimp(train)[6:7]
      dimnames(fold.coefs) <- dimnamesp(train)[6:7]
      fold.coefs
    },
    parallel_dim_i = 1L
  ) %>>%
    ## --- Fix up dim, dimnames: ---
    {
      old.dp = dimp(.)
      old.dnp = dimnamesp(.)
      new.dp = old.dp[-(2L + 6:7)]
      new.dnp = old.dnp[-(2L + 6:7)]
      exclude.dimension.flags = sapply(new.dnp, function(eltnames) {
        length(eltnames) == 1L && eltnames=="all"
      })
      ## . <- aperm(., c(2L + 1:3, 1:2, 6:7))
      new.dp <- new.dp[!exclude.dimension.flags]
      new.dnp <- new.dnp[!exclude.dimension.flags]
      dim(.) <- new.dp
      dimnames(.) <- new.dnp
      .
    } %>>%
    {.}
}

ensemble_and_components_linlog_plot = function(weightset, swgbf.component.target.multicasts, voxel.data, t.target.specs, m.forecast.types, s,w,g, Target) {
  swgbf.component.target.multicasts.slice = swgbf.component.target.multicasts[s,w,g,,,drop=FALSE]
  bftm.weights = weightset %>>%
    extract_partial_(c(dimnames(swgbf.component.target.multicasts.slice),
                       Target=Target,
                       dimnames(m.forecast.types)),
                     dimension.missing.behavior="ignore") %>>%
    select_dims(c("Backcaster","Forecaster","Target","Type"))
  bMtm.weights = bftm.weights
  names(dimnames(bMtm.weights))[[2L]] <- "Model"
  Mtm.weights = select_dims(bMtm.weights, 1L, "drop")
  plt = target_multicast_combined_target_plot(
    swgbf.component.target.multicasts.slice %>>% select_dims("Forecaster"),
    voxel.data,
    t.target.specs, m.forecast.types,
    Target,
    Mtm.weights
  ) +
    ggplot2::theme(legend.position="right")
  tbl = weights %>>%
    {structure(sprintf("%0.02f",.), dim=dim(.), dimnames=dimnames(.))} %>>%
    R.utils::wrap(list(1:2,3:4)) %>>%
    gridExtra::tableGrob()
  gridExtra::grid.arrange(plt, tbl, heights=c(2,1))
}

get_evaluation = function(forecast.value, observation.value, forecast.type) {
  forecast.type[["evaluate_forecast_value"]](forecast.value, observation.value)
}

save_spreadsheets =
  function(swg_.target.multicasts,
           swg.voxel.data,
           t.target.specs, m.forecast.types,
           spreadsheet.dir,
           subpath_or_NULL_for_save = function(swg.voxel.data,s,w,...) {
             spreadsheet.name = paste(s,w,...,sep=".") %>>%
               stringr::str_replace_all("/","-")
             paste0(spreadsheet.name,".csv")
           }) {
    invisible(map_join_(
      ## iterate over non-epigroup dimensions, flipping s and w for ordering purposes:
      ## arraylike.args=named_array_to_name_arrayvecs(swg_.target.multicasts)[-3L],
      ## f=function(s,w,...) {
      arraylike.args=named_array_to_name_arrayvecs(swg_.target.multicasts)[c(2L,1L,Seq(4L,ndimp(swg_.target.multicasts)))],
      f=function(w,s,...) {
        print(paste(s,w,...,sep="."))
        subpath = subpath_or_NULL_for_save(swg.voxel.data[s,w,,drop=FALSE], s,w,...)
        print(subpath)
        if (!is.null(subpath)) {
          filepath = file.path(spreadsheet.dir, subpath)
          if (!file.exists(filepath)) {
            ## get corresponding epigroup forecast tables and bind them together:
            spreadsheet =
              map_join(
                target_multicast_epigroup_forecast_table,
                swg_.target.multicasts[s,w,,...,drop=FALSE],
                swg.voxel.data[s,w,,drop=FALSE],
              no_join(t.target.specs), no_join(m.forecast.types),
              lapply_variant=lapply, shuffle=FALSE,
              progress.output=FALSE
              ) %>>%
              dplyr::bind_rows()
            dir = dirname(filepath) # allow 1 level of dir nesting within spreadsheet.dir
            if (!dir.exists(dir)) {
              dir.create(dir)
            }
            ## print(filepath)
            readr::write_csv(spreadsheet, filepath)
          }
        }
        NULL
      }, lapply_variant=lapply, shuffle=FALSE, progress.output=FALSE))
  }

save_linlog_plots =
  function(target_multicast_linlog_plotter,
           swg_.target.multicasts,
           swg.voxel.data,
           t.target.specs, m.forecast.types,
           linlog.plot.dir
           ) {
    linlog.plots = map_join(
      target_multicast_linlog_plotter,
      swg_.target.multicasts,
      swg.voxel.data,
      no_join(t.target.specs), no_join(m.forecast.types)
    )
    linlog.plot.names =
      dimnames(linlog.plots) %>>%
      expand.grid() %>>%
      {do.call(paste, c(as.list(.), list(sep=".")))} %>>%
      as.character() %>>%
      stringr::str_replace_all("/","-") %>>%
      structure(
        dim=dim(linlog.plots),
        dimnames=dimnames(linlog.plots)
      )
    if (!dir.exists(linlog.plot.dir)) {
      dir.create(linlog.plot.dir, recursive=TRUE)
    }
    invisible(map_join(
      function(linlog.plot, linlog.plot.name) {
        filepath = file.path(linlog.plot.dir, paste0(linlog.plot.name,".pdf"))
        print(filepath)
        ggplot2::ggsave(filepath, linlog.plot + ggplot2::ggtitle(linlog.plot.name))
        NULL
      },
      linlog.plots, linlog.plot.names,
      lapply_variant=lapply, shuffle=FALSE,
      progress.output=FALSE
    ))
  }

save_weighting_linlog_plots =
  function(weighset,
           swgbf.component.target.multicasts,
           swg.voxel.data,
           t.target.specs, m.forecast.types,
           plot.dir
           ) {
    if (!dir.exists(plot.dir)) {
      dir.create(plot.dir, recursive=TRUE)
    }
    swgt.plots = map_join(
      function(weightset, swgbf.component.target.multicasts, voxel.data, t.target.specs, m.forecast.types, s,w,g, Target) {
        filename = paste0("weighting_linlog_",paste(s,w,g,Target,sep=".")) %>>%
          stringr::str_replace_all("/","-")
        filepath = file.path(plot.dir, filename)
        print(filepath)
        pdf(filepath, width=10, height=10)
        ensemble_and_components_linlog_plot(weightset, swgbf.component.target.multicasts, voxel.data, t.target.specs, m.forecast.types, s,w,g, Target)
        dev.off()
      },
      no_join(weightset),
      no_join(swgbf.component.target.multicasts),
      swg.voxel.data,
      no_join(t.target.specs), no_join(m.forecast.types),
      named_array_to_name_arrayvecs(swg.voxel.data)[[1L]],
      named_array_to_name_arrayvecs(swg.voxel.data)[[2L]],
      named_array_to_name_arrayvecs(swg.voxel.data)[[3L]],
      t.target.specs %>>% named_arrayvec_to_name_arrayvec(),
      lapply_variant=lapply, shuffle=FALSE,
      progress.output=FALSE
    )
  }
