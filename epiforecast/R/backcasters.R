
chop_by_season = function(epidata.df) {
  epidata.df %>>%
    ## tidyr::chop(-dplyr::one_of("Season","season")) %>>%
    ## todo Above does not work as Season is not set properly for the unobserved portion of the latest season in the epidata df.  Should fix this and/or move to remove the concept of an epidata df with fill-in for missing portions of seasons.
    tidyr::chop(-season) %>>%
    ## xxx to try not to break old code, keep sorting from old implementation that used split:
    dplyr::arrange(season) %>>%
    {.}
}

chop_and_extend_by_season = function(first_epiweek_of_season, last_epiweek_of_season, custom_season_to_Season) function(epidata.df) {
    season.epiweek.df = epidata.df %>>%
        dplyr::group_by(season) %>>%
        dplyr::slice(1L) %>>%
        dplyr::select(season, first.model.week=model.week) %>>%
        dplyr::group_by(season, first.model.week) %>>%
        dplyr::mutate(epiweek=list(epiweek_Seq(first_epiweek_of_season(season), last_epiweek_of_season(season)))) %>>%
        dplyr::mutate(model.week = list(first.model.week + seq_along(epiweek[[1L]]) - 1L)) %>>%
        dplyr::ungroup() %>>%
        dplyr::select(-first.model.week) %>>%
        tidyr::unchop(-season) %>>%
        {.}
    epidata.df %>>%
        dplyr::select(-dplyr::one_of("Season","season","model.week")) %>>%
        dplyr::right_join(season.epiweek.df, by="epiweek") %>>%
        ## xxx to try not to break old code, keep sorting from old implementation that used split:
        dplyr::mutate(Season = custom_season_to_Season(season)) %>>%
        chop_by_season() %>>%
        {.}
}

backfill_ignorant_backsim = function(voxel.data, g.voxel.data, source.name, signal.name, epidata_df_to_chopped_trajectory_df=chop_by_season) {
  source.epidata.df = voxel.data[["epidata.dfs"]][[source.name]]
  source.chopped.trajectory.df = epidata_df_to_chopped_trajectory_df(source.epidata.df)
  old.dat = source.chopped.trajectory.df %>>%
    dplyr::filter(season != voxel.data[["season"]]) %>>%
    dplyr::filter(vapply(.[[signal.name]], function(trajectory) !any(is.na(trajectory)), logical(1L))) %>>%
    ## xxx regenerating Season labels, indicates suboptimal interface; see additional notes in season chopper above
    {stats::setNames(.[[signal.name]], season_to_Season(.[["season"]], voxel.data[["first.week.of.season"]]))}
  new.dat = source.chopped.trajectory.df %>>%
    dplyr::filter(season == voxel.data[["season"]]) %>>%
    {.[[signal.name]][[1L]]}
  new.dat.sim = match.new.dat.sim(new.dat)
  voxel.Season = season_to_Season(
    voxel.data[["season"]], voxel.data[["first.week.of.season"]])
  ## concatenate new.dat.sim onto old.dat, setting the :
  full.dat = c(old.dat,
               setNames(list(new.dat.sim), voxel.Season))
  return (full.dat)
}

backfill_ignorant_student_t2_nowcast_backcaster = function(n.sims=1000L) function(voxel.data, g.voxel.data, source.name, signal.name, epidata_df_to_chopped_trajectory_df=chop_by_season) {
    nowcast.epidata.df = voxel.data[["epidata.dfs"]][["nowcast"]]
    nowcast.chopped.trajectory.df = epidata_df_to_chopped_trajectory_df(nowcast.epidata.df)
    nowcast.new.trajectory.df = nowcast.chopped.trajectory.df %>>%
        dplyr::filter(season == voxel.data[["season"]]) %>>%
        tidyr::unchop(-season)
    nowcastless.full.dat = backfill_ignorant_backsim(voxel.data, g.voxel.data, source.name, signal.name, epidata_df_to_chopped_trajectory_df=epidata_df_to_chopped_trajectory_df)
    nowcastless.new.dat = nowcastless.full.dat %>>%
        {.[[length(.)]]}
    nowcast.time = which(
            is.na(nowcastless.new.dat[["ys"]]) &
            !is.na(nowcast.new.trajectory.df[["value"]])
        )
    if (length(nowcast.time) != 1L) {
      stop ('Expected nowcast to be available where signal was not in exactly one week.')
    }
    nowcast.value = nowcast.new.trajectory.df[["value"]][[nowcast.time]]
    nowcast.std = nowcast.new.trajectory.df[["std"]][[nowcast.time]]
    new.dat = nowcastless.new.dat %>>%
        upsample_sim(n.sims, inflate.weights=TRUE) %>>%
        {.[["ys"]][nowcast.time,] <- nowcast.value+nowcast.std*rt(n.sims, 2L); .}
    full.dat = nowcastless.full.dat %>>%
        {.[[length(.)]] <- new.dat; .}
    return (full.dat)
}

quantile_arx_pancaster = function(include.nowcast, max.weeks.ahead, lambda=1e-3, tol=1e-3, model.week.shift.range=c(-4L,+4L), include.intercept=TRUE, n.sims=200L) function(voxel.data, g.voxel.data, source.name, signal.name, epidata_df_to_chopped_trajectory_df=chop_by_season) {
  target.epigroup = voxel.data[["epigroup"]]
  target.source.name = voxel.data[["source.name"]]
  target.signal.name = voxel.data[["signal.name"]]
  target.season = voxel.data[["season"]]
  g.observed.latest = map_join(
    g.voxel.data,
    f=function(voxel.data) {
      lapply(voxel.data[["epidata.dfs"]], function(epidata.df) {
        epidata.df %>>%
          dplyr::group_by(epiweek, lag.group) %>>%
          dplyr::arrange(-lag) %>>%
          dplyr::filter(seq_along(lag)==1L) %>>%
          dplyr::ungroup()
      })
    },
    lapply_variant=lapply, shuffle=FALSE, show.progress=FALSE
  )
  g.observed.history = map_join(
      g.voxel.data,
      f=function(voxel.data) {
          lapply(voxel.data[["epidata.history.dfs"]], function(epidata.history.df) {
              epidata.history.df %>>%
                  dplyr::group_by(epiweek, lag.group) %>>%
                  dplyr::arrange(-lag) %>>%
                  dplyr::filter(seq_along(lag)==1L) %>>%
                  dplyr::ungroup() %>>%
                  dplyr::group_by(epiweek) %>>%
                  dplyr::mutate(lag.groups.from.latest=rank(-lag.group,,"min")-1L) %>>%
                  dplyr::ungroup()
          })
      },
      lapply_variant=lapply, shuffle=FALSE, show.progress=FALSE
  )
  request_availabilities = function(requests, g.latest, g.history, reference.epiweek) {
    requests %>>%
      dplyr::rowwise() %>>%
      dplyr::do({
        request = .
        availability =
          tibble::tibble(
                    epiweek = add_epiweek_integer(reference.epiweek, request[["relative.epiweek"]]),
                    lag.groups.from.latest = request[["lag.groups.from.latest"]]
                  ) %>>%
          {
            if (request[["use.lag.group"]]) {
              dplyr::left_join(., g.history[[request[["epigroup"]]]][[request[["source.name"]]]],
                               by=c("epiweek","lag.groups.from.latest"))
            } else {
              dplyr::left_join(., g.latest[[request[["epigroup"]]]][[request[["source.name"]]]],
                               by="epiweek") %>>%
                dplyr::mutate(lag.group=NA_integer_)
            }
          } %>>%
          dplyr::mutate(variable.name=request[["variable.name"]]) %>>%
          dplyr::rename_(.dots=c("request.value"=as.name(request[["signal.name"]]))) %>>%
          {
            if (!is.list(.[["request.value"]])) {
              dplyr::mutate(., request.value=lapply(request.value, function(val) {
                if (is.na(val)) NULL
                else val
              }))
            } else {
              .
            }
          } %>>%
          dplyr::mutate_at(dplyr::vars(request.value), as.list) %>>%
          dplyr::transmute(variable.name,
                           request.value,
                           request.available=
                             ## if (is.list(request.value)) !sapply(request.value, is.null)
                             ## else !is.na(request.value),
                             !sapply(request.value, is.null),
                           relative.epiweek=request[["relative.epiweek"]],
                           use.lag.group=request[["use.lag.group"]], lag.group,
                           signal.name=request[["signal.name"]], epigroup=request[["epigroup"]], source.name=request[["source.name"]]) %>>%
          {.}
        ## todo require only one target epiweek at a time here
        if (nrow(availability) != 1L) {
          stop ("Multiple matches for request.")
        }
        availability
      }) %>>%
      dplyr::ungroup()
  }
  ## todo grab observed data to replace missing sim
  description_train_tbl = function(descriptions,
                                   g.latest, g.history, simulation.descriptions,
                                   reference.epiweeks) {
    train.tbl = tibble::tibble(reference.epiweek=
                                 if (is.null(reference.epiweeks)) integer(0L)
                                 else reference.epiweeks)
    for (availability.i in seq_len(nrow(descriptions))) {
      availability = descriptions[availability.i,]
      obs.availability =
        if (availability[["source.name"]]=="simulations") {
          if (availability[["use.lag.group"]]) {
            stop ("use.lag.group=TRUE not supported with simulations source.")
          }
          availability %>>%
            dplyr::transmute(variable.name=signal.name, request.available) %>>%
            dplyr::left_join(simulation.descriptions, by="variable.name") %>>%
            dplyr::transmute(variable.name, request.available,
                             relative.epiweek = relative.epiweek + availability[["relative.epiweek"]],
                             use.lag.group, lag.group,
                             signal.name, epigroup, source.name)
          ## todo if sim is available, maybe should use it in training; but this
          ## requires retraining or more complex fitting methods
        } else {
          availability
        }
      var.tbl =
        (if (obs.availability[["use.lag.group"]]) {
           g.history[[obs.availability[["epigroup"]]]][[obs.availability[["source.name"]]]] %>>%
             dplyr::filter(lag.group==obs.availability[["lag.group"]])
         } else {
           g.latest[[obs.availability[["epigroup"]]]][[obs.availability[["source.name"]]]] %>>%
             dplyr::mutate(lag.group=NA_integer_)
         }) %>>%
        dplyr::rename_(.dots=stats::setNames(list(as.name(obs.availability[["signal.name"]])), availability[["variable.name"]])) %>>%
        dplyr::mutate(reference.epiweek=add_epiweek_integer(epiweek, -obs.availability[["relative.epiweek"]])) %>>%
        magrittr::extract(c("reference.epiweek",availability[["variable.name"]])) %>>%
        {.}
      train.tbl <- if (is.null(reference.epiweeks)) {
                     dplyr::full_join(train.tbl, var.tbl, by="reference.epiweek")
                   } else {
                     dplyr::left_join(train.tbl, var.tbl, by="reference.epiweek")
                   }
    }
    return (train.tbl)
  }
  covariate.requests =
    tibble::tribble(
              ~variable.name, ~relative.epiweek, ~use.lag.group, ~lag.groups.from.latest,
              ~signal.name, ~epigroup, ~source.name,
              "latest@s"   ,  0L,  TRUE,          0L, target.signal.name, target.epigroup, target.source.name,
              "nowcast@s"  ,  0L,  TRUE,          0L,            "value", target.epigroup,          "nowcast",
              "latest@s+1" , +1L,  TRUE,          0L, target.signal.name, target.epigroup, target.source.name,
              "latest@s-1" , -1L,  TRUE,          0L, target.signal.name, target.epigroup, target.source.name,
              "2ndlatest@s",  0L,  TRUE,          1L, target.signal.name, target.epigroup, target.source.name,
              "stable@s-1" , -1L, FALSE, NA_integer_,         "stable@s", target.epigroup,      "simulations",
              "stable@s-2" , -2L, FALSE, NA_integer_,         "stable@s", target.epigroup,      "simulations",
              "stable@s-3" , -3L, FALSE, NA_integer_,         "stable@s", target.epigroup,      "simulations",
              "stable@s-4" , -4L, FALSE, NA_integer_,         "stable@s", target.epigroup,      "simulations"
            ) %>>%
    dplyr::filter(include.nowcast | variable.name != "nowcast@s")
  ## todo especially if regularizing, may want explicit `latestcorrection@s-1` = `stable@s-1`-`latest@s-1` to model connection between corrections within same issue of nearby epiweeks
  ## todo add cross-region covariates
  simulation.descriptions =
    tibble::tribble(
              ~variable.name, ~relative.epiweek, ~use.lag.group, ~lag.group,
              ~signal.name, ~epigroup, ~source.name,
              "stable@s", 0L, FALSE, NA_integer_, target.signal.name, target.epigroup, target.source.name
            )
  pancast.epiweeks =
    g.voxel.data[[target.epigroup]][["epidata.dfs"]][[target.source.name]] %>>%
    epidata_df_to_chopped_trajectory_df() %>>%
    dplyr::filter(season == target.season) %>>%
    tidyr::unchop(-season) %>>%
    dplyr::filter(epiweek <= add_epiweek_integer(max(epiweek[!is.na(issue)]), max.weeks.ahead)) %>>%
    magrittr::extract2("epiweek") %>>%
    sort() %>>%
    {.}
  g.obs.sim.history = g.observed.history
  g.obs.sim.latest = g.observed.latest %>>%
    lapply(c, list(simulations=tibble::tibble(epiweek=integer(0L), lag.group=integer(0L))))
  for (response.i in seq_len(nrow(simulation.descriptions))) {
    simulation.description = simulation.descriptions[response.i,]
    if (simulation.description[["use.lag.group"]]) {
      stop ("use.lag.group=TRUE response requests not supported yet")
    } else {
      g.obs.sim.latest[[simulation.description[["epigroup"]]]][["simulations"]] <-
        g.obs.sim.latest[[simulation.description[["epigroup"]]]][["simulations"]] %>>%
        dplyr::bind_cols(stats::setNames(tibble::tibble(rep(list(NULL), nrow(.))),
                                         simulation.description[["variable.name"]]))
    }
  }
  fitting.formula = if (include.intercept) {
                        y ~ .
                    } else {
                        y ~ . + 0
                    }
  for (response.i in seq_len(nrow(simulation.descriptions))) {
    for (pancast.epiweek in pancast.epiweeks) {
      ## print(pancast.epiweek)
      covariate.availabilities =
        request_availabilities(covariate.requests, g.obs.sim.latest, g.obs.sim.history, pancast.epiweek) %>>%
        dplyr::filter(request.available) %>>%
        ## dplyr::group_by(relative.epiweek, use.lag.group, lag.group, signal.name, epigroup, source.name) %>>%
        ## dplyr::filter(seq_along(variable.name)==1L) %>>%
        ## dplyr::ungroup() %>>%
        {.}
      covariate.train.tbl =
        covariate.availabilities %>>%
        description_train_tbl(g.obs.sim.latest, g.obs.sim.history, simulation.descriptions, NULL) %>>%
        {.}
      response.description = simulation.descriptions[response.i,]
      response.train.tbl = response.description %>>%
        description_train_tbl(g.obs.sim.latest, g.obs.sim.history, simulation.descriptions, NULL)
      full.train.tbl = dplyr::full_join(covariate.train.tbl, response.train.tbl, by="reference.epiweek") %>>%
        dplyr::rename_(.dots=c("y"=as.name(response.description[["variable.name"]]))) %>>%
        {
            if (is.null(model.week.shift.range)) {
                .
            } else {
                . %>>%
                    dplyr::filter(dplyr::between(
                    (reference.epiweek %>>% epiweek_to_sunday() %>>% DateToYearWeekWdayDF(0L,3L) %>>% yearWeekDFToSeasonModelWeekDF(voxel.data[["first.week.of.season"]], 3L))[["model.week"]] -
                    (pancast.epiweek %>>% epiweek_to_sunday() %>>% DateToYearWeekWdayDF(0L,3L) %>>% yearWeekDFToSeasonModelWeekDF(voxel.data[["first.week.of.season"]], 3L))[["model.week"]],
                    model.week.shift.range[[1L]], model.week.shift.range[[2L]])) %>>%
                    {.}
            }
        } %>>%
        {.}
      train.data =
        full.train.tbl %>>%
        dplyr::select(-reference.epiweek) %>>%
        {
          ## filter out columns with too many NA's (also considering NA's in
          ## previously okay'd columns; earlier columns are favored for
          ## inclusion over later ones; "y" col must be included):
          original.nrow = nrow(.)
          min.nrow = max(10, original.nrow*0.10)
          running.dat = .[!is.na(.[["y"]]),]
          if (nrow(running.dat) < min.nrow) {
            stop ("Not enough non-NA y's for training.")
          }
          included.colnames = character(0L) # "y" should be grabbed at end
          for (colname in colnames(.)) {
              col.available = !is.na(running.dat[[colname]])
              if (sum(col.available) >= min.nrow) {
                  included.colnames <- c(included.colnames, colname)
                  running.dat <- running.dat[col.available,]
              }
          }
          ## print(included.colnames)
          running.dat[included.colnames]
        } %>>%
        na.omit() %>>%
        ## Throw out features that appear to be redundant (preventing singular
        ## matrix errors later):
        {
          ## todo lm is slower than, e.g., .lm.fit; cut formula processing out
          lm.fit = lm(fitting.formula, ., tol=tol)
          beta = coefficients(lm.fit)
          if (include.intercept) {
            stopifnot(names(beta)[[1L]] == "(Intercept)" &&
                      length(beta) - 1L == ncol(.) - 1L &&
                      names(.)[[ncol(.)]] == "y")
            .[c(!is.na(beta)[-1L],TRUE)]
          } else {
            stopifnot(length(beta) == ncol(.) - 1L &&
                      names(.)[[ncol(.)]] == "y")
            .[c(!is.na(beta),TRUE)]
          }
        } %>>%
        ## Try to avoid more insidious sources of singular matrices in quantile
        ## regression with jitter:
        dplyr::mutate_at(dplyr::vars(dplyr::everything()),
                         function(col) col+rnorm(length(col),,sd(col)*1e-3)) %>>%
        {.}
      covariate.test.tbl =
        covariate.availabilities %>>%
        {stats::setNames(.[["request.value"]], .[["variable.name"]])} %>>%
        lapply(function(val) {
          if (length(val)==1L) rep(val, n.sims)
          else val
        }) %>>%
        tibble::as_tibble() %>>%
        {.}
      ## fit = lm(fitting.formula, train.data)
      ## print(train.data)
      ## sigma = sqrt(mean(residuals(fit)^2))
      ## print(fit)
      ## print(anova(fit))
      ## print(sigma)
      ## lm.fit = lm(fitting.formula, train.data)
      ## print(anova(lm.fit))
      ## taus = runif(n.sims)
      ## simulated.values =
      ##   sapply(seq_len(n.sims), function(sim.i) {
      ##     ## fit = quantreg::rq(fitting.formula, taus[[sim.i]], train.data)
      ##     fit = quantreg::rq(fitting.formula, taus[[sim.i]], train.data, method="lasso", lambda=1e-6)
      ##     predict(fit, newdata=covariate.test.tbl[sim.i,,drop=FALSE])
      ##     ## predict(fit, newdata=covariate.test.tbl[sim.i,,drop=FALSE]) + rnorm(1L, sigma)
      ##   }) %>>% unname()
      ## xxx calculates n.sims^2 outputs instead f just n.sims!... but faster than above for current n.sims...
      ## todo base on deltas, etc.
      ## todo sliding window average / exponential average of observations / deltas
      taus = runif(n.sims)
      fit = tryCatch(
          ## todo quantreg::rq is slower than underlying fitting functions; cut formula processing out
          ## todo try hqreg
          ## fixme does rq w/ method lasso only use one tau?!?
          quantreg::rq(fitting.formula, taus, train.data, method="lasso", lambda=lambda),
          error=function(e) {
              ## todo iterate through taus, weight on recent data
              quantreg::rq(fitting.formula, taus, train.data, method="fn")
          }
      )
      simulated.values = predict(fit, newdata=covariate.test.tbl)[cbind(seq_len(n.sims),seq_len(n.sims))]
      g.obs.sim.latest[[response.description[["epigroup"]]]][["simulations"]] <-
        dplyr::bind_rows(g.obs.sim.latest[[response.description[["epigroup"]]]][["simulations"]],
                         tibble::tribble(~epiweek, ~lag.group, ~y,
                                         pancast.epiweek, NA_integer_, simulated.values) %>>%
                         dplyr::rename_(.dots=stats::setNames("y", response.description[["variable.name"]]))
                         )
    }
  }
  sim.obj.epiweeks =
    g.voxel.data[[target.epigroup]][["epidata.dfs"]][[target.source.name]] %>>%
    epidata_df_to_chopped_trajectory_df() %>>%
    dplyr::filter(season == target.season) %>>%
    tidyr::unchop(-season) %>>%
    magrittr::extract2("epiweek") %>>%
    {.}
  pancast.ys =
    tibble::tibble(epiweek=pancast.epiweeks) %>>%
    dplyr::left_join(g.obs.sim.latest[[target.epigroup]][["simulations"]], by="epiweek") %>>%
    magrittr::extract2("stable@s") %>>%
    sapply(identity) %>>%
    t() %>>%
    {.}
  sim.obj.ys = matrix(NA_real_, length(sim.obj.epiweeks), n.sims)
  sim.obj.ys[match(pancast.epiweeks, sim.obj.epiweeks),] <- pancast.ys
  new.dat.sim = match.new.dat.sim(sim.obj.ys)
  source.epidata.df = voxel.data[["epidata.dfs"]][[source.name]]
  source.chopped.trajectory.df = epidata_df_to_chopped_trajectory_df(source.epidata.df)
  old.dat = source.chopped.trajectory.df %>>%
      dplyr::filter(season != voxel.data[["season"]]) %>>%
      dplyr::filter(vapply(.[[signal.name]], function(trajectory) !any(is.na(trajectory)), logical(1L))) %>>%
      ## xxx regenerating Season labels, indicates suboptimal interface; see additional notes in season chopper above
      {stats::setNames(.[[signal.name]], season_to_Season(.[["season"]], voxel.data[["first.week.of.season"]]))}
  voxel.Season = season_to_Season(
    voxel.data[["season"]], voxel.data[["first.week.of.season"]])
  full.dat = c(old.dat,
               setNames(list(new.dat.sim), voxel.Season))
  return (full.dat)
}
