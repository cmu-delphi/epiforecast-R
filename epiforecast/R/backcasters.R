
## backfill_ignorant_backsim = function(voxel.data, signal.name) {
backfill_ignorant_backsim = function(voxel.data, g.voxel.data, source.name, signal.name) {
  ## old.dat = voxel.data[["epidata.df"]] %>>%
  old.dat = voxel.data[["epidata.dfs"]][[source.name]] %>>%
    dplyr::filter(season != voxel.data[["season"]]) %>>%
    split(.[["season"]]) %>>%
    magrittr::extract(
                sapply(., function(season.df) {
                  !any(is.na(season.df[[signal.name]]))
                })
              ) %>>%
    dplyr::bind_rows() %>>%
    {split(.[[signal.name]], .[["Season"]])}
  ## new.dat = voxel.data[["epidata.df"]] %>>%
  new.dat = voxel.data[["epidata.dfs"]][[source.name]] %>>%
    dplyr::filter(season == voxel.data[["season"]]) %>>%
    magrittr::extract2(signal.name)
  new.dat.sim = match.new.dat.sim(new.dat)
  voxel.Season = season_to_Season(
    voxel.data[["season"]], voxel.data[["first.week.of.season"]])
  ## concatenate new.dat.sim onto old.dat, setting the :
  full.dat = c(old.dat,
               setNames(list(new.dat.sim), voxel.Season))
  return (full.dat)
}

quantile_arx_pancaster = function(include.nowcast, max.weeks.ahead, lambda=1e-3, tol=1e-3, n.sims=200L) function(voxel.data, g.voxel.data, source.name, signal.name) {
  model.week.shift.range = c(-4L, +4L)
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
          dplyr::filter(order(-lag)==1L) %>>%
          dplyr::ungroup() %>>%
          dplyr::group_by(epiweek) %>>%
          dplyr::mutate(lag.groups.from.latest=order(-lag.group)-1L) %>>%
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
    ## DatesOfSeason(voxel.data[["season"]], usa.flu.first.week.of.season, 0L,3L) %>>%
    ## magrittr::extract2(1L) %>>%
    ## Date_to_epiweek()
    g.voxel.data[[target.epigroup]][["epidata.dfs"]][[target.source.name]] %>>%
    dplyr::filter(season == target.season) %>>%
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
        dplyr::filter(dplyr::between(
        (reference.epiweek %>>% epiweek_to_sunday() %>>% DateToYearWeekWdayDF(0L,3L) %>>% yearWeekDFToSeasonModelWeekDF(voxel.data[["first.week.of.season"]], 3L))[["model.week"]] -
        (pancast.epiweek %>>% epiweek_to_sunday() %>>% DateToYearWeekWdayDF(0L,3L) %>>% yearWeekDFToSeasonModelWeekDF(voxel.data[["first.week.of.season"]], 3L))[["model.week"]],
        model.week.shift.range[[1L]], model.week.shift.range[[2L]])) %>>%
        {.}
      train.data =
        full.train.tbl %>>%
        dplyr::select(-reference.epiweek) %>>%
        na.omit() %>>%
        ## Throw out features that appear to be redundant (preventing singular
        ## matrix errors later):
        {
          lm.fit = lm(y~., ., tol=tol)
          beta = coefficients(lm.fit)
          stopifnot(names(beta)[[1L]] == "(Intercept)" &&
                    length(beta) - 1L == ncol(.) - 1L &&
                    names(.)[[ncol(.)]] == "y")
          .[c(!is.na(beta)[-1L],TRUE)]
        } %>>%
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
      ## fit = lm(y~., train.data)
      ## print(train.data)
      ## sigma = sqrt(mean(residuals(fit)^2))
      ## print(fit)
      ## print(anova(fit))
      ## print(sigma)
      ## lm.fit = lm(y~., train.data)
      ## print(anova(lm.fit))
      ## taus = runif(n.sims)
      ## simulated.values =
      ##   sapply(seq_len(n.sims), function(sim.i) {
      ##     ## fit = quantreg::rq(y~., taus[[sim.i]], train.data)
      ##     fit = quantreg::rq(y~., taus[[sim.i]], train.data, method="lasso", lambda=1e-6)
      ##     predict(fit, newdata=covariate.test.tbl[sim.i,,drop=FALSE])
      ##     ## predict(fit, newdata=covariate.test.tbl[sim.i,,drop=FALSE]) + rnorm(1L, sigma)
      ##   }) %>>% unname()
      ## xxx calculates n.sims^2 outputs instead f just n.sims!... but faster than above for current n.sims...
      ## todo base on deltas, etc.
      ## todo sliding window average / exponential average of observations / deltas
      taus = runif(n.sims)
      fit = quantreg::rq(y~., taus, train.data, method="lasso", lambda=lambda)
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
    dplyr::filter(season == target.season) %>>%
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
  old.dat = voxel.data[["epidata.dfs"]][[source.name]] %>>%
    dplyr::filter(season != voxel.data[["season"]]) %>>%
    split(.[["season"]]) %>>%
    magrittr::extract(
                sapply(., function(season.df) {
                  !any(is.na(season.df[[signal.name]]))
                })
              ) %>>%
    dplyr::bind_rows() %>>%
    {split(.[[signal.name]], .[["Season"]])}
  voxel.Season = season_to_Season(
    voxel.data[["season"]], voxel.data[["first.week.of.season"]])
  full.dat = c(old.dat,
               setNames(list(new.dat.sim), voxel.Season))
  return (full.dat)
}
