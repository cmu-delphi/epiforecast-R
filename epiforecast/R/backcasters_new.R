
na_fill0_and_indicate_with_no_interactions = function(df, which.cols=names(df)) {
    df %>>%
        dplyr::mutate_at(which.cols, list(missing=function(col) as.numeric(is.na(col)))) %>>%
        dplyr::mutate_at(which.cols, dplyr::coalesce, 0) %>>%
        dplyr::rename_at(which.cols, paste0, "_fill0") %>>%
        {.}
}

fit_quantile_coefmat_thin_whiten = function(x, y, taus, method, weights=NULL, dtolrelmaxconstant=1e-6, whitening.approach=c("ud","u"), ...) {
    whitening.approach <- match.arg(whitening.approach)
    ## xxx not actually whitening... not centering...
    ## todo center xs instead of whitening w/ intercept?
    ## xxx order of y and maybe (Intercept) is important here:
    Xy = dplyr::bind_cols(y, `(Intercept)`=rep(1,nrow(y)), x)
    ## Remove instances with NA Xs or ys:
    Xyp = na.omit(Xy)
    if (is.null(weights)) {
        weightsp = NULL
    } else if (is.null(attr(Xyp,"na.action"))) {
        weightsp = weights
    } else {
        weightsp = weights[-attr(Xyp,"na.action")]
    }
    Xpmat = as.matrix(Xyp[,-1L])
    ypvec = Xyp[[1L]]
    ## todo look into different orthogonalization, whitening, and thinning approaches
    ## Prevent singular X matrix errors in quantreg::rq fits (based on rank provided by `qr`) by working off of singular vectors with significantly nonzero values instead of original features:
    svd.Xpmat = svd(Xpmat)
    Up = svd.Xpmat[["u"]]
    dp = svd.Xpmat[["d"]]
    Vp = svd.Xpmat[["v"]]
    dtolrelmax = max(dim(Xpmat)) * dtolrelmaxconstant
    ## `qr` may use a rule with .Machine[["double.eps"]]; expand this by at least 10x to try to ensure it will not complain about various fix-ups we might pass it based on the trimmed SVD:
    if (dtolrelmax < 10*.Machine[["double.eps"]]) {
        stop ('dtolrelmaxconstant led to dtolrelmax low enough that we might not expect to avoid singular design matrix errors for some thinning approaches.')
    }
    dp.inclusion.flags = dp >= dtolrelmax*max(dp)
    stopifnot(sum(dp.inclusion.flags) >= 1L)
    Upp = Up[,dp.inclusion.flags,drop=FALSE]
    dpp = dp[dp.inclusion.flags]
    Vpp = Vp[,dp.inclusion.flags,drop=FALSE]
    Phipp = switch(whitening.approach,
                   ud=t(dpp*t(Upp)),
                   u=Upp
                   )
    yppvec = ypvec
    weightspp = weightsp
    ##
    fit_for_tau =
        if (is.null(weights)) {
            function(tau) {
                quantreg::rq.fit(Phipp, yppvec, tau, method=method, ...)[["coefficients"]]
            }
        } else {
            function(tau) {
                quantreg::rq.wfit(Phipp, yppvec, tau, weightspp, method=method, ...)[["coefficients"]]
            }
        }
    coefmat =
        taus %>>%
        sapply(fit_for_tau) %>>%
        {
            switch(whitening.approach,
                   ud = Vpp %*% .,
                   u = Vpp %*% (./dpp)
                   )
        } %>>%
        magrittr::set_rownames(colnames(Xpmat)) %>>%
        {.}
    coefmat
}

augment_with_sirs_covariates = function(df, covariate.availabilities) {
  df %>>%
    {
      if (covariate.availabilities %>>%
          {match("stable@s-1", .[["variable.name"]])} %>>%
          {is.na(.) || !covariate.availabilities[["request.available"]][[.]]}
      ) {
        . # stable@s-1 not available; can't form any SIRS-inspired covariates
      } else {
        . %>>%
          ## add the SIRS-inspired covariates that are only conditioned on 1stable@s-11
          dplyr::mutate(
                 `(stable@s-1)(stable@s-1)` = `stable@s-1`*`stable@s-1`
          ) %>>%
          {
            if (covariate.availabilities %>>%
                {match("stable@s-2", .[["variable.name"]])} %>>%
                {is.na(.) || !covariate.availabilities[["request.available"]][[.]]}
            ) {
              . # stable@s-2 not available; can't form rest of SIRS-inspired covariates
            } else {
              ## add the rest of the SIRS-inspired covariates, with care to make relative delta NA when denominator is around 0 or below.
              . %>>%
                dplyr::mutate(
                       `(stable@s-1)(stable@s-2)` = `stable@s-1`*`stable@s-2`,
                       `(stable@s-1)(reldiff@s-1)` = `stable@s-1`*(`stable@s-1`-`stable@s-2`)/dplyr::if_else(`stable@s-2` > 1e-3, `stable@s-2`, NA_real_)
                )
            }
          }
      }
    }
}

quantile_arx_thinning_whitening_pancaster = function(max.weeks.ahead, include.nowcast, include.sirs.inspired, n.sims=200L, ...) function(voxel.data, g.voxel.data, source.name, signal.name) {
  max.weeks.ahead <- match.single.nonna.integer(max.weeks.ahead)
  if (!is.logical(include.nowcast) || !is.logical(include.sirs.inspired)) {
    stop ('include.* args should be logical')
  }
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
        {
          if (!include.sirs.inspired) {
            . # SIRS-inspired covariates not enabled; return existing df
          } else {
            . %>>% augment_with_sirs_covariates(covariate.availabilities)
          }
        } %>>%
        dplyr::select(-reference.epiweek) %>>%
        dplyr::filter(!is.na(y)) %>>%
        ## todo allow for different NA handling routines
        na_fill0_and_indicate_with_no_interactions(colnames(.)[colnames(.)!="y"]) %>>%
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
      test.covariate.data =
          covariate.test.tbl %>>%
          {
              if (!include.sirs.inspired) {
                  . # SIRS-inspired covariates not enabled; return existing df
              } else {
                  . %>>% augment_with_sirs_covariates(covariate.availabilities)
              }
          } %>>%
          na_fill0_and_indicate_with_no_interactions() %>>%
          {.}
      ## todo base on deltas, etc.
      ## todo sliding window average / exponential average of observations / deltas
      taus = runif(n.sims)
      coefmat = fit_quantile_coefmat_thin_whiten(train.data[names(train.data)!="y"], train.data["y"], taus, ...)
      intercept.row.i = 1L
      stopifnot(rownames(coefmat)[[intercept.row.i]] == "(Intercept)")
      stopifnot(isTRUE(all.equal(rownames(coefmat), c("(Intercept)",colnames(test.covariate.data)))))
      simulated.values =
          coefmat[intercept.row.i,,drop=TRUE] +
          colSums(coefmat[-intercept.row.i,,drop=FALSE]*t(as.matrix(test.covariate.data)))
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

## todo weight based on season, week, #nonmissing features, feature values, ...
