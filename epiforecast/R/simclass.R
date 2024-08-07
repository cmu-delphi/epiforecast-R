## author_header begin
## Copyright (C) 2017 Sangwon Hyun, Logan C. Brooks
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

##' Make transparent colors according to weights in (0,1).
##' @param weights Numeric vector with elements in [0,1].
##' @return Character vector containing color codes.
make_sim_ys_colors = function(weights){
    stopifnot(all(weights >= 0))
    mycols = c("#E41A1C", "#377EB8", "#4DAF4A")
    myrgbs = col2rgb(mycols[1])
    shades = round(weights/max(weights)*255)
    mycols = rgb(red = myrgbs[1],
                 green = myrgbs[2],
                 blue = myrgbs[3],
                 alpha = shades,
                 maxColorValue = 255)
    return(mycols)
}

##' Plotting function for "sim" class
##'
##' @import hexbin graphics
##' @method plot sim
##' @export
##' @export plot.sim
plot.sim = function(mysim, ylab = "Disease Intensity", xlab = "Time", lty = 1,
                    nplot = min(100,ncol(mysim$ys)), type = c("lineplot","density", "hexagonal"), overlay=FALSE, ...){

    type = match.arg(type)
    ## Make some colors
    mycols = make_sim_ys_colors(mysim$weights[1:nplot])

    desired.mfrow.or.null =
      if(mysim$control.list$model=="Empirical Bayes") {
        c(2,1)
      } else {
        NULL # (single-panel)
      }

    ## Save incoming mfrow par setting to restore if needed at the end, then
    ## override if the plot for the given sim object's model type will use
    ## multiple panels:
    if (!is.null(desired.mfrow.or.null)) {
      mfg.beforehand = par("mfg")
      ## mfg is of form (row, col, #rows, #cols); completed/fresh multipane plots have row==#rows, col==#cols; ensure that this is the case if changing mfrow:
      if (!identical(mfg.beforehand[1:2], mfg.beforehand[3:4])) {
        stop ("This type of sim object uses multiple panes in its plots, but it appears that R is in the middle of building a multi-pane plot; refusing to continue.")
      }
      mfrow.beforehand = par("mfrow")
      par(mfrow=desired.mfrow.or.null)
    }

    ## Make line plot of ys
    if(type=="lineplot"){

        ## Make line plot
        matplot(mysim$ys[,1:nplot], type = 'l', ylab = ylab, xlab = xlab, col = mycols,
            lty = lty, axes = F, lwd=.5,...)
        axis(1); axis(2);

        ## Plot improvements
        time.of.forecast = get.latest.time(mysim$new.dat)
        ticks = c(1, 5, 10, 40, 150, 500, 1000)
        abline(v = time.of.forecast, lty=2, col = 'grey50')
        axis(side = 1, at = time.of.forecast, lty = 0, font = 2)
        mtext(side=3, cex=1.5, font=2, text = paste("Simulated trajectories for season:", mysim$new.season.label))
        mtext(side=3, cex=1.25, line=-1.25, text = mysim$control.list$model)
    }


    ## Hexagonal plots (in progress)
    ## Easy way (wrong, since it plots /globally/; this is a problem of the line plots too!)
    if(type=="hexagonal"){
        nplot = 100
        xx= rep(1:nrow(mysim$ys), nplot )
        yy= do.call(c,lapply(1:nplot, function(ii){mysim$ys[,ii]}))
        hbin <- hexbin(x=xx,y=yy,xbins=nplot/3)
        plot(hbin)
    }

    ## Make density plot
    if(type=="density"){
        mtext(side=3, cex=1.5, font=2, text = paste("Simulated trajectories for season:", mysim$new.season.label))
        mtext(side=3, cex=1.25, line=-1.25, text = mysim$control.list$model)
        mycols = c("#E41A1C", "#377EB8", "#4DAF4A")
        myrgbs = col2rgb(mycols[1])
        pts.in.between = seq(from=1,to=52,by=2)
        pts.in.between = c(pts.in.between,52)
        intvls.in.between = sapply(1:(length(pts.in.between)-1), function(ii){pts.in.between[ii]:pts.in.between[(ii+1)]})
        means = sapply(1:52, function(ii) weighted.mean(mysim$ys[ii,],mysim$weights))
        ## plot(NA,xlim=c(0,52),ylim=c(0,3))
        plot(NA, xlim=c(0,52),ylim=range(mysim$ys)+c(-1,+4), ylab = ylab, xlab = xlab, col = mycols,
             lty = lty, axes = F, lwd=.5,...)
        for(jj in 1:length(intvls.in.between) ){
            this.intvl.means = means[intvls.in.between[[jj]]]
            mydens = density(apply(mysim$ys[intvls.in.between[[jj]],],2,mean), weights = mysim$weights/sum(mysim$weights))
            weights = mydens$y
            weights = weights-min(weights)
            shades = round(weights/max(weights)*255)
            mycols = rgb(red = myrgbs[1],
                         green = myrgbs[2],
                         blue = myrgbs[3],
                         alpha = shades,
                         maxColorValue = 255)
            centered.dens.x = mydens$x - mean(mydens$x)
            this.intvl.dens = list(y=sapply(1:length(this.intvl.means), function(ii){this.intvl.means[ii] + centered.dens.x}),
                                   x=intvls.in.between[[jj]],
                                   col=mycols)
            if(36 %in% intvls.in.between[[jj]]){
            sapply(1:3,function(ii)print(this.intvl.means[ii]))
            }
            abline(h=1.6)
            ## plot it
            for(ii in 1:length(this.intvl.means)){
                points(x=rep(this.intvl.dens$x[ii], nrow(this.intvl.dens$y)),
                       y=this.intvl.dens$y[,ii],
                       col=this.intvl.dens$col,
                       cex=2)
            }
        }

        ## plot actual mean
        lines(means, lwd=3)

        ## Overlay the lines
        if(overlay){
            mycols = make_sim_ys_colors(mysim$weights[1:nplot])
            matlines(mysim$ys[,1:nplot], col=mycols,type='l',lty=1)
        }
        axis(1);axis(2)
        mtext(side=3, cex=1.5, font=2, text = paste("Densities of trajectories for season:", mysim$new.season.label))
        mtext(side=3, cex=1.25, line=-1.25, text = mysim$control.list$model)
    }

    ## Make barplot
    if(mysim$control.list$model=="Empirical Bayes"){
        all.hist.seasons.in.sim = mysim$old.season.labels[mysim$sampled.row.major.parms$curve.i]
        all.hist.seasons.in.sim  =      as.factor(        all.hist.seasons.in.sim)
        levels(all.hist.seasons.in.sim) = mysim$old.season.labels
        mytable = table(all.hist.seasons.in.sim)/mysim$control.list$max.n.sims*100
        barplot(mytable,col="#8DD3C7")
        title("Contribution of historical seasons to simulated curves (%)")
        ## pie(mytable, main="Historical seasons")
    }

    ## Restore original mfrow par setting:
    if (!is.null(desired.mfrow.or.null)) {
      par(mfrow=mfrow.beforehand)
    }
}




##' Printing function for "sim" class
##' @param x Output from running an OO.sim() function.
##'
##' @method print sim
##' @export
##' @export print.sim
print.sim = function(x, verbose=TRUE,...){
    mysim = x

    ctrlist = mysim[["control.list"]]
    ctrmodel = ctrlist$model
    ctrlist = ctrlist[which(names(ctrlist)!="model")]

    cat(fill=TRUE)
    cat(paste0(ncol(mysim[["ys"]]), " simulated trajectories with weights produced by ", ctrmodel, " model"), fill=TRUE)

    print.each.control.list.item = function(x){
        stopifnot(is.list(x))
        if(!is.null(x[[1]])) { cat(names(x), "=", x[[1]], fill=T) }
    }

    cat(fill=TRUE)
    cat("===========================",fill=TRUE)
    cat("Stored simulation settings:",fill=TRUE)
    cat("---------------------------",fill=TRUE)
    if (length(ctrlist) == 0L) {
      cat("(none)",fill=TRUE)
    } else {
      sapply(seq_along(ctrlist),
             function(ii) print.each.control.list.item(ctrlist[ii]))
    }
    cat("===========================",fill=TRUE)

    cat(fill=TRUE)
    cat("Default sim object components are $ys, $weights, and $control.list.",fill=TRUE)
    additional.component.names = setdiff(names(mysim), c("ys", "weights", "control.list"))
    cat("===============================",fill=TRUE)
    cat("Names of additional components:",fill=TRUE)
    cat("-------------------------------",fill=TRUE)
    if(length(additional.component.names) == 0L) {
      cat("(none)",fill=TRUE)
    } else {
      cat(paste0(sapply(additional.component.names, as.name), collapse="\n"), fill=TRUE)
    }
    cat("===============================",fill=TRUE)

    if(verbose){
        cat(fill=TRUE)
        cat("The simulated trajectories look like this:",fill=T)
        cat(fill=TRUE)
        show.sample.trajectories(ys = mysim$ys,
                                 n.shown.rows = min(5L, nrow(mysim$ys)),
                                 n.shown.cols = min(5L, ncol(mysim$ys)))
    }
}


##' Showing sample trajectories of ys.
##' @param ys A matrix whose columns contain simulated curves.
##' @param n.shown.rows How many observations per curve to show.
##' @param n.shown.cols How many curves to show.
show.sample.trajectories =
  function(ys, n.shown.rows = 100L, n.shown.cols = 100L){

    ## basic error checking
    if(!(inherits(ys, "matrix"))) stop("ys is not a matrix!")
    if(any(apply(ys,2,function(mycol){any(is.na(mycol))}))) stop("ys has missing entries!")

    ## format simulated trajectories
    has.unshown.rows = n.shown.rows < nrow(ys)
    has.unshown.cols = n.shown.cols < ncol(ys)
    n.shown.rows <- min(n.shown.rows, nrow(ys))
    n.shown.cols <- min(n.shown.cols, ncol(ys))
    row.indices = c(seq_len(n.shown.rows),
                    c(NA_integer_)[has.unshown.rows])
    col.indices = c(seq_len(n.shown.cols),
                    c(NA_integer_)[has.unshown.cols])
    ys <- as.data.frame(round(ys[row.indices, col.indices, drop=FALSE],2L))
    if (has.unshown.rows) {
      ys[n.shown.rows+1L,] = rep("...", n.shown.cols+as.integer(has.unshown.cols))
    }
    if (has.unshown.cols) {
      ys[,n.shown.cols+1L] = rep("...", n.shown.rows+as.integer(has.unshown.rows))
    }
    rownames(ys) <- c(paste("time =", seq_len(n.shown.rows)),
                      c("...")[has.unshown.rows])
    colnames(ys) <- c(paste("sim", seq_len(n.shown.cols)),
                      c("...")[has.unshown.cols])

    ## show it
    print(ys)
}

##' Forecast distribution of targets such as peak height from model fits
##'
##' This S3 method generates forecasts of various targets (e.g., peak height of
##' a trajectory) derived from some type of model fit. For example, the fit
##' model could be a collection of simulations of the trajectory, and this
##' method will calculate the desired target for each simulation and aggregate
##' these results into a distributional forecast. Refer to the appropriate S3
##' implementation for a given model fit for additional details.
##' @param fit.model the fit model (trajectory forecast, simulations, regression
##'   fit, etc.), on which to base the target forecasts
##' @seealso \code{\link{target_forecast.sim}}
##' @export
target_forecast <- function(fit.model, ...) {
  UseMethod("target_forecast", fit.model)
}

wtd.quantile.or.na =
  function (x, weights = NULL, probs = c(0, 0.25, 0.5, 0.75, 1),
            type = c("quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n"),
            na.rm = TRUE) {
  if (na.rm && all(is.na(x))) {
    ## x is empty or all NA's; return appropriate NA:
    x[NA_integer_][[1L]]
  } else {
    Hmisc::wtd.quantile(x, weights=weights, probs=probs,
                        type=type, normwt=TRUE, na.rm=na.rm)
  }
}

##' Forecast a target (peak height, etc.) using a \code{sim} object
##'
##' @import rlist
##' @param mysim Output from running an OO.sim() function.
##' @method target_forecast sim
##' @export
##' @export target_forecast.sim
target_forecast.sim = function(mysim,
                               target = c("pwk","pht","ons","dur"),
                               target.name = target,
                               target.fun = target,
                               ## todo make a property of target types?:
                               target_trajectory_preprocessor = function(trajectory) trajectory,
                               target.spec = NULL,
                               target_value_formatter = identity,
                               ## todo make a property of forecast types?
                               target.multival.behavior = c("random.val", "closest.to.pred.val"),
                               compute.estimates = TRUE,
                               hist.bins = NULL,
                               ...){

    if (missing(target) && (missing(target.name) || missing(target.fun))) {
      stop("Must supply either (a) =target= (consider 'pwk', 'pht', 'ons', and 'dur'), or (b) =target.name= AND =target.fun=.")
    }

    target.multival.behavior <- match.arg(target.multival.behavior)
    if (target.multival.behavior != "random.val") {
      stop ('Only target.multival.behavior="random.val" is currently supported.')
    }

    target.fun <- match.fun(target.fun)
    target.values = apply(mysim[["ys"]], 2L, function(trajectory) {
      processed.trajectory = target_trajectory_preprocessor(trajectory)
      target.multival = target.fun(processed.trajectory, ...)
      ## ## Optimized version of sample(target.multival, 1L):
      ## switch(length(target.multival)+1L,
      ##        stop ("target function must produce at least one value"),
      ##        target.multival[[1L]],
      ##        target.multival[[sample.int(length(target.multival), 1L, replace=TRUE, useHash=TRUE)]]
      ##        )
      target.multival
    })
    target.weights = mysim[["weights"]]
    if (is.list(target.values)) {
      target.lengths = sapply(target.values, length)
      target.weights <- rep(target.weights, target.lengths)/rep(target.lengths, target.lengths)
      target.values <- dplyr::combine(target.values)
    } else if (is.matrix(target.values)) {
      target.weights <- rep(target.weights/nrow(target.values),
                            each=nrow(target.values))
      target.values <- as.vector(target.values)
    }

    ## Compute mean, median, two-sided 90% quantiles
    if (compute.estimates) {
      estimates = list(quantile = wtd.quantile.or.na(target.values, weights=target.weights, c(0.05,0.95)),
                       quartile = wtd.quantile.or.na(target.values, weights=target.weights, c(0.25,0.5,0.75)),
                       decile   = wtd.quantile.or.na(target.values, weights=target.weights, 1:9/10),
                       mean = stats::weighted.mean(target.values, target.weights),
                       median = matrixStats::weightedMedian(target.values, target.weights))
    }

    target.settings = list(
      target_trajectory_preprocessor = target_trajectory_preprocessor,
      target.multival.behavior = target.multival.behavior,
      ...
    )

    ## Return a list of things
    settings = rlist::list.remove(unclass(mysim), c("ys","weights"))
    return (structure(c(list(settings = settings,
                             target.values = stats::setNames(list(target.values), target.name),
                             target.weights = target.weights,
                             target_value_formatter = target_value_formatter,
                             target.settings = target.settings),
                        if (compute.estimates) list(estimates = estimates)
                        else list()
                        ),
                      class="target_forecast"))
}

##' Print a \code{forecast} object
##'
##' @param x output of a \code{forecast} method
##'
##' @method print target_forecast
##' @export
##' @export print.target_forecast
print.target_forecast = function(x, sig.digit = 2L, ...) {
  target.forecast = x
  target.name = names(target.forecast[["target.values"]])
  estimates = target.forecast[["estimates"]]
  target_value_formatter = target.forecast[["target_value_formatter"]]

  ## Print a bunch of things
  cat("Summary for", target.name, ":",fill=TRUE)
  cat("====================", fill=TRUE)
  if (is.null(estimates)) {
    cat("No estimate computations recorded.", fill=TRUE)
  } else {
    cat("The mean of", target.name, "is", target_value_formatter(round(estimates$mean,sig.digit)), fill=TRUE)
    cat("The median of", target.name, "is", target_value_formatter(round(estimates$median,sig.digit)), fill=TRUE)
    cat("And the 0.05, 0.95 quantiles are", target_value_formatter(round(estimates$quantile,sig.digit)), fill=TRUE)
    cat("And the quartiles are", target_value_formatter(round(estimates$quartile,sig.digit)), fill=TRUE)
    cat("And the deciles are", target_value_formatter(round(estimates$decile,sig.digit)), fill=TRUE)
  }
}

## ##' constructor should take in the bare minimum to create forecasts
## ##' i.e. full.dat, max.n.sims, control.list, baseline, etc.
## sim <- function(...) {
##   structure(list(...), class = "sim")
## }

c_for_named_lists = function(classed.list, ..., recursive=FALSE, use.names=TRUE) {
  ## check for attempts to override set arguments:
  if (recursive != FALSE) {
    stop ("recursive must be FALSE")
  }
  if (use.names != TRUE) {
    stop ("use.names must be TRUE")
  }
  ## check that each argument is nontrivially named and/or a list with all
  ## elements nontrivially named:
  dots = list(...)
  dots.namesish = if (is.null(names(dots))) {
                    rep("", length(dots))
                  } else {
                    dplyr::coalesce(names(dots), "")
                  }
  if (any(
    dots.namesish == "" &
    sapply(dots, function(arg) {
      !is.list(arg) || is.null(names(arg)) || any(is.na(names(arg)) | names(arg) == "")
    })
  )) {
    stop ('All components to add to the object must be (a) named (and/)or (b) lists with every component named.  Names of "" are treated as invalid.')
  }
  ## forward operation to c method used for lists:
  input.as.list = unclass(classed.list)
  stopifnot(!anyDuplicated(names(input.as.list)))
  result.as.list = c(input.as.list, ..., recursive=FALSE, use.names=TRUE)
  ## check for duplicate names:
  duplicate.name.flags = duplicated(names(result.as.list))
  if (any(duplicate.name.flags)) {
    stop (paste0("All components must be uniquely named; adding the specified components would result in the following duplicates:\n",
                 paste0(sapply(names(result.as.list)[duplicate.name.flags], as.name), collapse="\n")))
  }
  ## re-class list result:
  result = structure(result.as.list, class=class(classed.list))
  return (result)
}

##' Add information to a \code{sim} object
##'
##' Uses \code{recursive=FALSE} and \code{use.names=TRUE} when forwarding to the
##' list \code{c} method; any attempt to override these values will generate an
##' error.
##'
##' @param my.sim a \code{sim} object
##' @param ... list of components to add to the \code{sim} object; must lead to
##'   a resulting \code{sim} object with components that are all uniquely,
##'   nontrivially (\code{!=""}) named
##' @return \code{sim} object with the given components appended
##'
##' @method c sim
##' @export
##' @export c.sim
c.sim = function(my.sim, ...) {
  return (c_for_named_lists(my.sim, ...))
}

##' Add information to a \code{target_forecast} object
##'
##' Uses \code{recursive=FALSE} and \code{use.names=TRUE} when forwarding to the
##' list \code{c} method; any attempt to override these values will generate an
##' error.
##'
##' @param target.forecast a \code{target_forecast} object
##' @param ... list of components to add to the \code{target_forecast} object;
##'   must lead to a resulting \code{target_forecast} object with components that
##'   are all uniquely, nontrivially (\code{!=""}) named
##' @return \code{target_forecast} object with the given components appended
##'
##' @method c target_forecast
##' @export
##' @export c.target_forecast
c.target_forecast = function(target.forecast, ...) {
  return (c_for_named_lists(target.forecast, ...))
}

##' \code{plot} method for \code{target_forecast} objects
##'
##' @param x \code{target_forecast} object
##' @param add logical, length 1, non-\code{NA}: whether to plot on top of the
##'   currently active plot (vs. creating a new plot)
##'
##' @method plot target_forecast
##' @export
##' @export plot.target_forecast
plot.target_forecast = function(x, add=FALSE, ...) {
  target.forecast = x
  target.name = names(target.forecast[["target.values"]])
  target.values = target.forecast[["target.values"]][[1L]]
  target.weights = target.forecast[["target.weights"]]
  estimates = target.forecast[["estimates"]]

  ## Plot histogram
  par(mfrow=c(1,1))
  weights::wtd.hist(target.values, weight=target.weights, axes=FALSE, main="", xlab = target.name, col="skyblue")
  mtext(paste("Forecasts for target:", target.name),3,cex=2,padj=-.5)
  axis(1); axis(2);
  abline(v=estimates$mean, col = 'red', lwd=3)
  abline(v=estimates$median, col = 'blue', lwd=2)
  abline(v=estimates$quantile, col = 'grey50', lwd=1, lty=2)
  legend("topright", col = c('red','blue','grey50'), lty=c(1,1,2), lwd = c(2,2,1), legend = c("Mean", "Median", "5%, 95% quantiles"))

  ## add lines to original plot.sim output if necessary
  if(add) {
    stop("add=TRUE functionality not written yet")
  }
}

##' Resample a sim (if necessary) to get <= max.n.sims simulated curves
##'
##' @param sim.obj a sim object
##' @param max.n.sims a single non-NA non-negative integer; the inclusive upper
##'   bound on the number of simulated curves in the result
##' @return a sim object with <= max.n.sims simulated curves; the weights in the
##'   result are the same as the weights of the corresponding curves in
##'   \code{sim.obj}, implying that the sum of the weights in the result is less
##'   than or equal to the sum of the weights in the result (due to sampling
##'   without replacement), similar to reductions in "effective particles" in
##'   particle sampling contexts due to resampling.
##' @details Resampling is only performed if the number of simulations needs to
##'   change. Resampling is done without replacement to prevent unnecessary
##'   reductions in the number of "effective particles".
downsample_sim = function(sim.obj, max.n.sims) {
  max.n.sims <- match.single.nonna.integer(max.n.sims)
  if (max.n.sims < 0L) {
    stop ("max.n.sims must be a natural number.")
  }
  input.n.sims = length(sim.obj[["weights"]])
  if (input.n.sims <= max.n.sims) {
    result = sim.obj
  } else {
    inds = sample(seq_along(sim.obj$weights), max.n.sims,
                  prob=sim.obj$weights, replace=FALSE)
    result = sim.obj # to include other things besides ys and weights
    result[c("ys","weights")] <- list(
      sim.obj$ys[,inds,drop=FALSE],
      rep(mean(sim.obj$weights), length(inds))
    )
    result[["last.actually.sampled.to.n.sims.of"]] <- max.n.sims
  }
  return (result)
}

shuffle_sim_indices = function(sim.obj) {
  shuffled.inds = sample.int(length(sim.obj[["weights"]]))
  sim.obj[c("ys","weights")] <- list(
    sim.obj$ys[,shuffled.inds,drop=FALSE],
    sim.obj$weights[shuffled.inds]
  )
  return (sim.obj)
}

##' Resample a sim (if necessary) to get >= min.n.sims simulated curves
##'
##' @param sim.obj a sim object
##' @param max.n.sims a single non-NA non-negative integer; the inclusive lower
##'   bound on the number of simulated curves in the result
##' @param inflate.weights a single non-NA logical; TRUE indicated that the
##'   weights be inflated proportionally with the increase in the number of
##'   simulations; FALSE indicates that the total weight should be preserved
##'   instead.
##' @return a sim object with >= min.n.sims simulated curves; the sum of weights
##'   in the result will equal the sum of the results in the input.
##' @details Resampling is only performed if the number of simulations needs to
##'   change. Any resampling is (necessarily) done with replacement.
upsample_sim = function(sim.obj, min.n.sims, inflate.weights) {
  min.n.sims <- match.single.nonna.integer(min.n.sims)
  if (min.n.sims < 0L) {
    stop ("min.n.sims must be a natural number.")
  }
  input.n.sims = length(sim.obj[["weights"]])
  if (input.n.sims >= min.n.sims) {
    result = sim.obj
  } else {
    sample.inds = sample(seq_along(sim.obj$weights), min.n.sims-input.n.sims,
                         prob=sim.obj$weights, replace=TRUE)
    mean.result.weight =
      if (inflate.weights) {
        mean(sim.obj[["weights"]])
      } else {
        sum(sim.obj[["weights"]])/min.n.sims
      }
    result = sim.obj # to include other things besides ys and weights
    result[c("ys","weights")] <- list(
      sim.obj$ys[,c(seq_len(input.n.sims),
                    sample.inds), drop=FALSE],
      c(sim.obj[["weights"]]*mean.result.weight/mean(sim.obj[["weights"]]),
        rep(mean.result.weight, length(sample.inds)))
    )
    ## shuffle the results (e.g., to prevent dependencies in case this is paired
    ## with another upsampled sim object):
    result <- shuffle_sim_indices(result)
    result[["last.actually.sampled.to.n.sims.of"]] <- min.n.sims
  }
  return (result)
}

##' @export
simplified_simlike = function(fit.model) {
  UseMethod("simplified_simlike", fit.model)
}

##' @export
simplified_simlike.default = function(fit.model) {
  fit.model
}

##' @method simplified_simlike sim
##' @export
##' @export simplified_simlike.sim
simplified_simlike.sim = function(fit.model) {
  sim = fit.model
  stopifnot(ncol(sim[["ys"]]) == length(sim[["weights"]]))
  structure(
    list(
      ys = sim[["ys"]],
      weights = sim[["weights"]],
      control.list = list(
        model = "simplified",
        original.model = sim[["control.list"]][["model"]],
        original.n.sims = length(sim[["weights"]]),
        original.weight.sum = sum(sim[["weights"]])
      )
    ),
    class = "sim"
  ) %>>%
    downsample_sim(20L)
}
