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
##' @import RColorBrewer plotrix
##' @param weights Numeric vector with elements in [0,1].
##' @return Character vector containing color codes.
make_sim_ys_colors = function(weights){
    stopifnot(all(weights <= 1 & weights >= 0))
    mycols = RColorBrewer::brewer.pal(3,"Set1")
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

    if(mysim$control.list$model=="Empirical Bayes") par(mfrow=c(2,1))

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
        mycols = RColorBrewer::brewer.pal(3,"Set1")
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
        mytable = table(all.hist.seasons.in.sim)/mysim$control.list$n.sims*100
        barplot(mytable,col=RColorBrewer::brewer.pal(3,"Set3")[1])
        title("Contribution of historical seasons to simulated curves (%)")
        ## pie(mytable, main="Historical seasons")
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
    cat(paste0(ctrlist$n.sims, " simulated trajectories with weights produced by ", ctrmodel, " model"), fill=T)
               ## ", to produce ", ctrlist$n.sims), fill=TRUE)
    cat(fill=TRUE)
    cat("====================",fill=TRUE)
    cat("Simulation settings:",fill=TRUE)
    cat("--------------------",fill=TRUE)

    print.each.control.list.item = function(x){
        stopifnot(is.list(x))
        if(!is.null(x[[1]])) { cat(names(x), "=", x[[1]], fill=T) }
    }

    if(ctrmodel %in% c("Empirical Bayes", "Basis Regression")){
        invisible(sapply(seq_along(ctrlist),
                         function(ii) print.each.control.list.item(ctrlist[ii])))
    } else {
        stop(paste("print not written for --", ctrmodel, " -- yet!"))
    }
    cat("====================",fill=TRUE)

    if(verbose){
        cat(fill=TRUE)
        cat("The simulated trajectories look like this:",fill=T)
        cat(fill=TRUE)
        show.sample.trajectories(ys = mysim$ys, nshow = min(5,ncol(mysim$ys)))
    }

  additional.component.names = setdiff(names(mysim), c("ys", "weights", "control.list"))
  if(length(additional.component.names) != 0L){
    cat(fill=TRUE)
    cat("===============================",fill=TRUE)
    cat("Names of additional components:",fill=TRUE)
    cat("-------------------------------",fill=TRUE)
    cat(paste0(sapply(additional.component.names, as.name), collapse="\n"), fill=TRUE)
    cat("===============================",fill=TRUE)
  }
}


##' Showing sample trajectories of ys.
##' @param ys A matrix whose columns contain simulated curves.
##' @param nshow How many curves to show.
show.sample.trajectories = function(ys, nshow = 100){

    ## basic error checking
    if(!(class(ys) %in% c("matrix"))) stop("ys is not a matrix!")
    if(any(apply(ys,2,function(mycol){any(is.na(mycol))}))) stop("ys has missing entries!")

    ## format simulated trajectories
    ys = as.data.frame(round(ys[1:(nshow+1),1:(nshow+1)],2))
    ys[(nshow+1),] = ys[,(nshow+1)] = rep("...",(nshow+1))
    colnames(ys) = Map(paste,rep("sim",(nshow+1)), 1:(nshow+1))
    rownames(ys) = Map(paste0,rep("time = ",(nshow+1)),1:(nshow+1))

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
##' @seealso \code{\link{forecast.sim}}
##' @export
target_forecast <- function(fit.model, ...) {
  UseMethod("target_forecast", fit.model)
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
                               target.calculation.digits = Inf,
                               target.value.formatter = identity,
                               ## todo make a property of forecast types?
                               target.multival.behavior = c("random.val", "closest.to.pred.val"),
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
    target.values = apply(round(mysim[["ys"]], target.calculation.digits), 2L, function(trajectory) {
      target.multival = target.fun(trajectory, ...)
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
    estimates = list(quantile = Hmisc::wtd.quantile(target.values, weights=target.weights, c(0.05,0.95)),
                     quartile = Hmisc::wtd.quantile(target.values, weights=target.weights, c(0.25,0.5,0.75)),
                     decile   = Hmisc::wtd.quantile(target.values, weights=target.weights, 1:9/10),
                     mean = stats::weighted.mean(target.values, target.weights),
                     median = matrixStats::weightedMedian(target.values, target.weights))

    target.settings = list(
      target.calculation.digits = target.calculation.digits,
      target.multival.behavior = target.multival.behavior,
      target.value.formatter = target.value.formatter,
      ...
    )

    ## Return a list of things
    settings = rlist::list.remove(unclass(mysim), c("ys","weights"))
    return (structure(list(settings = settings,
                           target.values = stats::setNames(list(target.values), target.name),
                           target.weights = target.weights,
                           estimates = estimates,
                           target.settings = target.settings),
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
  target.value.formatter = target.forecast[["target.settings"]][["target.value.formatter"]]

  ## Print a bunch of things
  cat("Summary for", target.name, ":",fill=TRUE)
  cat("====================", fill=TRUE)
  cat("The mean of", target.name, "is", target.value.formatter(round(estimates$mean,sig.digit)), fill=TRUE)
  cat("The median of", target.name, "is", target.value.formatter(round(estimates$median,sig.digit)), fill=TRUE)
  cat("And the 0.05, 0.95 quantiles are", target.value.formatter(round(estimates$quantile,sig.digit)), fill=TRUE)
  cat("And the quartiles are", target.value.formatter(round(estimates$quartile,sig.digit)), fill=TRUE)
  cat("And the deciles are", target.value.formatter(round(estimates$decile,sig.digit)), fill=TRUE)
}

##' Calculate the (first) peak week in a vector of weekly observations
##'
##' @param trajectory a vector of weekly observations
##' @param ... ignored
##' @return integer vector, typically of length 1; indices of elements that are
##'   \code{>=} all other elements
##'
##' @export
pwk = function(trajectory, ...){
  if (any(is.na(trajectory))) {
    stop ("NA elements are not allowed when calculating =pwk=.")
  }
  return (which(trajectory==max(trajectory)))
}

##' Calculate the peak height of a trajectory
##'
##' @param trajectory a vector
##' @param ... ignored
##' @return a scalar --- the maximum value of the vector
##'
##' @export
pht = function(trajectory, ...){
  return (max(trajectory))
}

##' Calculate the onset of a trajectory
##'
##' @param trajectory a vector
##' @param baseline the onset threshold
##' @param is.inseason logical vector length-compatible with \code{trajectory},
##'   acting as a mask on possible output values
##' @param ... ignored
##' @return the first index which is part of the in-season and is part of a run
##'   of consecutive observations above the onset threshold that lasts at least
##'   for two additional indices
##'
##' @export
ons = function(trajectory, baseline, is.inseason, ...) {
  above.baseline = trajectory >= baseline
  next.three.above.baseline =
      above.baseline &
      dplyr::lead(above.baseline, 1L) &
      dplyr::lead(above.baseline, 2L)
  return (which(is.inseason & next.three.above.baseline)[1L][[1L]])
}

##' Calculate the onset of a trajectory
##'
##' @param trajectory a vector
##' @param baseline the onset threshold
##' @param is.inseason logical vector length-compatible with \code{trajectory},
##'   acting as a mask on possible output values
##' @param ... ignored
##' @return a single non-\code{NA} integer: the number of indices which are part
##'   of the in-season and part of a run of at least three consecutive
##'   observations above the onset threshold
##'
##' @export
dur = function(trajectory, baseline, is.inseason, ...) {
  above.baseline = trajectory >= baseline
  next.three.above.baseline =
    above.baseline &
    dplyr::lead(above.baseline, 1L) &
    dplyr::lead(above.baseline, 2L)
  part.of.run =
    (next.three.above.baseline |
     dplyr::lag(next.three.above.baseline, 1L) |
     dplyr::lag(next.three.above.baseline, 2L))
  return (sum(is.inseason & part.of.run, na.rm=TRUE))
}

## ##' constructor should take in the bare minimum to create forecasts
## ##' i.e. full.dat, n.sims, control.list, baseline, etc.
## sim <- function(...) {
##   structure(list(...), class = "sim")
## }

c_for_named_lists = function(sim, ..., recursive=FALSE, use.names=TRUE) {
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
  if (!all(
         rlang::names2(dots) != "" |
         sapply(dots, function(arg) {
           is.list(arg) && all(rlang::names2(arg) != "")
         })
       )) {
    stop ('All components to add to the sim object must be (a) named (and/)or (b) lists with every component named.  Names of "" are treated as invalid.')
  }
  ## forward operation to c method used for lists:
  input.as.list = unclass(sim)
  stopifnot(!anyDuplicated(names(input.as.list)))
  result.as.list = c(input.as.list, ..., recursive=FALSE, use.names=TRUE)
  ## check for duplicate names:
  duplicate.name.flags = duplicated(names(result.as.list))
  if (any(duplicate.name.flags)) {
    stop (paste0("All components must be uniquely named; adding the specified components would result in the following duplicates:\n",
                 paste0(sapply(names(result.as.list)[duplicate.name.flags], as.name), collapse="\n")))
  }
  ## re-class list result into sim object:
  result = structure(result.as.list, class="sim")
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

##' Add information to a \code{target_forcast} object
##'
##' Uses \code{recursive=FALSE} and \code{use.names=TRUE} when forwarding to the
##' list \code{c} method; any attempt to override these values will generate an
##' error.
##'
##' @param target.forcast a \code{target_forcast} object
##' @param ... list of components to add to the \code{target_forcast} object;
##'   must lead to a resulting \code{target_forcast} object with components that
##'   are all uniquely, nontrivially (\code{!=""}) named
##' @return \code{target_forcast} object with the given components appended
##'
##' @method c target_forcast
##' @export
##' @export c.target_forcast
c.target_forcast = function(target.forcast, ...) {
  return (c_for_named_lists(target.forcast, ...))
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
