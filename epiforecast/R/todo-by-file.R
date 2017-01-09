######################
## interface.R #######
######################

## todo length checks --- time.of.forecast not only with dat, but with fit
## todo remove n.out arg --- determine by min fit length

## ## source("../sample/sample_config_flu_1516.R")
## source("loaders.R")
## source("fitters.R")
## source("nidss/nidss_fetch_data.R", chdir=TRUE)

## nidss.datfit = fetchNIDSSDatFit("flu", "nationwide")

## olddatplus = nidss.datfit$olddatplus
## oldfit = nidss.datfit$fit
## include.ss = nidss.datfit$include.ss
## current.s = nidss.datfit$current.s
## first.year = nidss.datfit$first.year
## first.model.week = nidss.datfit$first.model.week

## fit.ss = include.ss[include.ss < current.s]
## exclude.2009.pandemic.season=TRUE
## if (exclude.2009.pandemic.season) {
##     fit.ss <- fit.ss[-2]
##     oldfit <- list(f=oldfit$f[,-2], tau=oldfit$tau[-2])
## }
## train.ss = fit.ss
## test.s = max(train.ss)+1

## newdat = olddatplus.to.newdat(olddatplus)
## newdat.attributes = attributes(newdat)
## newdat <- newdat[match(fit.ss, include.ss)]
## newdat.attributes$names <- names(newdat)
## attributes(newdat) <- newdat.attributes

## qwer = fit.eb.control.list(oldfit.to.newfit(oldfit), get.eb.control.list())
## asdf = eb.createForecasts(newdat, olddatplus$wili[olddatplus$season==test.s], oldfit.to.newfit(oldfit), 0)
## asdf = eb.createForecasts(newdat, olddatplus$wili[olddatplus$season==test.s], oldfit.to.newfit(oldfit), 1L)
## source("plotters.R")
## newfit = smooth.curves.to.newfit(eb.fitSmoothCurves(newdat))
## matplot.newdat(newdat)
## matplot.newfit(newdat, newfit)
## seriesplot.newfit(newdat, smooth.curves.to.newfit(eb.fitSmoothCurves(newdat)))

## xxx instead of n.out, allow NA's in the future trajectories, just fill in all; use !is.na as another ii.match mask?
## todo explicitly make object that represents a distribution of curves, corresponding fitting functions, then the conditioning method?
## todo rename forecast time to something with "ind"?
## todo documentation
## todo imports
## todo examples
