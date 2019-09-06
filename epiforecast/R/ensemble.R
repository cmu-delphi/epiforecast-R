##' @include simclass.R
NULL

##' @export
degenerate_em_weights = function(distr.cond.lkhds,
                                 init.weights=rep(1/dim(distr.cond.lkhds)[[2L]],dim(distr.cond.lkhds)[[2L]]),
                                 instance.weights=rep(1, dim(distr.cond.lkhds)[[1L]]),
                                 stop.eps = sqrt(.Machine[["double.eps"]])) {
  if (any(init.weights < 1e-10)) {
    stop("All init.weight's must be >=1e-10")
  }
  if (!isTRUE(all.equal(1, sum(init.weights)))) {
    stop("Sum of init.weight's must be all.equal to 1.")
  }

  ## Set some constants:
  n.obs = dim(distr.cond.lkhds)[[1L]]
  n.distr = dim(distr.cond.lkhds)[[2L]]
  t.distr.cond.lkhds = t(distr.cond.lkhds) # dim: n.distr x n.obs

  ## Set initial values of variables adjusted each step:
  weights = init.weights # length: n.distr
  t.lkhds = init.weights*t.distr.cond.lkhds # dim: n.distr x n.obs
  marginals = colSums(t.lkhds) # length: n.obs
  log.lkhd = weighted.mean(log(marginals), instance.weights) # scalar
  ## log.lkhds = list(log.lkhd)
  if (log.lkhd == -Inf) {
    stop ("All methods assigned a probability of 0 to at least one observed event.")
  } else {
    repeat {
      old.log.lkhd = log.lkhd # scalar
      weights <- matrixStats::colWeightedMeans(t(t.lkhds)/marginals, instance.weights) # length: n.distr
      t.lkhds <- weights*t.distr.cond.lkhds # dim: n.distr x n.obs
      marginals <- colSums(t.lkhds) # length: n.obs
      log.lkhd <- weighted.mean(log(marginals), instance.weights) # scalar
      ## xxx inefficient
      ## log.lkhds <- c(log.lkhds,list(log.lkhd))
      stopifnot (log.lkhd >= old.log.lkhd)
      if (log.lkhd-old.log.lkhd <= stop.eps || (log.lkhd-old.log.lkhd)/-log.lkhd <= stop.eps) {
        break
      }
    }
  }
  return (weights)
}
## Test for agreement:
##
## ## |cvMixtureCoeffLogLkhdTestValues|: choose the best mixing
## ## coefficient between two distributions from a set of test
## ## coefficients to minimize CV-estimated log-likelihood loss given the
## ## observed indicator values.
## ## |indicators|: nxm matrix, observed indicator values
## ## |distr1|: nxm matrix, first distr's probs/E[indicator values]
## ## |distr2|: nxm matrix, second distr's probs/E[indicator values]
## ## n: number of events
## ## m: number of folds
## ## |test.coeffs|: different mixing coefficients of distr1
## ## |safety.lambda|: coefficient with which to pre-mix a uniform distr in with distr1, distr2 to avoid 0 entries
## ## result: best test coefficient of distr1
## mixture_coef_loglkhd_test_values = function(indicators, distr1, distr2,
##     test.coeffs=c(0:1000/1000,10^(0:-120/10),1-10^(0:-120/10)),
##     safety.lambda=1e-15) {
##     n = nrow(indicators) # number of events per fold
##     m = ncol(indicators) # number of folds
##     distr1 <- (1-safety.lambda)*distr1 + safety.lambda/n # premix with tiny uniform
##     distr2 <- (1-safety.lambda)*distr2 + safety.lambda/n # premix with tiny uniform
##     elosses = sapply(test.coeffs, function(test.coeff)
##         sum(-indicators*log(test.coeff*distr1+(1-test.coeff)*distr2))
##                     ) # cv avg loss for each test coefficient
##     return (test.coeffs[which.min(elosses)])
## }
## {
##   print(mixture_coef_loglkhd_test_values(matrix(c(1,0,1,0,1,0,1,0),2),matrix(c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8),2),matrix(c(1,0,1,0,1,0,1,0),2)))
##   print(degenerate_em_weights(distrsToLkhds(matrix(c(1,0,1,0,1,0,1,0),2),matrix(c(0.2,0.8,0.2,0.8,0.2,0.8,0.2,0.8),2),matrix(c(1,0,1,0,1,0,1,0),2))))
##   print(degenerate_em_weights(matrix(c(0.2,0.2,0.2,0.2, 1,1,1,1),4)))
## }
##
## ## Test using instance.weights:
## degenerate_em_weights(matrix(c(0.5,0.2,0.2,0.2, 0.8,0.4,0.1,0.1),,2L))
## degenerate_em_weights(matrix(c(0.5,0.2,0.2, 0.8,0.4,0.1),,2L),instance.weights=c(1,1,2))
## ## Scale of instance.weights does not impact the result:
## degenerate_em_weights(matrix(c(0.5,0.2,0.2, 0.8,0.4,0.1),,2L),instance.weights=c(1,1,2)/5)

## Illustrative example/check:
## {
##     m = 1
##     a = -3
##     b = 4
##     par(ask=TRUE)
##     for (i in 1:10) {
##         hist(rnorm(800),col=rgb(1,0,0,0.5,1),freq=FALSE)
##         hist(rnorm(800,m),col=rgb(0,1,0,0.5,1),freq=FALSE,add=TRUE)
##         hist(runif(800,a,b),col=rgb(0,0,1,0.5,1),freq=FALSE,add=TRUE)
##     }
##     print('distr1 coeff:')
##     f = function(points) empiricalBucketMasses(points, (-5:5)*0.6, tack.neg.inf=TRUE, tack.pos.inf=TRUE)
##     cvMixtureCoeffLogLkhd(replicate(800,f(rnorm(1))),replicate(800,f(rnorm(5000,m))),replicate(800,f(runif(5000,a,b))))
## }

##' @export
lasso_lad_coef = function(y, X, include.intercept=TRUE) {
  y <- as.vector(y)
  coef(quantreg::rq(if(include.intercept) obs ~ . else obs ~ . + 0,
                    method="lasso",
                    data=cbind(obs=y, as.data.frame(X)),
                    lambda=mean(abs(y-X))))
}

##' @importFrom Matrix rBind cBind Diagonal Matrix
##' @export
simplex_lad_weights = function(y, X) {
  if (length(y) != nrow(X)) stop("length(y) != nrow(X)")
  n = nrow(X)
  p = ncol(X)

  ## p beta, n s+, n s-
  objective.in = c(rep(0,p),rep(1,n),rep(1,n))
  const.mat = rBind(cBind(                     X,   Matrix::Diagonal(n),  -Matrix::Diagonal(n) ),
                    cBind( Matrix::Matrix(0,n,p),   Matrix::Diagonal(n), Matrix::Matrix(0,n,n) ),
                    cBind( Matrix::Matrix(0,n,p), Matrix::Matrix(0,n,n),   Matrix::Diagonal(n) ),
                    cBind(   Matrix::Diagonal(p), Matrix::Matrix(0,p,n), Matrix::Matrix(0,p,n) ),
                        c(              rep(1,p),              rep(0,n),              rep(0,n) ))
  const.dir = c(rep("=" , n),
                rep(">=", n),
                rep(">=", n),
                rep(">=", p),
                    "="     )
  const.rhs = c(       y,
                rep(0,n),
                rep(0,n),
                rep(0,p),
                       1)

  structure(
    lpSolve::lp("min",objective.in,,const.dir,const.rhs,dense.const=as.matrix(summary(const.mat)))$solution[1:p],
    names=colnames(X)
  )
}

uniform_forecast = function(...) {
  return (structure(list(), class="uniform_forecast"))
}

##' @method target_forecast uniform_forecast
##' @export
##' @export target_forecast.uniform_forecast
target_forecast.uniform_forecast = function(uniform.forecast, ..., target.name, target.spec) {
  bin.info = target.spec[["bin_info_for"]](...)
  breaks = bin.info[["breaks"]]
  non.na.bin.representatives = breaks[-1L] - 0.5*diff(breaks)
  target.values = c(
    non.na.bin.representatives,
    get_na_value_or_empty_for_target(target.spec, ...)
  )
  target.weights = rep(1/length(target.values), length(target.values))
  return (structure(list(target.values = stats::setNames(list(target.values), target.name),
                         target.weights = target.weights,
                         method.settings = list(uniform.pseudoweight.total=0,smooth.sim.targets=FALSE)
                         ),
                    class="target_forecast"))
}

##' @method print uniform_forecast
##' @export
##' @export print.uniform_forecast
print.uniform_forecast = function(x, ...) {
  cat("uniform_forecast: produces uniform-probability forecasts for any target", fill=TRUE)
}

##' Add information to a \code{uniform_forecast} object
##'
##' Uses \code{recursive=FALSE} and \code{use.names=TRUE} when forwarding to the
##' list \code{c} method; any attempt to override these values will generate an
##' error.
##'
##' @param uniform.forecast a \code{uniform_forecast} object
##' @param ... list of components to add to the \code{uniform_forecast} object; must lead to
##'   a resulting \code{uniform_forecast} object with components that are all uniquely,
##'   nontrivially (\code{!=""}) named
##' @return \code{uniform_forecast} object with the given components appended
##'
##' @method c uniform_forecast
##'
##' @export
##' @export c.uniform_forecast
c.uniform_forecast = c_for_named_lists

