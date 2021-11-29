dquantile = function(dist, qs, bins) {
  d = c(0,cumsum(dist))
  idxs = apply(outer(qs,d,FUN="<"),1,function(x) { which(x)[1]})
  idxs[is.na(idxs)] = length(d)
  p = (qs - d[idxs-1]) / dist[idxs-1]
  result = bins[idxs-1] + p*(bins[idxs]-bins[idxs-1])
  return(result)
}

qdistribution = function(qs, qvals, bins) {
  result = (2:length(bins)) %>>%
    lapply(function(bi) {
      idxs = which(bins[bi-1] <= qs & qs < bins[bi])
      if (!(1 %in% idxs)) {
        left.wt = ((qs[idxs[1]] - bins[bi-1]) / (qs[idxs[1]] - qs[idxs[1]-1])) * (qvals[idxs[1]] - qvals[idxs[1]-1])
      } else {
        left.wt = 0
      }
      if (!(length(qs) %in% idxs)) {
        right.wt = (bins[bi] - qs[idxs[length(idxs)]]) / (qs[idxs[length(idxs)]+1] - qs[idxs[length(idxs)]]) * (qvals[idxs[length(idxs)]+1] - qvals[idxs[length(idxs)]])
      } else {
        right.wt = 0
      }
      if (length(idxs) > 0) {
        return(left.wt + right.wt + qvals[idxs[length(idxs)]] - qvals[idxs[1]])
      } else {
        idx = which(bins[bi-1] < qs)[1]
        wt = (bins[bi]-bins[bi-1]) / (qs[idx]-qs[idx-1]) * (qvals[idx]-qvals[idx-1])
        return(wt)
      }
    })
  mode(result) = "numeric"
  return(result)
}

get_observed_quantile = function(forecast, val, bins) {
  low.bins = which(bins <= val)
  right.idx = which(bins > val)[1]
  right.wt = forecast[right.idx-1]*(bins[right.idx]-val)/(bins[right.idx]-bins[right.idx-1])
  return(sum(forecast[low.bins]) - right.wt)
}

get_observed_quantile_cdc = function(forecast, val) {
  ##low.bins = which(bins <= val)
  ##right.idx = which(bins > val)[1]
  correct.bin = sum(1:length(val) * val)
  if (correct.bin > 1) {
    low.bins = forecast[1:(correct.bin-1)]
  }  else {
    low.bins = c()
  }
  ##if (right.idx == length(bins)) {
  ##  right.wt = forecast[right.idx-1]*(bins[right.idx]-val)/(bins[right.idx]-bins[right.idx-1])
  ##} else {
  ##  right.wt = 0.5*forecast[right.idx-1]
  ##}
  return(sum(low.bins) + 0.5*forecast[[correct.bin]])
}

quantile_bias = function(qs) { 
  sorted.qs = sort(qs)
  freq = seq(0,1,length.out=length(qs))
  return(sum(freq-sorted.qs)/length(qs))
}

quantile_spread = function(qs) {
  sorted.qs = sort(qs)
  freq = seq(0,1,length.out=length(qs))
  diffs = freq-sorted.qs
  return((sum(diffs[freq<=0.5]) - sum(diffs[freq>0.5]))/length(qs))
}

quantile_rmse = function(qs) {
  sorted.qs = sort(qs)
  freq = seq(0,1,length.out=length(qs))
  return(sqrt(sum((freq-sorted.qs)^2)/length(qs)))
}

forecast_mean = function(forecast, bins) {
  return(weighted.mean(bins,forecast))
}

forecast_std = function(forecast, bins) {
  m = forecast_mean(forecast, bins)
  return(sqrt(sum(forecast * (bins-m)^2)))
}

calibrate_forecast_null = function(forecast, qqs, bins) {
  forecast
}

calibrate_forecast = function(forecast, qqs, bins=-1, alpha=0, fit_spline=FALSE) {
  if (bins == -1) {
    bins = 1:(length(forecast)+1)
  }
  qqs = c(qqs, seq(0,1,length.out=(alpha/(1-alpha))*length(qqs)))
  old.quantiles = seq(0,1,length.out=1001)
  if (fit_spline) {
    sf = splinefun(seq(0,1,length.out=length(qqs)), sort(qqs), method="hyman")
    new.quantiles = sf(old.quantiles)
  } else {
    new.quantiles = quantile(qqs,old.quantiles)
  }
  new.quantiles = pmin(pmax(0,new.quantiles),1)
  new.quantiles[1] = 0
  new.quantiles[length(new.quantiles)] = 1
  forecast.quantiles = dquantile(forecast,new.quantiles,bins)
  forecast.quantiles[1] = bins[[1]]
  forecast.quantiles[length(forecast.quantiles)] = bins[[length(bins)]]
  forecast.quantiles = forecast.quantiles[order(forecast.quantiles)] # Sometimes precision problems
  calibrated.forecast = qdistribution(forecast.quantiles,old.quantiles,bins)
  return(calibrated.forecast)
}

calibrate_forecast_pseudocount = function(k) {
  function(forecast, qqs, bins) {
    qqs = c(qqs,seq(0,1,length.out=k))
    return(calibrate_forecast(forecast,qqs,bins))
  }
}

calibrate_forecast_beta = function(forecast, qqs, bins=-1, alpha=0) {
  if (bins == -1) {
    bins = 1:(length(forecast)+1)
  }
  qqs = pmin(pmax(qqs,0.001),0.999)
  f = function(ab) {
    sum(log(qqs))*(ab[[1]]-1) + sum(log(1-qqs))*(ab[[2]]-1) - length(qqs)*lbeta(ab[[1]],ab[[2]]) }
  ab = stats::optim(c(1,1),f,control=list(fnscale=-1))[["par"]]
  calibrated.forecast = diff((1-alpha)*pbeta(cumsum(c(0,forecast)),ab[[1]],ab[[2]])) + alpha*forecast
  return(calibrated.forecast)
}

calibrate_forecast_uniform_smoothing = function(cal, alpha, fit_spline=FALSE) {
  function(forecast, qqs, bins=-1) {
    return(cal(forecast, qqs, bins=bins, alpha=alpha, fit_spline=fit_spline))
  }
}

