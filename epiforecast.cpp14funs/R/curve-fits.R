
smooth.curves.to.fit = function(smooth.curves, type="Gaussian") {
    return (lapply(seq_along(smooth.curves$sigma.hat), function(fit.s.i) {
        list(f=smooth.curves$smooth.obj[[fit.s.i]], tau=smooth.curves$sigma.hat[fit.s.i], type=type)
    }))
}

fit.to.oldfit = function(fit) {
    f = lapply(fit, `[[`, "f")
    f <- sapply(f, `[`, seq_len(min(lengths(f))))
    tau = sapply(fit, `[[`, "tau")
    return (list(f=f, tau=tau))
}

oldfit.to.fit = function(oldfit, type="Gaussian") {
    return (lapply(seq_along(oldfit$tau), function(fit.s.i) {
        list(f=oldfit$f[,fit.s.i], tau=oldfit$tau[fit.s.i], type=type)
    }))
}
