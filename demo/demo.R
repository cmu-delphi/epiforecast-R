library(epiforecast)

## Make dummy csv
area.name = "hhs1"
hhs1.dat = fetchEpidataFullDat("fluview", area.name, "wili",
                               min.points.in.season=52L,
                               first.week.of.season = 21L,
                               cache.file=sprintf("fluview_%s_fetch.Rdata", area.name))

## Create a csv in the correct, desirable format
filename="./correct.csv"
hhs1.dat = lapply(hhs1.dat, function(myvec){ if(length(myvec)<53) return(c(myvec,NA))  else {return(myvec)}})
correct.hhs1.dat = do.call(cbind, hhs1.dat)
correct.hhs1.dat = correct.hhs1.dat[-53,]
write.csv(correct.hhs1.dat, file = filename,row.names=FALSE)
full.dat = read.from.file(filename)


## Try EB with basic options
mysim = eb.sim(full.dat, n.sims=1000)
epiforecast:::plot.sim(mysim)
targets = epiforecast:::forecast.sim(mysim,'pht')


## Try EB with more simulation settings
control.list = get_eb_control_list(sd.option="prior",max.match.length=5)
mysim = eb.sim(full.dat, n.sims=100, control.list=control.list)
epiforecast:::plot.sim(mysim)


## Try BR
mysim = br.sim(full.dat, n.sims=100, bootstrap=T)
epiforecast:::plot.sim(mysim)
epiforecast:::print.sim(mysim,verbose=TRUE)
targets = epiforecast:::forecast.sim(mysim, "pht")

## Try BR with more simulation settings
control.list = get_br_control_list(df=5, cv.rule="1se")
mysim = br.sim(full.dat, n.sims=100, bootstrap=T, control.list=control.list)
epiforecast:::plot.sim(mysim)
epiforecast:::print.sim(mysim,verbose=TRUE)
targets = epiforecast:::forecast.sim(mysim, "pht", plot.hist=TRUE)
targets = epiforecast:::forecast.sim(mysim, "pwk", plot.hist=TRUE)



## Try twkde
mysim = twkde.sim(full.dat)
epiforecast:::plot.sim(mysim)
