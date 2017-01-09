---
title: "Fetching / Loading Data"
author: "Justin"
date: "2016-07-17"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Fetch data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

First, setup by loading the package and setting working directory


```
## Loading epiforecast
```

Here is an example.


```r
## helper function to merge gft with flu
mymerge = function(df1, df2, merge.names = "", merge.by = ""){#colnames.to.merge
    merge.dat = merge(df1,df2,by=c("year","week","date"),all=T)
    merge.dat[,"week"] = as.numeric(merge.dat[,"week"])
    merge.dat[,"year"] = as.numeric(merge.dat[,"year"])
    ord.by.dat = order(merge.dat[,c("year")],merge.dat[,c("week")])
    merge.dat = merge.dat[ord.by.dat,]
    merge.dat = merge.dat[,c("season.y","year","week","date","gft.num","wili")]
    colnames(merge.dat)[colnames(merge.dat)=="season.y"] = "season"
    return(merge.dat)
}


## National
filename = file.path(outputdir, "fluview_nat_df.csv")
fluview.nat.all.df =   # If I need just 2003 and on, I do #fluview.nat.recent.df()
    trimPartialPastSeasons(fetchEpidataDF("fluview", "nat",
                                          first.week.of.season=21L,
                                          cache.file="fluview_nat_allfetch.Rdata"),
                           "wili",
                           min.points.in.season=33L)

fluview.nat.all.df = fluview.nat.all.df[,c("region","wili","year","week","date","season","model.week")]
gft.nat.df = fetchEpidataDF("gft","nat")
gft.nat.df = gft.nat.df[,c("location","year","week","model.week","date","season","num")]
colnames(gft.nat.df)[which(colnames(gft.nat.df) == "num")] = "gft.num"
nat.dat = mymerge(df1=gft.nat.df, df2=fluview.nat.all.df)
cat("writing national data",fill=T)
write.csv(nat.dat, file=filename,row.names=F)


## By region
regs = paste0("reg",1:10)
hhss = paste0("hhs",1:10)
for(jj in 1:10){
    cat("writing region",jj, fill=T)
    myreg = regs[jj]
    myhhs = hhss[jj]
    cachename = paste0("fluview_",myhhs,"_fetch.Rdata")
    filename =  file.path(outputdir, paste0("fluview_",myhhs,"_df.csv"))
    regdat = trimPartialPastSeasons(fetchEpidataDF("fluview", myhhs,
                                                   first.week.of.season=21L),
                                                   ## cache.file=cachename),
                                    "wili", min.points.in.season=33L)
    regdat2 = regdat[,c("region","wili","year","week","date","season","model.week")]
    gft.reg.dat = fetchEpidataDF("gft","hhs1")
    colnames(gft.reg.dat)[which(colnames(gft.reg.dat) == "num")] = "gft.num"
    reg.dat = mymerge( df1=gft.reg.dat, df2=regdat2)

    # write the file
    write.csv(reg.dat, file = filename, row.names=F)
}
```
