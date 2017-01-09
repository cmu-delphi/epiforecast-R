context("Testing the read.from.file() function..")

## ## Setup
## dummytable = matrix(c(runif(n=9,1,3),NA),ncol=5)
## colnames(dummytable) = c(1997:2001)


## ## Tests:
## test_that("Return type is list of vector", {
##     write.table(dummytable,
##                 file = "a.csv",
##                 col.names = T,
##                 row.names=F,
##                 sep=",")
##     mydat  = read.from.file("a.csv")
##     mydat
##     expect_equal("list", typeof(mydat))
##     file.remove("a.csv")
## })


## Fetch some data
area.name = "hhs1"
hhs1.dat = fetchEpidataFullDat("fluview", area.name, "wili",
                               min.points.in.season=52L,
                               first.week.of.season = 21L,
                               cache.file=sprintf("fluview_%s_fetch.Rdata", area.name))
alt.names1=paste("season",1:length(names(hhs1.dat)))
alt.names2 = alt.names1;    alt.names2[5] = "" 
alt.names3 = rep("", length(names(hhs1.dat)))

## When one label is missing
    

## Create csv with /no/ names
test_that("When alls labels are missing, error is thrown", {
    filename="./all.names.missing.file.csv"
    all.names.missing.hhs1.dat = hhs1.dat
    names(all.names.missing.hhs1.dat) = alt.names3
    all.names.missing.hhs1.dat = do.call(cbind, all.names.missing.hhs1.dat)
    all.names.missing.hhs1.dat = all.names.missing.hhs1.dat[-53,]
    write.csv(all.names.missing.hhs1.dat, file = filename, row.names=FALSE)
    expect_error(full.dat = read.from.file(filename))
})

## Create csv with one name missing
test_that("When one label is missing, error is thrown", {
    filename="./one.name.missing.file.csv"
    one.name.missing.hhs1.dat = hhs1.dat
    names(one.name.missing.hhs1.dat) = alt.names2
    one.name.missing.hhs1.dat = do.call(cbind, one.name.missing.hhs1.dat)
    one.name.missing.hhs1.dat = one.name.missing.hhs1.dat[-53,]
    write.csv(one.name.missing.hhs1.dat, file = filename, row.names=FALSE)
    expect_error(full.dat = read.from.file(filename))
})

## Last season is not partially observed
test_that("When last season is not partially observed, error is thrown", {
    filename="./last.column.weird.file.csv"
    last.column.weird.hhs1.dat = do.call(cbind, hhs1.dat)
    write.csv(last.column.weird.hhs1.dat, file = filename,row.names=FALSE)
    expect_error(full.dat = read.from.file(filename))
})

## CSV is correct, so everything should be fine here.
test_that("When last season is not partially observed, error is thrown", {
    filename="./correct.csv"
    correct.hhs1.dat = do.call(cbind, hhs1.dat)
    correct.hhs1.dat = correct.hhs1.dat[-53,]
    write.csv(correct.hhs1.dat, file = filename,row.names=FALSE)
    full.dat = read.from.file(filename)
    expect_equal(class(full.dat), "list")
}
