context("Testing the check.file.contents() function..")

## Setup
dummytable = matrix(c(runif(n=9,1,3),NA),ncol=5)
dummytable.full = matrix(runif(n=10,1,3),ncol=5)
colnames(dummytable.full) = c(1997:2001)
library(s

## Tests:

test_that("Missing column names returns error.", {
    wrong.names = c("a","b","c","","d")
    colnames(dummytable) = wrong.names
    write.table(dummytable,
                file = "a.csv",
                col.names = T,
                row.names=F,
                sep=",")
    expect_error(check.file.contents("a.csv"))
    file.remove("a.csv")
})


test_that("Column names containing NA returns error.", {
    wrong.names = c(NA,1:4)
    colnames(dummytable) = wrong.names
    write.table(dummytable,
                file = "a.csv",
                col.names = T,
                row.names= F,
                sep=",")
    expect_error(check.file.contents("a.csv"))
    file.remove("a.csv")
})



test_that("Last column being full returns error.", {
    write.table(dummytable.full,
                file = "a.csv",
                col.names = T,
                row.names=F,
                sep=",")
    print(fread("a.csv",header=TRUE))
    expect_error(check.file.contents("a.csv"))
    file.remove("a.csv")
})
