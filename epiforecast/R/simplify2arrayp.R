## These functions are modified R methods from base R, under GPL 2+. Original
## license information from license() is below. The modification is from Aaron
## Rumack.

## This software is distributed under the terms of the GNU General
## Public License, either Version 2, June 1991 or Version 3, June 2007.
## The terms of version 2 of the license are in a file called COPYING
## which you should have received with
## this software and which can be displayed by RShowDoc("COPYING").
## Version 3 of the license can be displayed by RShowDoc("GPL-3").

## Copies of both versions 2 and 3 of the license can be found
## at https://www.R-project.org/Licenses/.

## A small number of files (the API header files listed in
##                          R_DOC_DIR/COPYRIGHTS) are distributed under the
## LESSER GNU GENERAL PUBLIC LICENSE, version 2.1 or later.
## This can be displayed by RShowDoc("LGPL-2.1"),
## or obtained at the URI given.
## Version 3 of the license can be displayed by RShowDoc("LGPL-3").

## 'Share and Enjoy.'

simplify2arrayp = function (x, higher = TRUE)
{
  if (length(common.len <- unique(as.vector(lengths(x)))) > 1L)
    return(x)
  if (common.len == 1L)
    unlist(x, recursive = FALSE)
  else if (common.len > 1L) {
    n <- length(x)
    r <- unlist(x, recursive = FALSE, use.names = FALSE)
    if (higher && length(c.dim <- unique(lapply(x, dim))) ==
        1 && is.numeric(c.dim <- c.dim[[1L]]) && prod(d <- c(c.dim,
                                                             n)) == length(r)) {
      iN1 <- is.null(n1 <- dimnames(x[[1L]]))
      n2 <- names(x)
      dnam <- if (!(iN1 && is.null(n2)))
                c(if (iN1) rep.int(list(n1), length(c.dim)) else n1,
                  list(n2))
      array(r, dim = d, dimnames = dnam)
    }
    else if (prod(d <- c(common.len, n)) == length(r))
      array(r, dim = d, dimnames = if (!(is.null(n1 <- names(x[[1L]])) &
                                         is.null(n2 <- names(x))))
                                     list(n1, n2))
    else x
  }
  else x
}
