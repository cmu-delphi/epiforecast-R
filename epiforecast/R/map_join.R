## author_header begin
## Copyright (C) 2017 Logan C. Brooks
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

##' @include namesp.R
NULL

##' @export
array_proxy = function(prefix, original.dimnamesp, indices.into.original) {
  ## xxx avoid storing indices.into.original?
  ## original.named.dimp = sapply(original.dimnamesp, length)
  ## indices.into.original = array(seq_len(prod(original.named.dimp)), original.named.dimp)
  dimnames(indices.into.original) <- original.dimnamesp
  structure(
    list(prefix=prefix,
         original.dimnamesp=original.dimnamesp,
         indices.into.original=indices.into.original),
    class="array_proxy"
  )
}

##' @method print array_proxy
##' @export
##' @export print.array_proxy
print.array_proxy = function(x,...) {
  x.impl = unclass(x)
  original.dimnamesp = x.impl[["original.dimnamesp"]]
  current.dimnamesp = dimnamesp(x.impl[["indices.into.original"]])
  cat("array_proxy with filepath pattern ")
  cat(paste0(x.impl[["prefix"]],".",paste(names(x.impl[["original.dimnamesp"]]),collapse=".")))
  cat("\nOriginal dimnamesp:\n")
  print(original.dimnamesp)
  cat("Current dimnamesp:")
  if (identical(original.dimnamesp, current.dimnamesp)) {
    cat(" same as original\n")
  } else {
    cat("\n")
    print(current.dimnamesp)
  }
  cat("Current length: ")
  cat(length(x.impl[["indices.into.original"]]))
  cat("\nFirst few (up to six) indices into original:\n")
  print(x.impl[["indices.into.original"]][seq_len(min(6L, length(x.impl[["indices.into.original"]])))])
}

##' @method dim array_proxy
##' @export
##' @export dim.array_proxy
dim.array_proxy = function(x) {
  x.impl = unclass(x)
  return (dim(x.impl[["indices.into.original"]]))
}

##' @method dim<- array_proxy
##' @export
##' @export dim<-.array_proxy
`dim<-.array_proxy` = function(x, value) {
  x.impl = unclass(x)
  dim(x.impl[["indices.into.original"]]) <- value
  structure(x.impl, class="array_proxy")
}

##' @method dimnames array_proxy
##' @export
##' @export dimnames.array_proxy
dimnames.array_proxy = function(x) {
  x.impl = unclass(x)
  return (dimnames(x.impl[["indices.into.original"]]))
}

##' @method dimnames<- array_proxy
##' @export
##' @export dimnames<-.array_proxy
`dimnames<-.array_proxy` = function(x, value) {
  x.impl = unclass(x)
  dimnames(x.impl[["indices.into.original"]]) <- value
  structure(x.impl, class="array_proxy")
}

##' @method [ array_proxy
##' @export
##' @export [.array_proxy
`[.array_proxy` = function(x, ...) {
  result = unclass(x)
  result[["indices.into.original"]] <- result[["indices.into.original"]][...]
  class(result) <- "array_proxy"
  return (result)
}

##' @method [<- array_proxy
##' @export
##' @export [<-.array_proxy
`[<-.array_proxy` = function(x, ..., value) {
  stop ("Assignment to an array_proxy is forbidden.")
}

##' @method [[ array_proxy
##' @export
##' @export [[.array_proxy
`[[.array_proxy` = function(x, ...) {
  x.impl = unclass(x)
  index = x.impl[["indices.into.original"]][[...]]
  array.ind = arrayInd(index, sapply(x.impl[["original.dimnamesp"]], length))
  ## xxx pre-compute and save original.named.dimp?
  array.char.ind = mapply(`[[`, x.impl[["original.dimnamesp"]], array.ind)
  array.char.ind <- gsub("/","-",array.char.ind)
  filepath = paste0(x.impl[["prefix"]],".",paste(array.char.ind, collapse="."),".rds")
  result = readRDS(filepath)
  return (result)
}

##' @method [[<- array_proxy
##' @export
##' @export [[<-.array_proxy
`[[<-.array_proxy` = function(x, ..., value) {
  stop ("Assignment to an array_proxy is forbidden.")
}

##' @method as.list array_proxy
##' @export
##' @export as.list.array_proxy
as.list.array_proxy = function(X, FUN, ...) {
  stop ("as.list.array_proxy is disabled; storing results in a list may result in high memory usage.")
}

##' @method as.vector array_proxy
##' @export
##' @export as.vector.array_proxy
as.vector.array_proxy = function(X, FUN, ...) {
  stop ("as.vector.array_proxy is not implemented; storing results in a vector may result in high memory usage.")
}

##' @method length array_proxy
##' @export
##' @export length.array_proxy
length.array_proxy = function(x) {
  x.impl = unclass(x)
  length(x.impl[["indices.into.original"]])
}

##' @method names array_proxy
##' @export
##' @export names.array_proxy
names.array_proxy = function(x) {
  x.impl = unclass(x)
  names(x.impl[["indices.into.original"]])
}

##' @method aperm array_proxy
##' @export
##' @export aperm.array_proxy
aperm.array_proxy = function(a, perm, ...) {
  result = unclass(a)
  result[["indices.into.original"]] <- aperm(result[["indices.into.original"]], perm, ...)
  class(result) <- "array_proxy"
  return (result)
}

##' Map a function over the natural (Cartesian/other) join of array-like objects
##'
##' @param f the function to map
##' @param arraylike.args a list of array-like objects with named dimnames
##'   (including "scalars" with \code{\link{ndimp}} of 0) and/or
##'   \code{\link{no_join}} objects
##' @return an object of \code{class} \code{"array"} and \code{mode}
##'   \code{"list"} containing outputs from \code{f}, or a single output from
##'   \code{f} if all \code{arraylike.args} are scalars or \code{no_join}'s
##' @details The function \code{f} will be called with a number of arguments
##'   equal to the number of array-like arguments, with the argument names
##'   specified in the \code{arraylike.args} list; each argument will correspond
##'   to an element of one of the \code{arraylike.args}.
##'   \code{\link{dimnamesnamesp}} will be used on each array-like object to
##'   determine what indexing \code{colnames} it would have if melted; the
##'   result will have \code{\link{dimnamesnamesp}} containing all unique
##'   \code{\link{dimnamesnamesp}} from each of the arguments. For any two
##'   dimnames elements of with the same name selected from any of the arraylike
##'   arguments, either (a) the two should be identical, or (b) at least one
##'   should be trivial (NULL or repeated ""'s). Other types of objects can be
##'   wrapped in a list of class \code{"no_join"} using \code{\link{no_join}}
##'   and included in \code{arraylike.args}; they will be treated as scalars
##'   (constants), will not affect the number or naming of dimensions; the
##'   corresponding argument fed to \code{f} is always just the object wrapped
##'   inside the \code{no_join} object. This may be necessary if the same object
##'   should be used for all calls to \code{f} but the object appears to be
##'   array-like, to prevent new dimensions from being created or expected.
##'
##' @rdname map_join
##' @seealso no_join
##' @export
map_join_ = function(f, arraylike.args,
                     eltname.mismatch.behavior=c("stop","intersect"),
                     lapply_variant=parallel::mclapply, shuffle=TRUE,
                     show.progress=TRUE,
                     cache.prefix=NULL,
                     use.proxy=FALSE) {
  f <- match.fun(f)
  eltname.mismatch.behavior <- match.arg(eltname.mismatch.behavior)
  if (is.null(cache.prefix) && use.proxy) {
    stop ("Can only use a filesystem proxy if cache.prefix is specified (non-NULL).")
  }
  cache.dir =
    if (is.null(cache.prefix)) {
      NULL
    } else {
      dirname(cache.prefix)
    }
  if (!is.null(cache.dir) && !dir.exists(cache.dir)) {
    dir.create(cache.dir, recursive=TRUE)
  }
  ## todo allow for manually specifying index.dnnp?  make sure to drop exactly the index dimensions in each inputs
  index.dnp = list()
  for (arraylike.arg.i in seq_along(arraylike.args)) {
    arraylike.arg = arraylike.args[[arraylike.arg.i]]
    arraylike.arg.name = namesp(arraylike.args)[[arraylike.arg.i]]
    arg.dnp = dimnamesp(arraylike.arg)
    if (length(arg.dnp) == 0L && class(arraylike.arg) != "no_join") {
      ## warning (sprintf("arg %d had ndimp of 0L but was not marked no_join; marking as no_join now", arraylike.arg.i))
      arraylike.args[[arraylike.arg.i]] <- arraylike.arg <- no_join(arraylike.arg)
    }
    for (arg.dimension.i in seq_along(arg.dnp)) {
      dimension.name = names(arg.dnp)[[arg.dimension.i]]
      if (is.null(dimension.name) || dimension.name=="") {
        stop (sprintf('All dimensions must be (nontrivially) named.  Problem with arg %d (namep\'d "%s") dimension %d.',
                      arraylike.arg.i, arraylike.arg.name, arg.dimension.i))
      }
      dimension.eltnames = arg.dnp[[arg.dimension.i]]
      ## xxx may want to pre-allocate list with proper dnnp's
      existing.eltnames = index.dnp[[dimension.name]]
      if (is.null(existing.eltnames)) {
        ## no previous arg with this dimension; make dimension with these eltnames
        index.dnp[[dimension.name]] <- dimension.eltnames
      } else if (length(existing.eltnames) != length(dimension.eltnames) &&
                 (eltname.mismatch.behavior == "stop" ||
                  eltname.mismatch.behavior == "intersect" &&
                 (all(existing.eltnames=="") ||
                  all(dimension.eltnames==""))
                 )) {
          ## there is a length mismatch that we don't want to or can't fix by
          ## taking the eltnames intersection
          stop (sprintf('Inconsistent lengths found for dimension named "%s".  Length of dimension in arg %d (namep\'d "%s") (dimension %d): %d.  Previous length: %d.',
                        dimension.name, arraylike.arg.i, arraylike.arg.name, arg.dimension.i, length(dimension.eltnames), length(existing.eltnames)))
      } else if (length(existing.eltnames) != length(dimension.eltnames) &&
                 eltname.mismatch.behavior == "intersect" ||
                 any(existing.eltnames != dimension.eltnames)) {
        if (all(existing.eltnames=="")) {
          ## no pre-existing eltnames; assign
          index.dnp[[dimension.name]] <- dimension.eltnames
        } else if (all(dimension.eltnames=="")) {
          ## no need to do anything; should use existing eltnames
        } else {
          ## eltnames's don't match, and not because one was actually not named.
          if (eltname.mismatch.behavior == "stop") {
            stop (sprintf('dimnames associated with dimension named "%s" do not match.  Associated dimnames in arg %d (namep\'d "%s") (dimension %d): %s.  Previous associated dimnames: %s.',
                          dimension.name, arraylike.arg.i, arraylike.arg.name, arg.dimension.i, paste(utils::capture.output(dput(dimension.eltnames)), collapse=" "), paste(utils::capture.output(dput(existing.eltnames)), collapse=" ")))
          } else if (eltname.mismatch.behavior == "intersect") {
            index.dnp[[dimension.name]] <- intersect(existing.eltnames, dimension.eltnames)
          } else {
            stop (sprintf('Unrecognized/unhandled eltname.mismatch.behavior "%s".', eltname.mismatch.behavior))
          }
        }
      } else {
        ## no need to do anything; eltnames's match
      }
    }
  }

  index.dimension.lengths = vapply(index.dnp, length, 1L)
  result.from.arg.dimension.maps = lapply(
    arraylike.args, function(arraylike.arg) {
      match(dimnamesnamesp(arraylike.arg), names(index.dnp))
    }
  )
  result.to.arg.dimension.index.map.lists = lapply(
    seq_along(arraylike.args), function(arraylike.arg.i) {
      arraylike.arg.dim.indices = dimnames_or_inds(arraylike.args[[arraylike.arg.i]])
      corresponding.result.dim.indices = dimnamesp_to_dimindices(index.dnp[result.from.arg.dimension.maps[[arraylike.arg.i]]])
      Map(function(inds.a, inds.b) {
        if (class(inds.a)==class(inds.b)) {
          match(inds.a, inds.b)
        } else if (class(inds.a)=="character" && class(inds.b)=="integer" ||
                   class(inds.a)=="integer" && class(inds.b)=="integer") {
          stopifnot(length(inds.a)==length(inds.b))
          seq_along(inds.a)
        } else {
          stop ("Internal error: unable to translate dimension indices for the result to dimension indices of an argument.")
        }
      },
      corresponding.result.dim.indices, arraylike.arg.dim.indices)
    }
  )
  perm =
    if (shuffle) {
      sample.int(prod(index.dimension.lengths))
    } else {
      seq_len(prod(index.dimension.lengths))
    }
  length.perm = length(perm)
  result = lapply_variant(seq_along(perm), function(job.i) {
    result.elt.i = perm[[job.i]]
    if (show.progress && job.i == signif(job.i, 1L)) {
      print(paste0(job.i,"/",length.perm," ",Sys.time()))
    }
    indices = stats::setNames(as.vector(arrayInd(result.elt.i, index.dimension.lengths)), names(index.dimension.lengths))
    cache.file =
      if (is.null(cache.prefix)) {
        NULL
      } else {
        names.or.is = lapply(seq_along(indices), function(index.i) {
          index = indices[[index.i]]
          eltname = index.dnp[[index.i]][[index]]
          if (eltname != "") {
            gsub("/","-",eltname)
          } else {
            index
          }
        })
        paste0(cache.prefix,".",paste(names.or.is, collapse="."),".rds")
      }
    subresult =
      if (!is.null(cache.file) && file.exists(cache.file)) {
        if (use.proxy) {
          NULL
        } else {
          readRDS(cache.file)
        }
      } else {
        args =
          arraylike.args %>>%
          {stats::setNames(seq_along(.), names(.))} %>>%
          lapply(function(arraylike.arg.i) {
            arraylike.arg = arraylike.args[[arraylike.arg.i]]
            arraylike.arg.indices =
              mapply(
                magrittr::extract2,
                result.to.arg.dimension.index.map.lists[[arraylike.arg.i]],
                indices[result.from.arg.dimension.maps[[arraylike.arg.i]]]
              )
            arg =
              if (ndimp(arraylike.arg) == 0L) {
                stopifnot(class(arraylike.arg)=="no_join")
                arraylike.arg[[1L]]
              } else {
                arraylike.arg[t(as.matrix(arraylike.arg.indices))][[1L]]
              }
            arg
          })
        computed.subresult = do.call(f, args)
        if (!is.null(cache.file)) {
          saveRDS(computed.subresult, cache.file)
        }
        computed.subresult
      }
    if (use.proxy) {
      result.elt.i
    } else {
      subresult
    }
  })
  if (use.proxy) {
    ## xxx if there are errors, they are lost here...
    mode(result) <- "integer"
  }
  if (shuffle) {
    result[perm] <- result
  }
  if (length(index.dimension.lengths) > 0L) {
    dim(result) <- index.dimension.lengths
    dimnames(result) <- index.dnp
  } else {
    result <- result[[1L]]
  }
  if (use.proxy) {
    result <- array_proxy(cache.prefix, dimnamesp(result), result)
  }
  return (result)
}

##' @details \code{map_join} is provided as a potentially more convenient
##'   interface, eliminating the need to explicitly form a list of arraylike
##'   args; it converts the \code{...} arguments into a list and delegates to
##'   \code{map_join_}
##'
##' @param ... array-like objects with named dimnames, converted into an
##'   \code{arraylike.args} parameter using \code{list(...)}
##'
##' @examples
##' library("pipeR")
##' map_join(`*`, 2L,3L, lapply_variant=lapply)
##' map_join(`*`,
##'          with_dimnamesnames(2:3,"A"),
##'          with_dimnamesnames(1:3,"B")) %>>%
##'   {mode(.) <- "numeric"; .}
##' map_join(`*`,
##'          vector_as_named_array(2:3,"A",letters[1:2]),
##'          with_dimnamesnames(array(1:6,2:3), c("A","B"))) %>>%
##'   {mode(.) <- "numeric"; .}
##' cache.dir = tempfile()
##' map_join(`*`,
##'          vector_as_named_array(2:3,"A",letters[1:2]),
##'          with_dimnamesnames(1:3,"B"),
##'          cache.prefix=file.path(cache.dir,"outer_product")) %>>%
##' {mode(.) <- "numeric"; .}
##' arraylike.args = list(NULL
##'   , A=array(1:24,2:4) %>>%
##'     {dimnames(.) <- list(DA=paste0("S",1:2),DB=1:3,DC=1:4); .}
##'   , B=vector_as_named_array_(c(2.0,2.1), "DA", c("S1","S2"))
##'   , C=matrix(1:4, 2L,2L) %>>%
##'       {dimnames(.) <- list(DA=paste0("S",1:2),DA=paste0("S",1:2)); .}
##'   , D=142
##'   , E=1:5
##'   , F=vector_as_named_array_(11:14, "DC", 1:4)
##'   , G=c(S1=1,S2=2)
##'   , CP=matrix(1:4, 2L,2L) %>>%
##'       {dimnames(.) <- list(DA=paste0("S",1:2),DA=paste0("S",1:2)); .} %>>%
##'       magrittr::extract(1:2,2:1)
##'   , FP=vector_as_named_array_(11:15, "DC", 1:5)[c(5,3,4,2,1)]
##' )[-1L]
##' map_join_(list, arraylike.args[c("A","B","C","D","F")])[[DA="S1",DB=1L,DC=1L]]
##' map_join_(list, arraylike.args[c("A","B","C","D","F","CP","FP")],
##'           eltname.mismatch.behavior="intersect")[[DA="S1",DB=1L,DC=1L]]
##'
##' @rdname map_join
##' @export
map_join = function(f, ...,
                    eltname.mismatch.behavior=c("stop","intersect"),
                    lapply_variant=parallel::mclapply, shuffle=TRUE,
                    show.progress=TRUE,
                    cache.prefix=NULL,
                    use.proxy=FALSE) {
  arraylike.args = list(...)
  eltname.mismatch.behavior <- match.arg(eltname.mismatch.behavior)
  map_join_(f, arraylike.args,
            eltname.mismatch.behavior=eltname.mismatch.behavior,
            lapply_variant=lapply_variant, shuffle=shuffle,
            show.progress=show.progress,
            cache.prefix=cache.prefix,
            use.proxy=use.proxy)
}

##' Mark an array-like or other argument to be used like a constant/scalar in \code{\link{map_join}}
##'
##' @param x the argument to mark
##' @return an object of class \code{"map_join"} wrapping \code{x}
##'
##' @examples
##' map_join(`+`, with_dimnamesnames(1:3,"A"), with_dimnamesnames(1:3,"A"))
##' map_join(`+`, with_dimnamesnames(1:3,"A"), no_join(with_dimnamesnames(1:3,"A")))
##' @export
no_join = function(x) {
  structure(list(x), class="no_join")
}

extract_partial_= function(arraylike, index.sets,
                           dimension.missing.behavior=c("stop", "ignore"),
                           index.set.mismatch.behavior=c("stop", "intersect"),
                           drop=FALSE) {
  dimension.missing.behavior <- match.arg(dimension.missing.behavior)
  index.set.mismatch.behavior <- match.arg(index.set.mismatch.behavior)
  ## Handle cases where there are multiple index.set's for the same dimension:
  index_set_concord =
    switch(index.set.mismatch.behavior,
           stop=function(dimension.index.sets) {
             unique.dimension.index.sets = unique(dimension.index.sets)
             if (length(unique.dimension.index.sets) != 1L) {
               if (length(unique(sapply(unique.dimension.index.sets, class))) != 1L) {
                 stop ("Multiple classes of index sets provided for the same dimension; not supported.")
               } else {
                 stop ("Multiple index sets provided for the same dimension, but they are not all the same.")
               }
             }
             unique.dimension.index.sets[[1L]]
           },
           intersect=function(dimension.index.sets) {
             if (length(unique(sapply(dimension.index.sets, class))) != 1L) {
                 stop ("Multiple classes of index sets provided for the same dimension; not supported.")
             }
             Reduce(intersect, dimension.index.sets)
           },
           stop ("Unrecognized value of index.set.mismatch.behavior.")
           )
  simple.index.sets =
    split(index.sets, names(index.sets)) %>>%
    lapply(index_set_concord)
  dnp = dimnamesp(arraylike)
  extract.index.args = rep(list(TRUE), length(dnp))
  matched.dimension.is = match(names(simple.index.sets), names(dnp))
  if (any(is.na(matched.dimension.is))) {
    if (dimension.missing.behavior=="stop") {
      stop ("Dimension name specified in index.sets not present in dimnames(arraylike).")
    } else if (dimension.missing.behavior=="ignore") {
      retain.simple.index.set = !is.na(matched.dimension.is )
      simple.index.sets <- simple.index.sets[retain.simple.index.set]
      matched.dimension.is <- matched.dimension.is[retain.simple.index.set]
    } else {
      stop ("Unrecognized value for dimension.missing.behavior.")
    }
  }
  extract.index.args[match(names(simple.index.sets), names(dnp))] <- simple.index.sets
  do.call(`[`, c(list(arraylike), extract.index.args, list(drop=drop)))
}
## library("pipeR")
## array(1:24,2:4) %>>%
##   with_dimnames(list(A=letters[1:2], B=NULL, C=NULL)) %>>%
##   extract_partial_(list(A="a",B=2:3,B=2:3))
## array(1:24,2:4) %>>%
##   with_dimnames(list(A=letters[1:2], B=NULL, C=NULL)) %>>%
##   extract_partial_(list(A="a",B=2:3,B=1:2), index.set.mismatch.behavior="intersect")
## array(1:24,2:4) %>>%
##   with_dimnames(list(A=letters[1:2], B=NULL, C=NULL)) %>>%
##   extract_partial_(list(A="a",B=2:3,D=1), dimension.missing.behavior="ignore")

extract_along_= function(arraylike, along.arraylikes) {
  index.sets = dplyr::combine(lapply(along.arraylikes, dimnamesp))
  if (any(sapply(index.set, function(index.set) {
    all(index.set=="")
  }))) {
    stop ("All dimensions of along.arraylikes must be nontrivially named.")
  }
  return (extract_partial_(arraylike, index.sets))
}

select_dims = function(arraylike, dimensions,
                       select.behavior=c("retain", "drop"),
                       dimension.missing.behavior=c("stop", "ignore")) {
  select.behavior <- match.arg(select.behavior)
  dimension.missing.behavior <- match.arg(dimension.missing.behavior)
  if (is.character(dimensions)) {
    dimensions <- match(dimensions, dimnamesnamesp(arraylike))
    if (any(is.na(dimensions))) {
      if (dimension.missing.behavior=="stop") {
        stop ("Dimension name specified not present in dimnamesnamesp(arraylike).")
      } else if (dimension.missing.behavior=="ignore") {
        dimensions <- dimensions[!is.na(dimensions)]
      } else {
        stop ("Unrecognized value for dimension.missing.behavior.")
      }
    }
  } else if (is.logical(dimensions)) {
    dimensions <- which(dimensions)
  } else if (is.numeric(dimensions)){
    if (any(is.na(dimensions) | ! dimensions %in% seq_len(ndimp(arraylike)))) {
      stop ("All numeric dimension indices must be in seq_len(ndimp(arraylike)).")
    }
  } else {
    stop ("Unsupported class for argument =dimensions=.")
  }
  dimensions <- unique(dimensions)
  drop.dimensions =
    switch(select.behavior,
           retain=setdiff(seq_len(ndimp(arraylike)), dimensions),
           drop=dimensions,
           stop ("Unrecognized value for select.behavior."))
  if (any(dimp(arraylike)[drop.dimensions] != 1L)) {
    stop ("All dimensions to drop must be of size 1.")
  }
  original = arraylike
  new.dimp =
    if (length(drop.dimensions) == 0L) {
      dimp(original)
    } else {
      dim(original)[-drop.dimensions]
    }
  if (length(new.dimp) == 0L) {
    dim(arraylike) <- NULL
  } else {
    dim(arraylike) <- new.dimp
    if (length(drop.dimensions) == 0L) {
      dimnames(arraylike) <- dimnamesp(original)
    } else {
      dimnames(arraylike) <- dimnamesp(original)[-drop.dimensions]
    }
    if (select.behavior == "retain") {
      arraylike <- aperm(arraylike, order(dimensions))
    }
  }
  return (arraylike)
}
## select_dims(
##   with_dimnamesnames(structure(1, dim=rep(1,5)), letters[1:5])
## , letters[c(1,3,2,4,5)]) %>>%
##   dimnamesnamesp()

## todo make map_join work with data.frame inputs as well?
## todo side effect only variant
## todo multiple-cache version + special output class reading from cache rather than storing in memory
## todo check intersect behavior --- are the correct things actually selected and lined up?
## todo convenience variants of this for small vs. large operations?
## todo shorthand for . %>% apply(., seq_len(ndimp(.)), f, ...) plus special case for dim-less things? (and/or for vapply)
