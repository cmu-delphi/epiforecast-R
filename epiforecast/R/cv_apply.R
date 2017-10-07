##' @import parallel
##' @import R.utils
NULL

cv_apply_helper = function(train_data, test_data, indexer_list, fn, parallel_dim_i=0L, ...) {
  current_dim_i = length(indexer_list)
  stopifnot(dim(train_data)[current_dim_i] == dim(test_data)[current_dim_i])
  if (current_dim_i == 0L) {
    result = fn(train_data, test_data, ...)
    return (result)
  } else {
    current_dim_size = dim(train_data)[current_dim_i]
    current_dim_seq = seq_len(current_dim_size)
    current_dim_inpnames = dimnames(train_data)[[current_dim_i]]
    indexer_name = names(indexer_list)[current_dim_i]
    indexer_val = indexer_list[[current_dim_i]]
    switch(indexer_name,
           each={
             train_inds = current_dim_seq
             test_inds = train_inds
             current_dim_outnames = current_dim_inpnames
           },
           all={
             train_inds = list(TRUE)
             test_inds = train_inds
             current_dim_outnames = "all"
           },
           smear={
             train_inds = lapply(current_dim_seq, function(center_ind) {
               window_inds = center_ind + indexer_val
               window_inds <- window_inds[1L <= window_inds & window_inds <= current_dim_size]
             })
             test_inds = current_dim_seq
             current_dim_outnames = current_dim_inpnames
           },
           ablation={
             train_inds = lapply(current_dim_seq, function(left_out_ind) {
               current_dim_seq[-left_out_ind]
             })
             test_inds = train_inds
             current_dim_outnames = current_dim_inpnames
           },
           subsets={
             train_inds = indexer_val
             test_inds = train_inds
             current_dim_outnames = names(indexer_val)
           },
           loo={
             train_inds = lapply(current_dim_seq, function(train_out_ind) {
               current_dim_seq[-train_out_ind]
             })
             test_inds = current_dim_seq
             current_dim_outnames = current_dim_inpnames
           },
           oneahead={
             test_start_ind = indexer_val
             if (test_start_ind > current_dim_size) {
               stop ("oneahead argument greater than corresponding dimension width")
             }
             test_inds = test_start_ind-1L + seq_len(current_dim_size-test_start_ind+1L)
             train_inds = lapply(test_inds-1L, seq_len)
             current_dim_outnames = current_dim_inpnames[test_inds]
           },
           loo_oneahead={
             oneahead_start_ind = indexer_val
             if (oneahead_start_ind > current_dim_size) {
               stop ("oneahead argument greater than corresponding dimension width")
             }
             loo_test_inds = seq_len(oneahead_start_ind-1L)
             loo_train_inds = lapply(loo_test_inds, function(train_out_ind) {
               loo_test_inds[-train_out_ind]
             })
             oneahead_test_inds = oneahead_start_ind-1L + seq_len(current_dim_size-oneahead_start_ind+1L)
             oneahead_train_inds = lapply(oneahead_test_inds-1L, seq_len)
             train_inds = c(loo_train_inds, oneahead_train_inds)
             test_inds = c(loo_test_inds, oneahead_test_inds)
             current_dim_outnames = current_dim_inpnames[test_inds]
           },
           {
             stop("Unrecognized indexer name.")
           }
           )
    stopifnot(length(train_inds) == length(test_inds))
    current_dim_lapply = if (current_dim_i == parallel_dim_i) {
                           print("mclapply")
                           mclapply
                         } else {
                           lapply
                         }
    subresult.list =
      setNames(current_dim_lapply(seq_along(train_inds), function(indset_i) {
        train_data_inds = as.list(rep(TRUE, length(dim(train_data))))
        train_data_inds[[current_dim_i]] <- train_inds[[indset_i]]
        inner_train_data = do.call(`[`, c(list(train_data, drop=FALSE), train_data_inds))
        test_data_inds = as.list(rep(TRUE, length(dim(test_data))))
        test_data_inds[[current_dim_i]] <- test_inds[[indset_i]]
        inner_test_data = do.call(`[`, c(list(test_data, drop=FALSE), test_data_inds))
        cv_apply_helper(inner_train_data, inner_test_data, head(indexer_list,-1L), fn, parallel_dim_i=parallel_dim_i, ...)
      }), current_dim_outnames)
    result = simplify2array(subresult.list)
    return (result)
  }
}

##' \code{apply}-like function applying binary functions on training and test set selections
##'
##' @param data an array
##' @param indexer_list a named list with one entry per dimension of \code{data}
##'   specifying how to select training and test indices for that dimension; for
##'   the \code{i}th entry:

##'   * \code{each=NULL} slices both training and test data into
##'   \code{dim(data)[[i]]} pieces by selecting each index along the \code{i}th
##'   dimension; the corresponding output dimension width of
##'   \code{dim(data)[[i]]}

##'   * \code{all=NULL} performs no indexing along the \code{i}th dimension for
##'   either training or test data; the corresponding output dimension has width
##'   1 and name "all"

##'   * \code{smear=relative.indices} acts like \code{each=NULL}, but instead of
##'   each slices corresponding to a single index, allows for nearby indices to
##'   be included as well; for the \code{j}th slice, includes data for valid
##'   indices in \code{j+relative.indices}

##'   * \code{ablation=NULL} acts like \code{each=NULL}, but for the \code{j}th
##'   "slice", excludes, rather than selects, data corresponding to index
##'   \code{j}

##'   * \code{subsets=subset.list} indexes both training and test data based on
##'   the subsets in \code{subset.list}; each entry in \code{subset.list} should
##'   be a vector of indices (logical, integer, or character) into the
##'   \code{i}th dimension of \code{data}; the output dimension has width
##'   \code{length(subset.list)} and names \code{names(subset.list)}

##'   * \code{loo=NULL} performs leave-one-out cross-validation indexing: the
##'   training set like \code{ablation=NULL}, while the test set is indexed like
##'   \code{each=NULL}

##'   * \code{oneahead=test_start_ind} slices the test data by taking each index
##'   greater than or equal to the specified single \emph{integer} index
##'   \code{test_start_ind}; the test data corresponding to index i is paired
##'   with training data from indices strictly preceding i

##'   * \code{loo_oneahead=oneahead_start_ind} slices the test data by each
##'   index (similar to each); if the index i is less than the specified single
##'   \emph{integer} index \code{oneahead_start_ind}, then it is paired with
##'   training data from indices strictly preceding \code{oneahead_start_ind},
##'   excluding i; if i >= \code{oneahead_start_ind}, then the training data
##'   contains all indices strictly preceding i

##' @param fn a \code{function(training.slice, test.slice)} returning a scalar,
##'   vector, matrix, array, or list, with the same class and fixed structure
##'   for all inputs

##' @return an array with dimensionality equal to the sum of the
##'   dimensionless of the output of \code{fn} and of \code{data}

##' @md
##' @export
cv_apply = function(data, indexer_list, fn, parallel_dim_i=0L, ...) {
  ## If =data= is 1-D, convert it to an array so it will have non-NULL dim:
  if (is.null(dim(data))) {
    data <- as.array(data)
  }
  if (length(dim(data)) != length(indexer_list)) {
    stop ("Need exactly one indexer_list entry per dimension of data (or exactly 1 entry for vector data).")
  }
  if (is.null(names(indexer_list))) {
    stop ("Indexer types must be specified using names in indexer_list; indexer_list has no names.")
  }
  ## Use recursive cv_apply_helper to compute the entries of the result:
  result = cv_apply_helper(data, data, indexer_list, fn, parallel_dim_i=parallel_dim_i, ...)
  ## --- Adjust the class and dimnames: ---
  ## Make sure the result is an array:
  result <- as.array(result)
  ## Make sure the dimnames are a list, not NULL:
  if (is.null(dimnames(result))) {
    dimnames(result) <- rep(list(NULL), length(dim(result)))
  }
  ## Make sure the dimnames names are a character vector, no NULL:
  if (is.null(names(dimnames(result)))) {
    names(dimnames(result)) <- rep("",length(dim(result)))
  }
  ## If the original input =data= had dimnames names, assign them to the
  ## corresponding dimensions in the result:
  if (!is.null(names(dimnames(data)))) {
    names(dimnames(result))[tail(seq_along(dim(result)), length(dim(data)))] <- names(dimnames(data))
  }
  return (result)
}
