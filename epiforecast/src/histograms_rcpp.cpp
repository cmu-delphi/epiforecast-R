// author_header begin
// Copyright (C) 2016 Logan C. Brooks
//
// This file is part of epiforecast.  Algorithms included in epiforecast were developed by Logan C. Brooks, David C. Farrow, Sangwon Hyun, Shannon Gallagher, Ryan J. Tibshirani, Roni Rosenfeld, and Rob Tibshirani (Stanford University), members of the Delphi group at Carnegie Mellon University.
//
// Research reported in this publication was supported by the National Institute Of General Medical Sciences of the National Institutes of Health under Award Number U54 GM088491. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health. This material is based upon work supported by the National Science Foundation Graduate Research Fellowship Program under Grant No. DGE-1252522. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the National Science Foundation. David C. Farrow was a predoctoral trainee supported by NIH T32 training grant T32 EB009403 as part of the HHMI-NIBIB Interfaces Initiative. Ryan J. Tibshirani was supported by NSF grant DMS-1309174.
// author_header end
// license_header begin
// epiforecast is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2 of the License.
//
// epiforecast is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with epiforecast.  If not, see <http://www.gnu.org/licenses/>.
// license_header end

#include <Rcpp.h>

//' Weighted, more \code{nbins}-restrictive version of \code{base::tabulate}
//'
//' @param bin integer-compatible vector; entries must be non-NA and between 1
//' and \code{nbins}; these indices denote entries in the result vector to which
//' the corresponding weights in \code{w} should be added
//' @param nbins single non-NA, non-negative integer; length of the vector to
//' return
//' @param w numeric-compatible vector of the same length as \code{bin}; weights
//' corresponding to the indices in \code{bin}
//' @return numeric vector of length \code{nbins}; the \code{i}th entry is like
//' \code{sum(w[bin==i])}, but with a naive summation algorithm
//'
//' @useDynLib epiforecast
//' @export
// [[Rcpp::export(name="weighted_tabulate")]]
Rcpp::NumericVector WeightedTabulateRcpp
(
 Rcpp::IntegerVector bin,
 Rcpp::IntegerVector nbins,
 Rcpp::NumericVector w
) {
  if (bin.size() != w.size()) {
    ::Rf_error("Invalid input: length(bin) != length(w)");
  }
  if (nbins.size() != 1 ||
      Rcpp::IntegerVector::is_na(* nbins.begin()) ||
      * nbins.begin() < 0) {
    ::Rf_error("Invalid input: nbins is not a single, non-NA, non-negative integer.");
  }
  Rcpp::IntegerVector::stored_type nbins0 = * nbins.begin();
  Rcpp::NumericVector result(nbins0); // 0-initialized
  Rcpp::IntegerVector::iterator bin_it = bin.begin();
  Rcpp::NumericVector::iterator w_it = w.begin();
  for(; bin_it != bin.end(); ++bin_it, ++w_it) {
    if (Rcpp::IntegerVector::is_na(*bin_it)) {
      ::Rf_error("Invalid input: any(is.na(bin))");
    }
    Rcpp::IntegerVector::stored_type index = *bin_it - 1;
    if (index < 0 || index >= nbins0) {
      ::Rf_error("Invalid input: !all(1L <= bin & bin <= nbins)");
    }
    result[index] += *w_it;
  }
  return result;
}

// clang-format off
// todo system includes
/* Local Variables: */
/* clang-format-style: "Google" */
/* flycheck-clang-language-standard: "c++14" */
/* flycheck-gcc-language-standard: "c++14" */
/* flycheck-clang-include-path: ("/usr/local/lib/R/site-library/Rcpp/include/" "/usr/share/R/include/") */
/* flycheck-gcc-include-path: ("/usr/local/lib/R/site-library/Rcpp/include/" "/usr/share/R/include/") */
/* End: */
// clang-format on
