// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// WeightedTabulateRcpp
Rcpp::NumericVector WeightedTabulateRcpp(Rcpp::IntegerVector bin, Rcpp::IntegerVector nbins, Rcpp::NumericVector w);
RcppExport SEXP _epiforecast_WeightedTabulateRcpp(SEXP binSEXP, SEXP nbinsSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type bin(binSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type nbins(nbinsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(WeightedTabulateRcpp(bin, nbins, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_epiforecast_WeightedTabulateRcpp", (DL_FUNC) &_epiforecast_WeightedTabulateRcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_epiforecast(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
