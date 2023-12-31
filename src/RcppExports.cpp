// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rmnGibbsc
arma::mat rmnGibbsc(int m, arma::vec x0, arma::vec mu, arma::mat sigma);
RcppExport SEXP _SA23204191_rmnGibbsc(SEXP mSEXP, SEXP x0SEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmnGibbsc(m, x0, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA23204191_rmnGibbsc", (DL_FUNC) &_SA23204191_rmnGibbsc, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA23204191(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
