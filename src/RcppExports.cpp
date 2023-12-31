// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gram_schimdtCpp
arma::mat gram_schimdtCpp(arma::mat A);
RcppExport SEXP _fpta_gram_schimdtCpp(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(gram_schimdtCpp(A));
    return rcpp_result_gen;
END_RCPP
}
// Schur_wrapperCpp
Rcpp::List Schur_wrapperCpp(arma::mat projection);
RcppExport SEXP _fpta_Schur_wrapperCpp(SEXP projectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type projection(projectionSEXP);
    rcpp_result_gen = Rcpp::wrap(Schur_wrapperCpp(projection));
    return rcpp_result_gen;
END_RCPP
}
// predictCpp
arma::vec predictCpp(arma::mat B1, arma::mat B2, arma::mat Coef);
RcppExport SEXP _fpta_predictCpp(SEXP B1SEXP, SEXP B2SEXP, SEXP CoefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B2(B2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Coef(CoefSEXP);
    rcpp_result_gen = Rcpp::wrap(predictCpp(B1, B2, Coef));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fpta_gram_schimdtCpp", (DL_FUNC) &_fpta_gram_schimdtCpp, 1},
    {"_fpta_Schur_wrapperCpp", (DL_FUNC) &_fpta_Schur_wrapperCpp, 1},
    {"_fpta_predictCpp", (DL_FUNC) &_fpta_predictCpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fpta(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
