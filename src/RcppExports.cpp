// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// FOCuS_offline
int FOCuS_offline(NumericVector Y, double thres);
RcppExport SEXP _FOCuS_FOCuS_offline(SEXP YSEXP, SEXP thresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type thres(thresSEXP);
    rcpp_result_gen = Rcpp::wrap(FOCuS_offline(Y, thres));
    return rcpp_result_gen;
END_RCPP
}
// test
void test();
RcppExport SEXP _FOCuS_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test();
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FOCuS_FOCuS_offline", (DL_FUNC) &_FOCuS_FOCuS_offline, 2},
    {"_FOCuS_test", (DL_FUNC) &_FOCuS_test, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_FOCuS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
