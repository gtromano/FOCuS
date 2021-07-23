// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// FOCuS
List FOCuS(Rcpp::Function dataGen, const double thres, const double& mu0, std::list<double>& grid, const double& K);
RcppExport SEXP _FOCuS_FOCuS(SEXP dataGenSEXP, SEXP thresSEXP, SEXP mu0SEXP, SEXP gridSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Function >::type dataGen(dataGenSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< std::list<double>& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const double& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(FOCuS(dataGen, thres, mu0, grid, K));
    return rcpp_result_gen;
END_RCPP
}
// FOCuS_offline
List FOCuS_offline(NumericVector Y, const double thres, const double& mu0, std::list<double>& grid, const double& K);
RcppExport SEXP _FOCuS_FOCuS_offline(SEXP YSEXP, SEXP thresSEXP, SEXP mu0SEXP, SEXP gridSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< std::list<double>& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const double& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(FOCuS_offline(Y, thres, mu0, grid, K));
    return rcpp_result_gen;
END_RCPP
}
// FOCuS_melk
List FOCuS_melk(NumericVector Y, const double thres, const double& mu0, std::list<double>& grid, const double& K);
RcppExport SEXP _FOCuS_FOCuS_melk(SEXP YSEXP, SEXP thresSEXP, SEXP mu0SEXP, SEXP gridSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< std::list<double>& >::type grid(gridSEXP);
    Rcpp::traits::input_parameter< const double& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(FOCuS_melk(Y, thres, mu0, grid, K));
    return rcpp_result_gen;
END_RCPP
}
// MOSUM_offline_kirch
List MOSUM_offline_kirch(NumericVector Y, const double thres, std::vector<int> W);
RcppExport SEXP _FOCuS_MOSUM_offline_kirch(SEXP YSEXP, SEXP thresSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(MOSUM_offline_kirch(Y, thres, W));
    return rcpp_result_gen;
END_RCPP
}
// MOSUM_offline_kirch2
List MOSUM_offline_kirch2(NumericVector Y, const double thres, std::vector<int> W);
RcppExport SEXP _FOCuS_MOSUM_offline_kirch2(SEXP YSEXP, SEXP thresSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(MOSUM_offline_kirch2(Y, thres, W));
    return rcpp_result_gen;
END_RCPP
}
// PageCUSUM_offline
List PageCUSUM_offline(NumericVector Y, const double thres, const double& mu0, std::vector<double>& grid);
RcppExport SEXP _FOCuS_PageCUSUM_offline(SEXP YSEXP, SEXP thresSEXP, SEXP mu0SEXP, SEXP gridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type grid(gridSEXP);
    rcpp_result_gen = Rcpp::wrap(PageCUSUM_offline(Y, thres, mu0, grid));
    return rcpp_result_gen;
END_RCPP
}
// CUSUM_offline
List CUSUM_offline(NumericVector Y, const double thres, const double& mu0);
RcppExport SEXP _FOCuS_CUSUM_offline(SEXP YSEXP, SEXP thresSEXP, SEXP mu0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu0(mu0SEXP);
    rcpp_result_gen = Rcpp::wrap(CUSUM_offline(Y, thres, mu0));
    return rcpp_result_gen;
END_RCPP
}
// simpleMelkman
List simpleMelkman(NumericVector x, bool onlyPrune, bool exportInR);
RcppExport SEXP _FOCuS_simpleMelkman(SEXP xSEXP, SEXP onlyPruneSEXP, SEXP exportInRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyPrune(onlyPruneSEXP);
    Rcpp::traits::input_parameter< bool >::type exportInR(exportInRSEXP);
    rcpp_result_gen = Rcpp::wrap(simpleMelkman(x, onlyPrune, exportInR));
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
// YuCUSUM_offline
List YuCUSUM_offline(NumericVector Y, const double thres);
RcppExport SEXP _FOCuS_YuCUSUM_offline(SEXP YSEXP, SEXP thresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type thres(thresSEXP);
    rcpp_result_gen = Rcpp::wrap(YuCUSUM_offline(Y, thres));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FOCuS_FOCuS", (DL_FUNC) &_FOCuS_FOCuS, 5},
    {"_FOCuS_FOCuS_offline", (DL_FUNC) &_FOCuS_FOCuS_offline, 5},
    {"_FOCuS_FOCuS_melk", (DL_FUNC) &_FOCuS_FOCuS_melk, 5},
    {"_FOCuS_MOSUM_offline_kirch", (DL_FUNC) &_FOCuS_MOSUM_offline_kirch, 3},
    {"_FOCuS_MOSUM_offline_kirch2", (DL_FUNC) &_FOCuS_MOSUM_offline_kirch2, 3},
    {"_FOCuS_PageCUSUM_offline", (DL_FUNC) &_FOCuS_PageCUSUM_offline, 4},
    {"_FOCuS_CUSUM_offline", (DL_FUNC) &_FOCuS_CUSUM_offline, 3},
    {"_FOCuS_simpleMelkman", (DL_FUNC) &_FOCuS_simpleMelkman, 3},
    {"_FOCuS_test", (DL_FUNC) &_FOCuS_test, 0},
    {"_FOCuS_YuCUSUM_offline", (DL_FUNC) &_FOCuS_YuCUSUM_offline, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_FOCuS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
