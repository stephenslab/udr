// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ed_iid_rcpp
arma::mat ed_iid_rcpp(const arma::mat& X, const arma::mat& U, const arma::vec& p);
RcppExport SEXP _udr_ed_iid_rcpp(SEXP XSEXP, SEXP USEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(ed_iid_rcpp(X, U, p));
    return rcpp_result_gen;
END_RCPP
}
// loglik_ud_iid_rcpp
double loglik_ud_iid_rcpp(const arma::mat& X, const arma::vec& w, const arma::cube& U, const arma::mat& V);
RcppExport SEXP _udr_loglik_ud_iid_rcpp(SEXP XSEXP, SEXP wSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_ud_iid_rcpp(X, w, U, V));
    return rcpp_result_gen;
END_RCPP
}
// loglik_ud_notiid_rcpp
double loglik_ud_notiid_rcpp(const arma::mat& X, const arma::vec& w, const arma::cube& U, const arma::cube& V);
RcppExport SEXP _udr_loglik_ud_notiid_rcpp(SEXP XSEXP, SEXP wSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_ud_notiid_rcpp(X, w, U, V));
    return rcpp_result_gen;
END_RCPP
}
// compute_loglik_matrix_iid_rcpp
arma::mat compute_loglik_matrix_iid_rcpp(const arma::mat& X, const arma::cube& U, const arma::mat& V);
RcppExport SEXP _udr_compute_loglik_matrix_iid_rcpp(SEXP XSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loglik_matrix_iid_rcpp(X, U, V));
    return rcpp_result_gen;
END_RCPP
}
// compute_loglik_matrix_notiid_rcpp
arma::mat compute_loglik_matrix_notiid_rcpp(const arma::mat& X, const arma::cube& U, const arma::cube& V);
RcppExport SEXP _udr_compute_loglik_matrix_notiid_rcpp(SEXP XSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_loglik_matrix_notiid_rcpp(X, U, V));
    return rcpp_result_gen;
END_RCPP
}
// compute_posterior_probs_iid_rcpp
arma::mat compute_posterior_probs_iid_rcpp(const arma::mat& X, const arma::vec& w, const arma::cube& U, const arma::mat& V);
RcppExport SEXP _udr_compute_posterior_probs_iid_rcpp(SEXP XSEXP, SEXP wSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_posterior_probs_iid_rcpp(X, w, U, V));
    return rcpp_result_gen;
END_RCPP
}
// compute_posterior_probs_notiid_rcpp
arma::mat compute_posterior_probs_notiid_rcpp(const arma::mat& X, const arma::vec& w, const arma::cube& U, const arma::cube& V);
RcppExport SEXP _udr_compute_posterior_probs_notiid_rcpp(SEXP XSEXP, SEXP wSEXP, SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_posterior_probs_notiid_rcpp(X, w, U, V));
    return rcpp_result_gen;
END_RCPP
}
// ted_rcpp
arma::mat ted_rcpp(const arma::mat& X, const arma::vec& p, double minval, unsigned int r);
RcppExport SEXP _udr_ted_rcpp(SEXP XSEXP, SEXP pSEXP, SEXP minvalSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type minval(minvalSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(ted_rcpp(X, p, minval, r));
    return rcpp_result_gen;
END_RCPP
}
// update_resid_covariance_rcpp
arma::mat update_resid_covariance_rcpp(const arma::mat& X, const arma::cube& U, const arma::mat& V, const arma::mat& P);
RcppExport SEXP _udr_update_resid_covariance_rcpp(SEXP XSEXP, SEXP USEXP, SEXP VSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(update_resid_covariance_rcpp(X, U, V, P));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_udr_ed_iid_rcpp", (DL_FUNC) &_udr_ed_iid_rcpp, 3},
    {"_udr_loglik_ud_iid_rcpp", (DL_FUNC) &_udr_loglik_ud_iid_rcpp, 4},
    {"_udr_loglik_ud_notiid_rcpp", (DL_FUNC) &_udr_loglik_ud_notiid_rcpp, 4},
    {"_udr_compute_loglik_matrix_iid_rcpp", (DL_FUNC) &_udr_compute_loglik_matrix_iid_rcpp, 3},
    {"_udr_compute_loglik_matrix_notiid_rcpp", (DL_FUNC) &_udr_compute_loglik_matrix_notiid_rcpp, 3},
    {"_udr_compute_posterior_probs_iid_rcpp", (DL_FUNC) &_udr_compute_posterior_probs_iid_rcpp, 4},
    {"_udr_compute_posterior_probs_notiid_rcpp", (DL_FUNC) &_udr_compute_posterior_probs_notiid_rcpp, 4},
    {"_udr_ted_rcpp", (DL_FUNC) &_udr_ted_rcpp, 4},
    {"_udr_update_resid_covariance_rcpp", (DL_FUNC) &_udr_update_resid_covariance_rcpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_udr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
