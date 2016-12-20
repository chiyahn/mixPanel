// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// LikelihoodMDP
SEXP LikelihoodMDP(Rcpp::NumericMatrix y_rcpp, Rcpp::NumericMatrix y_lagged_rcpp, Rcpp::NumericMatrix x_rcpp, Rcpp::NumericMatrix z_rcpp, Rcpp::NumericVector alpha_rcpp, Rcpp::NumericVector mu_rcpp, Rcpp::NumericVector sigma_rcpp, Rcpp::NumericMatrix rho_rcpp, Rcpp::NumericMatrix beta_rcpp, Rcpp::NumericVector gamma_rcpp);
RcppExport SEXP mixPanel_LikelihoodMDP(SEXP y_rcppSEXP, SEXP y_lagged_rcppSEXP, SEXP x_rcppSEXP, SEXP z_rcppSEXP, SEXP alpha_rcppSEXP, SEXP mu_rcppSEXP, SEXP sigma_rcppSEXP, SEXP rho_rcppSEXP, SEXP beta_rcppSEXP, SEXP gamma_rcppSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y_rcpp(y_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y_lagged_rcpp(y_lagged_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x_rcpp(x_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type z_rcpp(z_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha_rcpp(alpha_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu_rcpp(mu_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma_rcpp(sigma_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type rho_rcpp(rho_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type beta_rcpp(beta_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gamma_rcpp(gamma_rcppSEXP);
    rcpp_result_gen = Rcpp::wrap(LikelihoodMDP(y_rcpp, y_lagged_rcpp, x_rcpp, z_rcpp, alpha_rcpp, mu_rcpp, sigma_rcpp, rho_rcpp, beta_rcpp, gamma_rcpp));
    return rcpp_result_gen;
END_RCPP
}
// LikelihoodMDPInitialFixed
SEXP LikelihoodMDPInitialFixed(Rcpp::NumericMatrix y_rcpp, Rcpp::NumericMatrix y_lagged_rcpp, Rcpp::NumericMatrix x_rcpp, Rcpp::NumericMatrix z_rcpp, Rcpp::NumericVector alpha_rcpp, Rcpp::NumericVector mu_rcpp, Rcpp::NumericVector sigma_rcpp, Rcpp::NumericMatrix rho_rcpp, Rcpp::NumericMatrix beta_rcpp, Rcpp::NumericVector gamma_rcpp);
RcppExport SEXP mixPanel_LikelihoodMDPInitialFixed(SEXP y_rcppSEXP, SEXP y_lagged_rcppSEXP, SEXP x_rcppSEXP, SEXP z_rcppSEXP, SEXP alpha_rcppSEXP, SEXP mu_rcppSEXP, SEXP sigma_rcppSEXP, SEXP rho_rcppSEXP, SEXP beta_rcppSEXP, SEXP gamma_rcppSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y_rcpp(y_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type y_lagged_rcpp(y_lagged_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x_rcpp(x_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type z_rcpp(z_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha_rcpp(alpha_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu_rcpp(mu_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma_rcpp(sigma_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type rho_rcpp(rho_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type beta_rcpp(beta_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type gamma_rcpp(gamma_rcppSEXP);
    rcpp_result_gen = Rcpp::wrap(LikelihoodMDPInitialFixed(y_rcpp, y_lagged_rcpp, x_rcpp, z_rcpp, alpha_rcpp, mu_rcpp, sigma_rcpp, rho_rcpp, beta_rcpp, gamma_rcpp));
    return rcpp_result_gen;
END_RCPP
}
