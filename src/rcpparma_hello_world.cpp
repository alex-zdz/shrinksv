// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "utils_latent_states.h"
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);
	                     
    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}


// utils_latent_states:


// [[Rcpp::export]]
Rcpp::List cholesky_tridiagonal(
        const arma::vec& omega_diag,
        const double omega_offdiag) {
    const int T = omega_diag.n_elem - 1;
    arma::vec chol_diag(T+1);
    arma::vec chol_offdiag(T+1);
    chol_diag[0] = std::sqrt(omega_diag[0]);
    for (int j = 1; j < T+1; j++) {
        chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
        chol_diag[j] = std::sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
    }
    // return {std::move(chol_diag), std::move(chol_offdiag)};
    Rcpp::List chol_ret;
    chol_ret["chol_diag"] = chol_diag;
    chol_ret["chol_offdiag"] = chol_offdiag;
    return chol_ret;
}


// [[Rcpp::export]]
arma::vec forward_algorithm(
        const arma::vec& chol_diag,
        const arma::vec& chol_offdiag,
        const arma::vec& covector) {
    const int T = chol_diag.n_elem - 1;
    arma::vec htmp(T+1);
    htmp[0] = covector[0]/chol_diag[0];
    for (int j = 1; j < T+1; j++) {
        htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
    }
    return htmp;
}
// [[Rcpp::export]]
arma::vec backward_algorithm(
        const arma::vec& chol_diag,
        const arma::vec& chol_offdiag,
        const arma::vec& htmp) {
    const int T = chol_diag.size() - 1;
    arma::vec h(T+1);
    h[T] = htmp[T] / chol_diag[T];
    for (int j = T-1; j >= 0; j--) {
        h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
    }
    return h;
}