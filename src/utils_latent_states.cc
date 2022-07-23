// /*
//  * R package stochvol by
//  *     Gregor Kastner Copyright (C) 2013-2018
//  *     Gregor Kastner and Darjus Hosszejni Copyright (C) 2019-
//  *  
//  *  This file is part of the R package stochvol: Efficient Bayesian
//  *  Inference for Stochastic Volatility Models.
//  *  
//  *  The R package stochvol is free software: you can redistribute it
//  *  and/or modify it under the terms of the GNU General Public License
//  *  as published by the Free Software Foundation, either version 2 or
//  *  any later version of the License.
//  *  
//  *  The R package stochvol is distributed in the hope that it will be
//  *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//  *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  *  General Public License for more details.
//  *  
//  *  You should have received a copy of the GNU General Public License
//  *  along with the R package stochvol. If that is not the case, please
//  *  refer to <http://www.gnu.org/licenses/>.
//  */
// 
// /*
//  * utils_latent_states.cc
//  * 
//  * Definitions of the functions declared in utils_latent_states.h.
//  * Documentation can also be found in utils_latent_states.h.
//  */
// 
// #include <RcppArmadillo.h>
// // #include "utils.h"
// #include "utils_latent_states.h"
// // #include "densities.h"
// #include <cmath>
// 
// namespace stochvol {
// 
// namespace fast_sv {
// 
// CholeskyTridiagonal cholesky_tridiagonal(
//     const arma::vec& omega_diag,
//     const double omega_offdiag) {
//   const int T = omega_diag.n_elem - 1;
//   arma::vec chol_diag(T+1);
//   arma::vec chol_offdiag(T+1);
//   chol_diag[0] = std::sqrt(omega_diag[0]);
//   for (int j = 1; j < T+1; j++) {
//     chol_offdiag[j-1] = omega_offdiag/chol_diag[j-1];
//     chol_diag[j] = std::sqrt(omega_diag[j]-chol_offdiag[j-1]*chol_offdiag[j-1]);
//   }
//   return {std::move(chol_diag), std::move(chol_offdiag)};
// }
// 
// arma::vec forward_algorithm(
//     const arma::vec& chol_diag,
//     const arma::vec& chol_offdiag,
//     const arma::vec& covector) {
//   const int T = chol_diag.n_elem - 1;
//   arma::vec htmp(T+1);
//   htmp[0] = covector[0]/chol_diag[0];
//   for (int j = 1; j < T+1; j++) {
//     htmp[j] = (covector[j] - chol_offdiag[j-1]*htmp[j-1])/chol_diag[j];
//   }
//   return htmp;
// }
// 
// arma::vec backward_algorithm(
//     const arma::vec& chol_diag,
//     const arma::vec& chol_offdiag,
//     const arma::vec& htmp) {
//   const int T = chol_diag.size() - 1;
//   arma::vec h(T+1);
//   h[T] = htmp[T] / chol_diag[T];
//   for (int j = T-1; j >= 0; j--) {
//     h[j] = (htmp[j] - chol_offdiag[j] * h[j+1]) / chol_diag[j];
//   }
//   return h;
// }
// 
// arma::uvec inverse_transform_sampling(
//     const arma::vec& mixprob,
//     const int T) {
//   arma::uvec r(T);
//   for (int j = 0; j < T; j++) {
//     int index = (10-1)/2;  // start searching in the middle
//     const double unnorm_cdf_value = R::unif_rand()*mixprob[9 + 10*j];  // current (non-normalized) value
//     bool larger = false;  // indicates that we already went up
//     bool smaller = false; // indicates that we already went down
//     while(true) {
//       if (unnorm_cdf_value > mixprob[index +  10*j]) {
//         index++;
//         if (smaller) {
//           break;
//         } else {
//           larger = true;
//         }
//       } else if (larger || index == 0) {
//         break;
//       } else {
//         index--;
//         smaller = true;
//       }
//     }
//     r[j] = index;
//   }
//   return r;
// }
// 
// arma::vec find_mixture_indicator_cdf(
//     const arma::vec& datanorm)  {
//   const int T = datanorm.n_elem;
//   arma::vec mixprob(10 * T);
//   for (int j = 0; j < T; j++) {  // TODO slow (10*T calls to exp)!
//     const int first_index = 10*j;
//     mixprob[first_index] = std::exp(mix_pre[0]-(datanorm[j]-mix_mean[0])*(datanorm[j]-mix_mean[0])*mix_2varinv[0]);
//     for (int r = 1; r < 10; r++) {
//       mixprob[first_index+r] = mixprob[first_index+r-1] + std::exp(mix_pre[r]-(datanorm[j]-mix_mean[r])*(datanorm[j]-mix_mean[r])*mix_2varinv[r]);
//     }
//   }
//   return mixprob;
// }
// 
// }  // END namespace fast_sv
// 
// 
// 
// }
