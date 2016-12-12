/*
// Compute likelihood of MDP models.
// Author: Chiyoung Ahn
*/

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

const double LOG2PI_OVERTWO = 0.91893853320467274178; // (log(2*pi) / 2)

// Returns likelihood of univariate MS-AR model where rho is switching.
// Note that even if rho is non-switching, setting rho as a s by M matrix with
// repeated column of the original rho will give you the likelihood for
// MS-AR model with non-switching rho.
// [[Rcpp::export]]
SEXP LikelihoodMDP (Rcpp::NumericMatrix y_rcpp,
					Rcpp::NumericMatrix y_lagged_rcpp,
					Rcpp::NumericMatrix x_rcpp,
					Rcpp::NumericMatrix z_rcpp,
					Rcpp::NumericVector alpha_rcpp,
					Rcpp::NumericVector mu_rcpp,
					Rcpp::NumericVector sigma_rcpp,
					Rcpp::NumericMatrix rho_rcpp,
					Rcpp::NumericMatrix beta_rcpp,
					Rcpp::NumericVector gamma_rcpp)
{
	arma::mat    y(y_rcpp.begin(),
                 y_rcpp.nrow(),
                 y_rcpp.ncol(), false);
	arma::mat    y_lagged(y_lagged_rcpp.begin(),
								y_lagged_rcpp.nrow(),
								y_lagged_rcpp.ncol(), false);
	arma::mat    x(x_rcpp.begin(),
								x_rcpp.nrow(),
								x_rcpp.ncol(), false);
	arma::mat    z(z_rcpp.begin(),
								z_rcpp.nrow(),
								z_rcpp.ncol(), false);
	arma::colvec alpha(alpha_rcpp.begin(),
								alpha_rcpp.size(), false);
	arma::colvec mu(mu_rcpp.begin(), mu_rcpp.size(), false);
	arma::colvec sigma(sigma_rcpp.begin(), sigma_rcpp.size(), false);
	arma::mat    rho(rho_rcpp.begin(),
                   rho_rcpp.nrow(), rho_rcpp.ncol(), false);

	arma::mat    beta(beta_rcpp.begin(),
								beta_rcpp.nrow(),
								beta_rcpp.ncol(), false);
	arma::colvec gamma(gamma_rcpp.begin(),
								gamma_rcpp.size(), false);

	int T = y.n_rows;
	int N = y.n_cols;
	int M = alpha.n_rows;
	int s = rho.n_rows;
	int p = beta.n_rows;
	int q = gamma.n_rows;

  double likelihood = 0;

	// partition blocks
	arma::mat* y_lagged_blocks = new arma::mat[N];
	arma::mat* x_blocks = new arma::mat[N];
	arma::mat* z_blocks = new arma::mat[N];
	for (int i = 0; i < N; i++)
	{
		int y_block_first = i * s;
		int x_block_first = i * p;
		int z_block_first = i * q;
		y_lagged_blocks[i] = y_lagged.cols(y_block_first, y_block_first + s - 1);
		x_blocks[i] = x.cols(x_block_first, x_block_first + p - 1);
		z_blocks[i] = z.cols(z_block_first, z_block_first + q - 1);
	}

	for (int i = 0; i < N; i++)
	{
		// initial setting; keep track of minimum value and its index
		// to divide everything by the min. value in order to prevent
		// possible numerical errors when computing posterior probs.
		int min_index = -1;
		double min_value = std::numeric_limits<double>::infinity();
		double* ratios = new double[M];
		double* likelihoods_i = new double[M];
		double likelihood_i = 0;

		// compute likelihoods for each block
		for (int j = 0; j < M; j++)
		{
			likelihoods_i[j] = 0;
			for (int t = 0; t < T; t++)
			{
				arma::colvec likelihood_ijt = y(t,i) -
					y_lagged_blocks[i].row(t) * rho.col(j) -
		      x_blocks[i].row(t) * beta.col(j) -
		      z_blocks[i].row(t) * gamma - mu(j); // explicit gluing
		    likelihood_ijt *= likelihood_ijt;

				// you can subtract SQRT2PI in the final log-likelihood.
				likelihoods_i[j] += likelihood_ijt(0);
			}

			likelihoods_i[j] /= (2 * (sigma(j) * sigma(j)));

			if (min_value > likelihoods_i[j])
			{
				min_value = likelihoods_i[j];
				min_index = j;
			}

			ratios[j] = alpha(j) / sigma(j);
		}

		// minimum back
		for (int j = 0; j < M; j++)
		{
			if (j == min_index)
				likelihoods_i[j] = 1.0;
			else
				likelihoods_i[j] = (ratios[j] / ratios[min_index]) *
											exp(min_value - likelihoods_i[j]);
			likelihood_i += likelihoods_i[j];
		}

		likelihood += log(likelihood_i) - min_value + log(ratios[min_index]);

		delete[] ratios;
		delete[] likelihoods_i; // clear memory
	}

	// subtract SQRT2PIs that have been missed in iteration
	likelihood -= T * N * LOG2PI_OVERTWO;

	// clear memory for blocks
	delete[] y_lagged_blocks;
	delete[] x_blocks;
	delete[] z_blocks;

	return (wrap(likelihood));
}
