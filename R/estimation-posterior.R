#' Estimates posterior probability of being in jth component for each ith agent.
#' @export
#' @title EstimateMDPPosteriors
#' @name EstimateMDPPosteriors
#' @param theta A list that represents the parameters of a model.
#' @param y (T-s) by N matrix that represents observed y variables 
#' @param y.lagged (T-s) by (N*min(1,s)) matrix that represents lagged y variables
#' @param x (T-s) by (N*p) matrix that represents observed exogeneous variables
#' that have coefficients that are different across components.
#' @param z (T-s) by (N*q) matrix that represents observed exogeneous variables
#' that have coefficients that are uniform across components.
#' @return N times M matrix whose (i,j)th element represents the posterior
#' probability of ith agent belonging to jth component.
EstimateMDPPosteriors <- function(theta, 
                                      y, y.lagged, x, z, 
                                   initial.fixed = FALSE)
{
  M <- length(theta$alpha)
  N <- ncol(y)
  T <- nrow(y)
  
  is.rho.switching <- TRUE
  if (!is.null(theta$rho))
    is.rho.switching <- (ncol(as.matrix(theta$rho)) > 1)
  is.sigma.switching <- (length(theta$sigma) > 1)
  
  rho <- matrix(rep(0, M), ncol = M)
  sigma <- theta$sigma
  beta <- matrix(rep(0, M), ncol = M)
  gamma <- as.matrix(0)
  
  # even if rho/sigma is not switching, make it like a switching parameter
  # by creating a matrix with a duplicated column so that we can use a single
  # code to estimate posterior probabilities.
  if (!is.rho.switching)
    rho <- matrix(replicate(M, rho), ncol = M)
  if (!is.sigma.switching)
    sigma <- replicate(M, sigma)
  
  if (!is.null(theta$rho))
    rho <- theta$rho
  
  if (!is.null(theta$beta))
  {
    x <- as.matrix(x)
    beta <- theta$beta
  }
  else
    x <- matrix(0, ncol = N, nrow = T)
  
  if (!is.null(theta$gamma))
  {
    z <- as.matrix(z)
    gamma <- theta$gamma
  }
  else
    z <- matrix(0, ncol = N, nrow = T)
  
  if (initial.fixed)
    return (PosteriorMDPInitialFixed(y, y.lagged, x, z, 
                         theta$alpha, theta$mu, sigma, rho, beta, gamma))
  return (PosteriorMDP(y, y.lagged, x, z, 
                       theta$alpha, theta$mu, sigma, rho, beta, gamma))
}