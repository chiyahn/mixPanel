#' Estimates a Fisher information matrix given parameters and a sample.
#' @export
#' @title EstimateMDPFisher
#' @name EstimateMDPFisher
#' @param y (T-s) by N matrix that represents observed y variables 
#' @param y.lagged (T-s) by (N*min(1,s)) matrix that represents lagged y variables
#' @param x (T-s) by (N*p) matrix that represents observed exogeneous variables
#' that have coefficients that are different across components.
#' @param z (T-s) by (N*q) matrix that represents observed exogeneous variables
#' that have coefficients that are uniform across components.
#' @param eps epsilon value used to compute a numerical score function
#' @param initial.fixed Determines whether initial values are treated as fixed.
#' @return size(theta) times size(theta) matrix that represents a Fisher
#' information matrix
EstimateMDPFisher <- function(theta, y, y.lagged,
                                      x, z,
                                      eps = 1e-6, initial.fixed = FALSE)
{
  # assume that rho/sigma is different across components.
  is.rho.switching <- TRUE
  is.sigma.switching <- TRUE
  
  T <- nrow(y)
  N <- ncol(y)
  n <- length(y)
  M <- length(theta$alpha)
  s <- 0
  p <- 0
  q <- 0
  if (!is.null(theta$rho))
    s <- nrow(as.matrix(theta$rho))
  
  if (is.null(theta$beta))
    x <- matrix(rep(0,(T * N)), ncol = N)
  else
    p <- nrow(theta$beta)
  
  if (is.null(theta$gamma))
    z <- matrix(rep(0,(T * N)), ncol = N)
  else
    q <- nrow(as.matrix(theta$gamma))
  
  # this holds only for univariate observation cases.
  mu.index <- M  
  sigma.index <- M + M
  
  rho.index <- ifelse(is.sigma.switching, M, 1) + sigma.index
  beta.index <- s * ifelse(is.rho.switching, M, 1) + rho.index
  gamma.index <- p * M + beta.index
  
  LikelihoodsFunction <- LikelihoodsMDP
  if (initial.fixed)
    LikelihoodsFunction <- LikelihoodsMDPInitialFixed
  
  ObjectiveLogLikelihoods <- function(theta.vectorized)
  {
    alpha <- 1
    if (M > 1)
    {
      alpha <- c(theta.vectorized[1:(M - 1)])
      alpha <- c(alpha, (1-sum(alpha)))
    }
    
    mu <- theta.vectorized[mu.index:(sigma.index - 1)]
    sigma <- theta.vectorized[sigma.index:(rho.index - 1)]
    if (!is.sigma.switching) # make it as a switching parameter if not.
      sigma <- rep(sigma, M)
    
    rho <- t(rep(0,M))
    beta <- t(rep(0,M))
    gamma <- 0
    
    # i.e. rho exists
    if (s > 0)
    {
      rho <- theta.vectorized[rho.index:(beta.index - 1)]
      if (!is.rho.switching) # make it as a switching parameter if not.
        rho <- rep(rho, M)
      rho <- matrix(rho, ncol = M)
    }
    
    # i.e. beta exists
    if (p > 0)
      beta <- matrix(theta.vectorized[beta.index:(gamma.index - 1)],
                     ncol = M)
    # i.e. gamma exists
    if (q > 0)
      gamma <- theta.vectorized[gamma.index:length(theta.vectorized)]
    
    return (LikelihoodsFunction(y, y.lagged, x, z,
                                alpha, mu, sigma, rho, beta, gamma))
  }
  
  
  arg <- MDPThetaToReducedColumn(theta)
  # Defines a step (make sure it does not bind with the ub/lb)
  h <- pmax(eps, abs(arg)) * eps ^ {2/3}
  argh <- arg + h
  h <- argh - arg
  h.diag <- diag(h)
  
  G <- matrix(0, nrow = length(arg), ncol = N)
  H <- matrix(0, nrow = length(arg), ncol = length(arg))
  
  for (i in 1:length(arg))
    G[i,] <- (ObjectiveLogLikelihoods(arg + h.diag[i,]) - 
                ObjectiveLogLikelihoods(arg - h.diag[i,])) /
    (2 * h[i])
  
  for (k in 1:N)
    H <- H + G[,k] %*% t(G[,k])
  
  return (H / N)
}