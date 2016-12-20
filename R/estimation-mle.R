EstimateMDP.MLE.nloptr <- function(thetas, y, y.lagged, 
                               x, z,
                               epsilon = 0.001, maxit = 500,
                               sigma.min = 0.05, alpha.eps = 10e-4,
                               initial.fixed = FALSE)
{
  # assume that rho/sigma is different across components.
  is.rho.switching <- TRUE
  is.sigma.switching <- TRUE
  
  theta <- thetas[[1]]
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
  
  # hard constraints
  alpha.lb <- NULL
  alpha.ub <- NULL
  if (M > 1)
  {  
    alpha.lb <- rep(alpha.eps, (M-1))
    alpha.ub <- rep((1-alpha.eps), (M-1))
  }
  sigma.lb <- rep(sigma.min, ifelse(is.sigma.switching, M, 1))
  
  # dynamically defines a function that transforms a vectorized theta to a list
  ReducedColumnToTheta <- function(theta.vectorized)
  {
    alpha <- 1
    if (M > 1)
    {
      alpha <- theta.vectorized[1:(M - 1)]
      alpha <- c(alpha, (1 - sum(alpha)))
    }

    mu <- theta.vectorized[mu.index:(sigma.index - 1)]
    sigma <- theta.vectorized[sigma.index:(rho.index - 1)]
    
    rho <- NULL
    beta <- NULL
    gamma <- NULL
    
    if (s > 0)
    {
      rho <- theta.vectorized[rho.index:(beta.index - 1)]
      if (is.rho.switching)
        rho <- matrix(rho, ncol = M)
    }
    if (p > 0)
      beta <- matrix(theta.vectorized[beta.index:(gamma.index - 1)],
                                ncol = M)
    if (q > 0)
      gamma <- theta.vectorized[gamma.index:length(theta.vectorized)]
    
    return (list(alpha = alpha, mu = mu, sigma = sigma,
            rho = rho, beta = beta, gamma = gamma))
  }
  
  LikelihoodFunction <- LikelihoodMDP
  if (initial.fixed)
    LikelihoodFunction <- LikelihoodMDPInitialFixed
    
  ObjectiveLogLikelihood <- function(theta.vectorized)
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
    
    # slsqp solves a minimization problem;
    # take a negative value to turn the problem into max. problem
    return (-LikelihoodFunction(y, y.lagged, x, z,
                  alpha, mu, sigma, rho, beta, gamma))
  }
  
  # defines a constraint on alpha.
  ConstraintAlpha <- function(theta.vectorized)
  {
    return (alpha.ub - sum(theta.vectorized[1:(M-1)]))
  }
  
  SLSQPStep <- function(theta.vectorized)
  {
    # sanity check; if a candidate contains a singularity, you must not use it.
    if (anyNA(theta.vectorized) || is.null(theta.vectorized))
      return (list(convergence = -3, value = -Inf))
    
    mu <- theta.vectorized[mu.index:(sigma.index - 1)]
    sigma <- theta.vectorized[sigma.index:(rho.index - 1)]
    
    
    mu.lb <- pmin(-1, mu * (1 - 0.8 * sign(mu)))
    mu.ub <- pmax(1, mu * (1 + 0.8 * sign(mu)))
    sigma.ub <- sigma * 8
    rho.lb <- NULL
    rho.ub <- NULL
    if (s > 0) 
    {
      rho <- theta.vectorized[rho.index:(gamma.index - 1)]
      rho.lb <- pmin(-1, rho * (1 - 0.5 * sign(rho)))
      rho.ub <- pmax(1, rho * (1 + 0.5 * sign(rho)))
    }
    
    lb <- c(alpha.lb, mu.lb, sigma.lb, rho.lb)
    lb <- c(lb, rep(-Inf, (length(theta.vectorized) - length(lb))))
    ub <- c(alpha.ub, mu.ub, sigma.ub, rho.ub)
    ub <- c(ub, rep(Inf, (length(theta.vectorized) - length(ub))))
    
    # sanity check for derivatives.
    if (!NLOPTRSanityCheck(x0 = theta.vectorized, fn = ObjectiveLogLikelihood))
      return (list (convergence = -Inf, value = -Inf))
    
    result <- nloptr::slsqp(theta.vectorized,
                            fn = ObjectiveLogLikelihood,
                            lower = lb, upper = ub, hin = ConstraintAlpha,
                            control = list(maxeval = maxit, ftol_abs = epsilon))
    result$value <- -result$value # take negative back to make it actual log-lik.
    return (result)
  }
  
  nloptr.thetas.matrix <- sapply(thetas,
                                 function (theta) MDPThetaToReducedColumn(theta))
  nloptr.results <- apply(nloptr.thetas.matrix, 2, SLSQPStep)
  nloptr.convergence <- unlist(lapply(nloptr.results, "[[", "convergence"))
  nloptr.likelihoods <- unlist(lapply(nloptr.results, "[[", "value"))
  nloptr.likelihoods[!is.finite(nloptr.likelihoods)] <- -Inf # abnormal values
  nloptr.likelihoods[nloptr.convergence < 0] <- -Inf # non-convergence
  
  # extract the one that returns the best log.likelihood
  nloptr.result <- nloptr.results[[(which(nloptr.likelihoods==
                                        max(nloptr.likelihoods))[1])]]
  # if (!is.finite(long.result$value))
  #   return (list (succeeded = FALSE))
  return (list(theta = ReducedColumnToTheta(nloptr.result$par),
               likelihood = nloptr.result$value,
               nloptr.results = nloptr.results,
               succeeded = TRUE))
}

NLOPTRSanityCheck <- function (x0, fn)
{
  # sanity check for derivatives
  heps <- .Machine$double.eps^(1/3) # epsilon used in nloptr package
  n <- length(x0)
  hh <- diag(heps, n)
  gr <- numeric(n)
  for (i in 1:n)
    if (is.na(fn(x0 + hh[,i])) || is.na(fn(x0 - hh[,i])))
      return (FALSE)
  return (TRUE)
}