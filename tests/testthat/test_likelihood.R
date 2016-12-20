context("Likelihood estimation")

test_that("LikelihoodMDP", {
  
  set.seed(1234)
  s1 <- 0
  s2 <- 1
  N.fixed <- 50
  T.fixed <- 5
  p <- 2
  q <- 3
  x <- matrix(runif(T.fixed * N.fixed * p), nrow = T.fixed)
  z <- matrix(runif(T.fixed * N.fixed * q, 1, 3), nrow = T.fixed)
  x.empty <- matrix(0, nrow = T.fixed, ncol = N.fixed) # for p = 0 case
  z.empty <- matrix(0, nrow = T.fixed, ncol = N.fixed) # for q = 0 case
  
  theta1 <- GenerateMDPTheta(M = 1, s = s1)
  theta2 <- GenerateMDPTheta(M = 1, s = s2)
  theta3 <- GenerateMDPTheta(M = 2, s = s1)
  theta4 <- GenerateMDPTheta(M = 2, s = s2)
  theta5 <- GenerateMDPTheta(M = 3, s = s1, p = p)
  theta6 <- GenerateMDPTheta(M = 3, s = s2, p = p)
  theta7 <- GenerateMDPTheta(M = 3, s = s1, q = q)
  theta8 <- GenerateMDPTheta(M = 3, s = s2, q = q)
  theta9 <- GenerateMDPTheta(M = 3, s = s1, p = p, q = q)
  theta10 <- GenerateMDPTheta(M = 3, s = s2, p = p, q = q)
  sample1 <- GenerateMDPSample(theta1, N = N.fixed, T = T.fixed)
  sample2 <- GenerateMDPSample(theta2, N = N.fixed, T = T.fixed)
  sample3 <- GenerateMDPSample(theta3, N = N.fixed, T = T.fixed)
  sample4 <- GenerateMDPSample(theta4, N = N.fixed, T = T.fixed)
  sample5 <- GenerateMDPSample(theta5, N = N.fixed, T = T.fixed, x = x)
  sample6 <- GenerateMDPSample(theta6, N = N.fixed, T = T.fixed, x = x)
  sample7 <- GenerateMDPSample(theta7, N = N.fixed, T = T.fixed, z = z)
  sample8 <- GenerateMDPSample(theta8, N = N.fixed, T = T.fixed, z = z)
  sample9 <- GenerateMDPSample(theta9, N = N.fixed, T = T.fixed, x = x, z = z)
  sample10 <- GenerateMDPSample(theta10, N = N.fixed, T = T.fixed, x = x, z = z)
  
  LikelihoodMDP.wrapper <- function(sample, x, z)
  {
    theta <- (sample$MDP.model)$theta
    M <- length(theta$alpha)
    rho <- t(as.matrix(rep(0,M)))
    beta <- t(as.matrix(rep(0,M)))
    gamma <- as.matrix(0)
    if (!is.null(theta$rho))
      rho <- theta$rho
    if (!is.null(theta$beta))
      beta <- theta$beta
    if (!is.null(theta$gamma))
      gamma <- theta$gamma
    
    likelihood <- LikelihoodMDP(sample$y.sample, 
                                sample$y.lagged, x, z,
                                c(theta$alpha), c(theta$mu), c(theta$sigma),
                                rho, beta, c(gamma))
    likelihood.initial.fixed <- NULL
    if (is.null(sample$MDP.model$theta$rho)) # if s = 0
      likelihood.initial.fixed <- LikelihoodMDPInitialFixed(sample$y.sample, 
                                  sample$y.lagged, x, z,
                                  c(theta$alpha), c(theta$mu), c(theta$sigma),
                                  rho, beta, c(gamma))
    list(likelihood = likelihood, 
         likelihood.initial.fixed = likelihood.initial.fixed,
         sample = sample)
  }
  
  SanityCheck <- function(wrapped)
  {
    # check if likelihood is valid.
    expect_less_than(wrapped$likelihood, 0)
    
    # check if likelihoods are same for the case when s = 0
    if (is.null(wrapped$sample$MDP.model$theta$rho))
      expect_equal(wrapped$likelihood, wrapped$likelihood.initial.fixed)
  }
  
  SanityCheck(LikelihoodMDP.wrapper(sample1, x.empty, z.empty))
  SanityCheck(LikelihoodMDP.wrapper(sample2, x.empty, z.empty))
  SanityCheck(LikelihoodMDP.wrapper(sample3, x.empty, z.empty))
  SanityCheck(LikelihoodMDP.wrapper(sample4, x.empty, z.empty))
  SanityCheck(LikelihoodMDP.wrapper(sample5, x, z.empty))
  SanityCheck(LikelihoodMDP.wrapper(sample6, x, z.empty))
  SanityCheck(LikelihoodMDP.wrapper(sample7, x.empty, z))
  SanityCheck(LikelihoodMDP.wrapper(sample8, x.empty, z))
  SanityCheck(LikelihoodMDP.wrapper(sample9, x, z))
  SanityCheck(LikelihoodMDP.wrapper(sample10, x, z))
  
})
