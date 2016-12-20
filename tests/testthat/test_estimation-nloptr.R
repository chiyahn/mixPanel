context("Parameter estimation - nloptr")

test_that("using nloptr", {
  
  set.seed(1234)
  s1 <- 0
  s2 <- 1
  N.fixed <- 30
  T.fixed <- 20
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
  
  nloptr.wrapper <- function(sample, x = NULL, z = NULL)
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
    
    nloptr.result <- EstimateMDP.MLE.nloptr(thetas = list(theta),
                                            y = sample$y.sample, y.lagged = sample$y.lagged, 
                                            x = x, z = z, initial.fixed = FALSE)
    nloptr.result.fixed <- EstimateMDP.MLE.nloptr(thetas = list(theta),
                          y = sample$y.sample, y.lagged = sample$y.lagged, 
                          x = x, z = z, initial.fixed = TRUE)
    
    if (is.null(x))
      x <- x.empty
    if (is.null(z))
      z <- z.empty
    
    likelihood <- LikelihoodMDP(sample$y.sample, 
                                sample$y.lagged, x, z,
                                c(theta$alpha), c(theta$mu), c(theta$sigma),
                                rho, beta, c(gamma))
    likelihood.fixed <- LikelihoodMDPInitialFixed(sample$y.sample, 
                                                  sample$y.lagged, x, z,
                                                  c(theta$alpha), c(theta$mu), c(theta$sigma),
                                                  rho, beta, c(gamma))
    
    list(initial.guess = theta,
         nloptr.result = nloptr.result,
         nloptr.result.fixed = nloptr.result.fixed,
         likelihood = likelihood,
         likelihood.fixed = likelihood.fixed)
    
  }
  
  SanityCheck <- function(result)
  {
    initial.guess <- result$initial.guess
    
    # check if likelihood makes sense
    expect_less_than(result$nloptr.result$likelihood, 0)
    expect_less_than(result$nloptr.result.fixed$likelihood, 0)
    expect_less_than(result$likelihood, 0)
    expect_less_than(result$likelihood.fixed, 0)
    
    # check if likelihood increased after optimization
    expect_true(result$nloptr.result$likelihood > result$likelihood)
    expect_true(result$nloptr.result.fixed$likelihood > result$likelihood.fixed)
    
    # check if all returned parameters are in the same format
    expect_equal(length(result$nloptr.result$theta$alpha), length(initial.guess$alpha))
    expect_equal(length(result$nloptr.result$theta$mu), length(initial.guess$mu))
    expect_equal(length(result$nloptr.result$theta$sigma), length(initial.guess$sigma))
    if (!is.null(initial.guess$rho))
    {
      expect_equal(dim(result$nloptr.result$theta$rho), dim(initial.guess$rho))
      expect_equal(dim(result$nloptr.result$theta$rho), dim(initial.guess$rho))
    }
    if (!is.null(initial.guess$beta))
    {
      expect_equal(dim(result$nloptr.result$theta$beta), dim(initial.guess$beta))
      expect_equal(dim(result$nloptr.result$theta$beta), dim(initial.guess$beta))
    }
    if (!is.null(initial.guess$gamma))
    {
      expect_equal(length(result$nloptr.result$theta$gamma), length(initial.guess$gamma))
      expect_equal(length(result$nloptr.result$theta$gamma), length(initial.guess$gamma))
    }
    
    expect_equal(length(result$nloptr.result.fixed$theta$alpha), length(initial.guess$alpha))
    expect_equal(length(result$nloptr.result.fixed$theta$mu), length(initial.guess$mu))
    expect_equal(length(result$nloptr.result.fixed$theta$sigma), length(initial.guess$sigma))
    if (!is.null(initial.guess$rho))
    {
      expect_equal(dim(result$nloptr.result.fixed$theta$rho), dim(initial.guess$rho))
      expect_equal(dim(result$nloptr.result.fixed$theta$rho), dim(initial.guess$rho))
    }
    if (!is.null(initial.guess$beta))
    {
      expect_equal(dim(result$nloptr.result.fixed$theta$beta), dim(initial.guess$beta))
      expect_equal(dim(result$nloptr.result.fixed$theta$beta), dim(initial.guess$beta))
    }
    if (!is.null(initial.guess$gamma))
    {
      expect_equal(length(result$nloptr.result.fixed$theta$gamma), length(initial.guess$gamma))
      expect_equal(length(result$nloptr.result.fixed$theta$gamma), length(initial.guess$gamma))
    }
    
    # check if parameters are valid
    expect_equal(sum(result$nloptr.result$theta$alpha), 1)
    expect_equal(sum(result$nloptr.result.fixed$theta$alpha), 1)
    
  }
  
  SanityCheck(nloptr.wrapper(sample1))
  SanityCheck(nloptr.wrapper(sample2))
  SanityCheck(nloptr.wrapper(sample3))
  SanityCheck(nloptr.wrapper(sample4))
  SanityCheck(nloptr.wrapper(sample5, x = x))
  SanityCheck(nloptr.wrapper(sample6, x = x))
  SanityCheck(nloptr.wrapper(sample7, z = z))
  SanityCheck(nloptr.wrapper(sample8, z = z))
  SanityCheck(nloptr.wrapper(sample9, x = x, z = z))
  SanityCheck(nloptr.wrapper(sample10, x = x, z = z))
  
})
