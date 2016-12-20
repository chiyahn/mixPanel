context("Fisher matrix estimation")

test_that("EstimateMDPFisher", {
  
  set.seed(1234)
  s1 <- 0
  s2 <- 3
  N.fixed <- 60
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
  
  EstimateMDPFisher.wrapper <- function(sample, x = NULL, z = NULL)
  {
    theta <- (sample$MDP.model)$theta

    fisher <- EstimateMDPFisher(theta = theta, 
                                y = sample$y.sample, y.lagged = sample$y.lagged,
                                x = x, z = z, initial.fixed = FALSE)
    fisher.fixed <- EstimateMDPFisher(theta = theta, 
                                y = sample$y.sample, y.lagged = sample$y.lagged,
                                x = x, z = z, initial.fixed = TRUE)
    list(fisher = fisher,
         fisher.fixed = fisher.fixed,
         sample = sample)
  }
  
  SanityCheck <- function(wrapped)
  {
    N <- ncol(wrapped$sample$y.sample)
    ses <- sqrt(diag(solve(wrapped$fisher)))/sqrt(N)
    ses.fixed <- sqrt(diag(solve(wrapped$fisher)))/sqrt(N)
    
    size <- length(ses)
    for (i in 1:size) # check if computed standard errors make sense
    {
      expect_gt(ses[size], 0)
      expect_gt(ses.fixed[size], 0)
    }
  }
  
  SanityCheck(EstimateMDPFisher.wrapper(sample1))
  SanityCheck(EstimateMDPFisher.wrapper(sample2))
  SanityCheck(EstimateMDPFisher.wrapper(sample3))
  SanityCheck(EstimateMDPFisher.wrapper(sample4))
  SanityCheck(EstimateMDPFisher.wrapper(sample5, x = x))
  SanityCheck(EstimateMDPFisher.wrapper(sample6, x = x))
  SanityCheck(EstimateMDPFisher.wrapper(sample7, z = z))
  SanityCheck(EstimateMDPFisher.wrapper(sample8, z = z))
  SanityCheck(EstimateMDPFisher.wrapper(sample9, x = x, z = z))
  SanityCheck(EstimateMDPFisher.wrapper(sample10, x = x, z = z))
  
})
