context("Posterior likelihood estimation")

test_that("PosteriorMDP", {
  
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
  
  PosteriorMDP.wrapper <- function(sample, x, z)
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
    
    posterior <- PosteriorMDP(sample$y.sample, 
                                sample$y.lagged, x, z,
                                c(theta$alpha), c(theta$mu), c(theta$sigma),
                                rho, beta, c(gamma))
    posterior.initial.fixed <- PosteriorMDPInitialFixed(sample$y.sample, 
                                sample$y.lagged, x, z,
                                c(theta$alpha), c(theta$mu), c(theta$sigma),
                                rho, beta, c(gamma))
    posterior.R <- EstimateMDPPosteriors(theta = theta, 
                                             y = sample$y.sample, y.lagged = sample$y.lagged,
                                             x = x, z = z, initial.fixed = FALSE)
    posterior.initial.fixed.R <- EstimateMDPPosteriors(theta = theta, 
                                 y = sample$y.sample, y.lagged = sample$y.lagged,
                                x = x, z = z, initial.fixed = TRUE)
    list(posterior = posterior, 
         posterior.initial.fixed = posterior.initial.fixed,
         posterior.R = posterior.R,
         posterior.initial.fixed.R = posterior.initial.fixed.R,
         sample = sample)
  }
  
  SanityCheck <- function(wrapped)
  {
    M <- length(wrapped$sample$MDP.model$theta$alpha)
    # check if posteriors are valid
    for (i in 1:N.fixed)
    {
      expect_equal(sum(wrapped$posterior[i,]), 1)
      expect_equal(sum(wrapped$posterior.initial.fixed[i,]), 1)
      for (j in 1:M)
      {
        expect_gt(wrapped$posterior[i,j], 0)
        expect_gt(wrapped$posterior.initial.fixed[i,j], 0)
        expect_equal(wrapped$posterior[i,j], wrapped$posterior.R[i,j])
        expect_equal(wrapped$posterior.initial.fixed[i,j], 
                     wrapped$posterior.initial.fixed.R[i,j])
      }
    }
    
  }
  
  SanityCheck(PosteriorMDP.wrapper(sample1, x.empty, z.empty))
  SanityCheck(PosteriorMDP.wrapper(sample2, x.empty, z.empty))
  SanityCheck(PosteriorMDP.wrapper(sample3, x.empty, z.empty))
  SanityCheck(PosteriorMDP.wrapper(sample4, x.empty, z.empty))
  SanityCheck(PosteriorMDP.wrapper(sample5, x, z.empty))
  SanityCheck(PosteriorMDP.wrapper(sample6, x, z.empty))
  SanityCheck(PosteriorMDP.wrapper(sample7, x.empty, z))
  SanityCheck(PosteriorMDP.wrapper(sample8, x.empty, z))
  SanityCheck(PosteriorMDP.wrapper(sample9, x, z))
  SanityCheck(PosteriorMDP.wrapper(sample10, x, z))
  
})
