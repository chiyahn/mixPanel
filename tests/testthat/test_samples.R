context("Sample generation")

test_that("GenerateMDPSample", {

  s1 <- 0
  s2 <- 3
  N.fixed <- 30
  T.fixed <- 20
  p <- 2
  q <- 3
  x <- matrix(runif(T.fixed * N.fixed * p), nrow = T.fixed)
  z <- matrix(runif(T.fixed * N.fixed * q, 1, 3), nrow = T.fixed)

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
  
  expect_equal(nrow(sample1$y), T.fixed)
  expect_equal(nrow(sample1$y.sample), T.fixed)
  expect_equal(nrow(sample1$y.lagged), T.fixed)
  
  expect_equal(nrow(sample2$y), T.fixed + s2)
  expect_equal(nrow(sample2$y.sample), T.fixed)
  expect_equal(nrow(sample2$y.lagged), T.fixed)
  
  expect_equal(nrow(sample3$y), T.fixed)
  expect_equal(nrow(sample3$y.sample), T.fixed)
  expect_equal(nrow(sample3$y.lagged), T.fixed)
  
  expect_equal(nrow(sample4$y), T.fixed + s2)
  expect_equal(nrow(sample4$y.sample), T.fixed)
  expect_equal(nrow(sample4$y.lagged), T.fixed)
  
  expect_equal(nrow(sample5$y), T.fixed)
  expect_equal(nrow(sample5$y.sample), T.fixed)
  expect_equal(nrow(sample5$y.lagged), T.fixed)
  expect_equal(nrow(sample5$x), T.fixed)
  expect_equal(ncol(sample5$x), (p*N.fixed))
  
  expect_equal(nrow(sample6$y), T.fixed + s2)
  expect_equal(nrow(sample6$y.sample), T.fixed)
  expect_equal(nrow(sample6$y.lagged), T.fixed)
  expect_equal(nrow(sample6$x), T.fixed)
  expect_equal(ncol(sample6$x), (p*N.fixed))
  
  expect_equal(nrow(sample7$y), T.fixed)
  expect_equal(nrow(sample7$y.sample), T.fixed)
  expect_equal(nrow(sample7$y.lagged), T.fixed)
  expect_equal(nrow(sample7$z), T.fixed)
  expect_equal(ncol(sample7$z), (q*N.fixed))
  
  expect_equal(nrow(sample8$y), T.fixed + s2)
  expect_equal(nrow(sample8$y.sample), T.fixed)
  expect_equal(nrow(sample8$y.lagged), T.fixed)
  expect_equal(nrow(sample8$z), T.fixed)
  expect_equal(ncol(sample8$z), (q*N.fixed))
  
  expect_equal(nrow(sample9$y), T.fixed)
  expect_equal(nrow(sample9$y.sample), T.fixed)
  expect_equal(nrow(sample9$y.lagged), T.fixed)
  expect_equal(nrow(sample9$x), T.fixed)
  expect_equal(ncol(sample9$x), (p*N.fixed))
  expect_equal(nrow(sample9$z), T.fixed)
  expect_equal(ncol(sample9$z), (q*N.fixed))
  
  expect_equal(nrow(sample10$y), T.fixed + s2)
  expect_equal(nrow(sample10$y.sample), T.fixed)
  expect_equal(nrow(sample10$y.lagged), T.fixed)
  expect_equal(ncol(sample10$y.lagged), (N.fixed * s2))
  expect_equal(nrow(sample10$x), T.fixed)
  expect_equal(ncol(sample10$x), (p*N.fixed))
  expect_equal(nrow(sample10$z), T.fixed)
  expect_equal(ncol(sample10$z), (q*N.fixed))
  
  
  
})
