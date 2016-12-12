context("Sample generation")

test_that("GenerateMDPSample", {

  s1 <- 0
  s2 <- 2
  T.fixed <- 5
  theta1 <- GenerateMDPTheta(M = 2, s = s1)
  theta2 <- GenerateMDPTheta(M = 2, s = s2)
  
  sample1 <- GenerateMDPSample(theta1, T = T.fixed)
  sample2 <- GenerateMDPSample(theta2, T = T.fixed)
  
  expect_equal(nrow(sample1$y), T.fixed)
  expect_equal(nrow(sample1$y.sample), T.fixed)
  expect_equal(nrow(sample1$y.lagged), T.fixed)
  
  expect_equal(nrow(sample2$y), T.fixed + s2)
  expect_equal(nrow(sample2$y.sample), T.fixed)
  expect_equal(nrow(sample2$y.lagged), T.fixed)
  
})
