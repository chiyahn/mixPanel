context("Lagging columns")

test_that("GetLaggedColumn", {
  expect_identical(GetLaggedColumn(0, 1:5, s = 3), 4:5)
  expect_identical(GetLaggedColumn(1, 1:5, s = 3), 3:4)
  expect_identical(GetLaggedColumn(2, 1:5, s = 3), 2:3)
  expect_identical(GetLaggedColumn(3, 1:5, s = 3), 1:2)

  expect_identical(GetLaggedColumn(0, 1:5, s = 2), 3:5)
  expect_identical(GetLaggedColumn(1, 1:5, s = 2), 2:4)
  expect_identical(GetLaggedColumn(2, 1:5, s = 2), 1:3)

})

test_that("GetLaggedAndSample", {
  lag.and.sample <- GetLaggedAndSample(y = 1:5, s = 3)
  expect_equal(lag.and.sample$y.sample, as.matrix(4:5))
  expect_equal((lag.and.sample$y.lagged)[,1], 3:4)
  expect_equal((lag.and.sample$y.lagged)[,2], 2:3)
  expect_equal((lag.and.sample$y.lagged)[,3], 1:2)
  
  lag.and.sample <- GetLaggedAndSample(y = 1:5, s = 0)
  expect_equal(lag.and.sample$y.sample, as.matrix(1:5))
  expect_equal((lag.and.sample$y.lagged)[,1], rep(0,5))
})
