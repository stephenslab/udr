context("mvebnm")

test_that("A few basic tests of mvebnm",{
 
  # Simulate data.
  set.seed(1)
  dat <- simulate_mvebnm_data_2d(100)

  # First, compare the R and C++ implementations of the Bovy et al
  # (2013) EM algorithm.
  capture.output(
    fit1 <- mvebnm(dat$X,U = dat$U,
                   control = list(version = "R",update.S = "em",
                                  update.U = "em",maxiter = 20)))
  capture.output(
    fit2 <- mvebnm(dat$X,U = dat$U,
                   control = list(version = "Rcpp",update.S = "em",
                       update.U = "em",maxiter = 20)))

  # The likelihood should be monotonically increasing at each
  # iteration.
  expect_nondecreasing(fit1$progress$loglik)
  expect_nondecreasing(fit2$progress$loglik)
  
  # The parameter estimates and likelihoods produced by the 
  # R and C++ implementations should be the same.
  expect_equal(fit1$w,fit2$w,scale = 1,tolerance = 1e-15)
  expect_equal(fit1$U,fit2$U,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$S,fit2$S,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$loglik,fit2$loglik,scale = 1,tolerance = 1e-12)
  expect_equal(fit1$progress[,-6],fit2$progress[,-6],
               scale = 1,tolerance = 1e-12)

  # Check that the dimension names are assigned correctly.
  expect_equal(attributes(fit1$w),attributes(dat$w))
  expect_equal(attributes(fit2$w),attributes(dat$w))
  expect_equal(attributes(fit1$S),attributes(dat$S))
  expect_equal(attributes(fit2$S),attributes(dat$S))
  expect_equal(attributes(simplify2array(fit1$U)),
               attributes(simplify2array(dat$U)))
  expect_equal(attributes(simplify2array(fit2$U)),
               attributes(simplify2array(dat$U)))
})

test_that("A few basic tests of mvebnm for univariate (m = 1) data",{
  set.seed(1)

})
