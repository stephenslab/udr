context("ud_fit")

test_that("Error is thrown for univariate case",{
  set.seed(1)
  X <- matrix(rnorm(20))
  expect_error(ud_init(X))
})

test_that("R and C++ implementations of ud_fit give the same result",{
 
  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X
   
  # Compare the R and C++ implementations with different invocations
  # of ud_init and ud_fit.
  fit0 <- ud_init(X)
  capture.output({
    fit1 <- ud_fit(fit0,control = list(maxiter = 20,version = "R"))
    fit2 <- ud_fit(fit0,control = list(maxiter = 20,version = "Rcpp"))
    fit3 <- ud_fit(fit0,control = list(unconstrained.update = "teem",
                                       maxiter = 20,version = "R"))
    fit4 <- ud_fit(fit0,control = list(unconstrained.update = "teem",
                                       maxiter = 20,version = "Rcpp"))
  })
  # fit0 <- ud_init(X,V = rep(list(dat$V),n))
  # capture.output({
  #   fit5 <- ud_fit(fit0,control = list(maxiter = 20,version = "R"))
  #   fit6 <- ud_fit(fit0,control = list(maxiter = 20,version = "Rcpp"))
  #  })
  
  # The likelihoods should be the same.
  expect_equal(fit1$progress$loglik,fit2$progress$loglik)
  expect_equal(fit3$progress$loglik,fit4$progress$loglik)
  
  # The parameter estimates should also be the same.
  expect_equal(fit1$w,fit2$w,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$w,fit4$w,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$U,fit2$U,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$U,fit4$U,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$V,fit2$V,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$V,fit4$V,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$progress[,-6],fit2$progress[,-6],scale = 1,tolerance=1e-12)
  expect_equal(fit3$progress[,-6],fit4$progress[,-6],scale = 1,tolerance=1e-12)
})
