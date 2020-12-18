context("ud_fit")

test_that("A few basic tests of ud_fit",{
 
  # Simulate data.
  set.seed(1)
  dat <- simulate_ud_data_2d(100)

  # First, compare the R and C++ implementations of the Bovy et al
  # (2013) EM algorithm.
  capture.output(
    fit1 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "R",update.S = "em",
                                  update.U = "ed",maxiter = 20)))
  capture.output(
    fit2 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "Rcpp",update.S = "em",
                                  update.U = "ed",maxiter = 20)))

  # The likelihood should increase at each iteration.
  expect_nondecreasing(fit1$progress$loglik)
  expect_nondecreasing(fit2$progress$loglik)
  
  # The parameter estimates and likelihoods produced by the 
  # R and C++ implementations should be the same.
  expect_equal(fit1$w,fit2$w,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$U,fit2$U,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$V,fit2$V,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$loglik,fit2$loglik,scale = 1,tolerance = 1e-12)
  expect_equal(fit1$progress[,-6],fit2$progress[,-6],
               scale = 1,tolerance = 1e-12)

  # Check that the row and column names are assigned correctly.
  expect_equal(attributes(fit1$w),attributes(dat$w))
  expect_equal(attributes(fit2$w),attributes(dat$w))
  expect_equal(attributes(fit1$V),attributes(dat$V))
  expect_equal(attributes(fit2$V),attributes(dat$V))
  expect_equal(attributes(simplify2array(fit1$U)),
               attributes(simplify2array(dat$U)))
  expect_equal(attributes(simplify2array(fit2$U)),
               attributes(simplify2array(dat$U)))

  # Next, compare the R and C++ implementations of "truncated
  # eigenvalue" EM updates.
  capture.output(
    fit3 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "R",update.S = "em",
                                  update.U = "teem",maxiter = 20)))
  capture.output(
    fit4 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "Rcpp",update.S = "em",
                                  update.U = "teem",maxiter = 20)))

  # The likelihood should increase at each iteration.
  expect_nondecreasing(fit3$progress$loglik)
  expect_nondecreasing(fit4$progress$loglik)

  # The parameter estimates and likelihoods produced by the 
  # R and C++ implementations should be the same.
  expect_equal(fit3$w,fit4$w,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$U,fit4$U,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$V,fit4$V,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$loglik,fit4$loglik,scale = 1,tolerance = 1e-12)
  expect_equal(fit3$progress[,-6],fit4$progress[,-6],
               scale = 1,tolerance = 1e-12)

  # Check that the row and column names are assigned correctly.
  expect_equal(attributes(fit3$w),attributes(dat$w))
  expect_equal(attributes(fit4$w),attributes(dat$w))
  expect_equal(attributes(fit3$V),attributes(dat$V))
  expect_equal(attributes(fit4$V),attributes(dat$V))
  expect_equal(attributes(simplify2array(fit3$U)),
               attributes(simplify2array(dat$U)))
  expect_equal(attributes(simplify2array(fit4$U)),
               attributes(simplify2array(dat$U)))
})

test_that("A few basic tests of ud_fit for univariate (m = 1) data",{

  # Simulate univariate data.
  set.seed(1)
  dat <- simulate_ebnm_data(100)

  # Compare the R and C++ implementations of the Bovy et al (2013) EM
  # algorithm.
  capture.output(
    fit1 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "R",update.S = "em",
                                  update.U = "ed",maxiter = 20)))
  capture.output(
    fit2 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "Rcpp",update.S = "em",
                                  update.U = "ed",maxiter = 20)))

  # The likelihood should increase at each iteration.
  expect_nondecreasing(fit1$progress$loglik)
  expect_nondecreasing(fit2$progress$loglik)
  
  # The parameter estimates and likelihoods produced by the 
  # R and C++ implementations should be the same.
  expect_equal(fit1$w,fit2$w,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$U,fit2$U,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$V,fit2$V,scale = 1,tolerance = 1e-14)
  expect_equal(fit1$loglik,fit2$loglik,scale = 1,tolerance = 1e-12)
  expect_equal(fit1$progress[,-6],fit2$progress[,-6],
               scale = 1,tolerance = 1e-12)

  # Check that the names are assigned correctly, and that the
  # variances are scalars.
  expect_equal(attributes(fit1$w),attributes(dat$w))
  expect_equal(attributes(fit2$w),attributes(dat$w))
  expect_equal(attributes(unlist(fit1$U)),attributes(unlist(dat$U)))
  expect_equal(attributes(unlist(fit2$U)),attributes(unlist(dat$U)))
  expect_false(is.matrix(fit1$V))
  expect_false(is.matrix(fit2$V))
  expect_false(any(sapply(fit1$U,is.matrix)))
  expect_false(any(sapply(fit2$U,is.matrix)))

  # Compare the R and C++ implementations of "truncated eigenvalue" EM
  # updates.
  capture.output(
    fit3 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "R",update.S = "em",
                                  update.U = "teem",maxiter = 20)))
  capture.output(
    fit4 <- ud_fit(dat$X,U = dat$U,
                   control = list(version = "Rcpp",update.S = "em",
                                  update.U = "teem",maxiter = 20)))

  # The likelihood should increase at each iteration.
  expect_nondecreasing(fit3$progress$loglik)
  expect_nondecreasing(fit4$progress$loglik)

  # The parameter estimates and likelihoods produced by the 
  # R and C++ implementations should be the same.
  expect_equal(fit3$w,fit4$w,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$U,fit4$U,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$V,fit4$V,scale = 1,tolerance = 1e-14)
  expect_equal(fit3$loglik,fit4$loglik,scale = 1,tolerance = 1e-12)
  expect_equal(fit3$progress[,-6],fit4$progress[,-6],
               scale = 1,tolerance = 1e-12)

  # Check that the names are assigned correctly, and that the
  # variances are scalars.
  expect_equal(attributes(fit3$w),attributes(dat$w))
  expect_equal(attributes(fit4$w),attributes(dat$w))
  expect_equal(attributes(unlist(fit3$U)),attributes(unlist(dat$U)))
  expect_equal(attributes(unlist(fit4$U)),attributes(unlist(dat$U)))
  expect_false(is.matrix(fit3$V))
  expect_false(is.matrix(fit4$V))
  expect_false(any(sapply(fit3$U,is.matrix)))
  expect_false(any(sapply(fit4$U,is.matrix)))
})
