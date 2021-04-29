context("ud_fit")

test_that(paste("R and C++ implementations of ud_fit produce same result",
                "when V is a matrix or list (and all Vs are the same)"),{

  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Run the R and C++ implementations of ud_fit when V is a matrix or
  # a list, and V is not updated. The updates of the prior covariance
  # matrices U are skipped in this test.
  control  <- list(maxiter = 20,resid.update = "none",scaled.update = "none",
                   rank1.update = "none",unconstrained.update = "none")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  set.seed(1); fit1 <- ud_init(X,V = dat$V,control = control1)
  set.seed(1); fit2 <- ud_init(X,V = dat$V,control = control2)
  set.seed(1); fit3 <- ud_init(X,V = rep(list(dat$V),n),control = control1)
  set.seed(1); fit4 <- ud_init(X,V = rep(list(dat$V),n),control = control2)
  capture.output(fit1 <- ud_fit(fit1,control = control1))
  capture.output(fit2 <- ud_fit(fit2,control = control2))
  capture.output(fit3 <- ud_fit(fit3,control = control1))
  capture.output(fit4 <- ud_fit(fit4,control = control2))
  
  # The likelihoods should be non-decreasing.
  expect_nondecreasing(fit1$progress$loglik)
  
  # The outputs of all four variants of ud_fit should be the same
  # (except for V and the timings).
  fit1$progress$timing <- 0
  fit2$progress$timing <- 0
  fit3$progress$timing <- 0
  fit4$progress$timing <- 0
  fit1["V"] <- NULL
  fit2["V"] <- NULL
  fit3["V"] <- NULL
  fit4["V"] <- NULL
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-12)
  expect_equal(fit1,fit3,scale = 1,tolerance = 1e-12)
  expect_equal(fit1,fit4,scale = 1,tolerance = 1e-12)

  # Now compare the R and C++ implementations when V is updated.
  control1$resid.update = "em"
  control2$resid.update = "em"
  set.seed(1); fit1 <- ud_init(X,V = dat$V,control = control1)
  set.seed(1); fit2 <- ud_init(X,V = dat$V,control = control2)
  capture.output(fit1 <- ud_fit(fit1,control = control1))
  capture.output(fit2 <- ud_fit(fit2,control = control2))
  expect_nondecreasing(fit1$progress$loglik)

  # The ud_fit outputs should be the same (except for the timings),
  # and the likelihoods should be non-decreasing.
  fit1$progress$timing <- 0
  fit2$progress$timing <- 0
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-12)
})

test_that(paste("Check R and C++ implementations of prior covariance",
                "matrix (U) updates when V is a matrix"),{
  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Run ud_fit with unconstrained.update = "ed".
  control  <- list(maxiter = 20,resid.update = "em",scaled.update = "none",
                   rank1.update = "none",unconstrained.update = "ed")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  set.seed(1); fit1 <- ud_init(X,V = dat$V,control = control1)
  set.seed(1); fit2 <- ud_init(X,V = dat$V,control = control2)
  capture.output(fit1 <- ud_fit(fit1,control = control1))
  capture.output(fit2 <- ud_fit(fit2,control = control2))

  # The likelihoods should be non-decreasing, and both ud_fit outputs
  # should be the same (except for the timings).
  fit1$progress$timing <- 0
  fit2$progress$timing <- 0
  expect_nondecreasing(fit1$progress$loglik)
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-12)

  # Run ud_fit with unconstrained.update = "teem".
  control  <- list(maxiter = 20,resid.update = "em",scaled.update = "none",
                   rank1.update = "none",unconstrained.update = "teem")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  set.seed(1); fit1 <- ud_init(X,V = dat$V,control = control1)
  set.seed(1); fit2 <- ud_init(X,V = dat$V,control = control2)
  capture.output(fit1 <- ud_fit(fit1,control = control1))
  # capture.output(fit2 <- ud_fit(fit2,control = control2))
  
  # The likelihoods should be non-decreasing, and both ud_fit outputs
  # should be the same (except for the timings).
  fit1$progress$timing <- 0
  # fit2$progress$timing <- 0
  expect_nondecreasing(fit1$progress$loglik)
  # expect_equal(fit1,fit2,scale = 1,tolerance = 1e-12)
})

test_that(paste("Check R and C++ implementations of prior covariance",
                "matrix (U) updates when the Vs are not all the same"),{
  # TO DO.
})
