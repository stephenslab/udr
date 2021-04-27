context("ud_fit")

test_that(paste("R and C++ implementations of ud_fit produce same result",
                "when V is a matrix or list"),{

  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Run the R and C++ implementations of ud_fit when V is a matrix or
  # a list, and V is not updated.
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

  # The outputs of all four variants of ud_fit should be the same
  # (except for the timings).
  fit1$progress$timing <- 0
  fit2$progress$timing <- 0
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-12)
})

# test_that("R and C++ implementations of ud_fit produce same result",{
 
  # Compare the R and C++ implementations with different invocations
  # of ud_init and ud_fit.
  #
  # fit0 <- ud_init(X)
  # fit1 <- ud_fit(fit0,control = list(maxiter = 20,version = "R"))
  # fit2 <- ud_fit(fit0,control = list(maxiter = 20,version = "Rcpp"))
  # fit3 <- ud_fit(fit0,control = list(unconstrained.update = "teem",
  #                                    maxiter = 20,version = "R"))
  # fit4 <- ud_fit(fit0,control = list(unconstrained.update = "teem",
  #                                    maxiter = 20,version = "Rcpp"))
  
  # The likelihoods should be the same.
  # expect_equal(fit1$progress$loglik,fit2$progress$loglik)
  # expect_equal(fit3$progress$loglik,fit4$progress$loglik)
  
  # The parameter estimates should also be the same.
  # expect_equal(fit1$w,fit2$w,scale = 1,tolerance = 1e-14)
  # expect_equal(fit3$w,fit4$w,scale = 1,tolerance = 1e-14)
  # expect_equal(fit1$U,fit2$U,scale = 1,tolerance = 1e-14)
  # expect_equal(fit3$U,fit4$U,scale = 1,tolerance = 1e-14)
  # expect_equal(fit1$V,fit2$V,scale = 1,tolerance = 1e-14)
  # expect_equal(fit3$V,fit4$V,scale = 1,tolerance = 1e-14)
  # expect_equal(fit1$progress[,-6],fit2$progress[,-6],scale = 1,tolerance=1e-12)
  # expect_equal(fit3$progress[,-6],fit4$progress[,-6],scale = 1,tolerance=1e-12)
# })
