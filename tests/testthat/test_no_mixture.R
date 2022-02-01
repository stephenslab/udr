context("no_mixture")

test_that("ud_fit works as expected when k = 1",{

  # Simulate data from an Ultimate Deconvolution model in which the
  # prior consists of a single mixture component.
  set.seed(1)
  n <- 100
  U <- rbind(c(1,0.5),
             c(0.5,1))
  X <- simulate_ud_data(n,U = U)

  # Fit an Ultimate Deconvolution model with a single unconstrained U.
  fit0 <- ud_init(X,V = V,n_unconstrained = 1,n_rank1 = 0,U_scaled = NULL)
  capture.output(
    fit1 <- ud_fit(fit0,control = list(unconstrained.update = "ted",
                                       version = "R")))
  capture.output(
    fit2 <- ud_fit(fit0,control = list(unconstrained.update = "ted",
                                       version = "Rcpp")))

  # The likelihoods should be non-decreasing.
  expect_nondecreasing(fit1$progress$loglik)
  expect_nondecreasing(fit2$progress$loglik)

  # The EM should automatically terminate after two iterations since
  # we are only updating a single prior covariance matrix, U.
  expect_equal(fit1$progress$iter,1:2)
  expect_equal(fit2$progress$iter,1:2)
  
  # The outputs from both calls to ud_fit should be the same (except
  # for the timings).
  fit1$progress$timing <- 0
  fit2$progress$timing <- 0
  expect_equal(fit1,fit2)
})
