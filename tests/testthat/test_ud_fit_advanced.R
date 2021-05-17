context("ud_fit_advanced")

test_that(paste("R and C++ implementations of advanced model fitting ",
                "interface produce the same resultud_fit produce same",
                "result"),{

  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Perform 4 EM updates.
  fit <- ud_init(X)
  capture.output(fit <- ud_fit(fit1,control = list(maxiter = 4)))

  # Update the responsibilities matrix.
  control <- list(scaled.update = "none",rank1.update = "none")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  fit1 <- compute_posterior_probs(fit,version = "R")
  fit2 <- compute_posterior_probs(fit,version = "Rcpp")
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-15)

  # Update the mixture weights.
  fit1 <- update_mixture_weights(fit1)
  fit2 <- update_mixture_weights(fit2)
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-15)
  
  # Update the residual covariance, V.
  fit1 <- update_resid_covariance(fit1,version = "R")
  fit2 <- update_resid_covariance(fit2,version = "Rcpp")
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-15)

  # Update the prior covariance matrices, U.
  updates1 <- assign_prior_covariance_updates(fit1,control1)$covupdates
  updates2 <- assign_prior_covariance_updates(fit2,control2)$covupdates
  fit1 <- update_prior_covariances(fit1,updates1)
  fit2 <- update_prior_covariances(fit2,updates2)
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-15)
})

