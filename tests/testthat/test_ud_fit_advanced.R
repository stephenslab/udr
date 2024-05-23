context("ud_fit_advanced")

test_that(paste("R and C++ implementations of advanced model fitting ",
                "interface produce the same result for i.i.d. case"),{

  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Perform 4 EM updates.
  fit <- ud_init(X)
  control <- list(scaled.update = "none",rank1.update = "none",
                  unconstrained.update = "none")
  capture.output(fit <- ud_fit(fit,control = c(control,list(maxiter = 4))))

  # Update the responsibilities matrix.
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

  # Update the prior covariance matrices, U, using the TED algorithm.
  control <- list(scaled.update = "none",rank1.update = "ted",
                  unconstrained.update = "ted")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  updates1 <- assign_prior_covariance_updates(fit1,control1)$covupdates
  updates2 <- assign_prior_covariance_updates(fit2,control2)$covupdates
  fit1 <- update_prior_covariances(fit1,updates1)
  fit2 <- update_prior_covariances(fit2,updates2)
  removeQ <- function (U)
    lapply(U,function (x) { x["Q"] <- NULL; return(x) })
  fit1$U <- removeQ(fit1$U)
  fit2$U <- removeQ(fit2$U)
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-14)

  # Update the prior covariance matrices, U, using the ED updates.
  control <- list(scaled.update = "none",rank1.update = "none",
                  unconstrained.update = "ed")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  updates1 <- assign_prior_covariance_updates(fit1,control1)$covupdates
  updates2 <- assign_prior_covariance_updates(fit2,control2)$covupdates
  fit1   <- update_prior_covariances(fit1,updates1)
  fit2   <- update_prior_covariances(fit2,updates2)
  fit1$U <- removeQ(fit1$U)
  fit2$U <- removeQ(fit2$U)
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-14)
})

test_that(paste("R and C++ implementations of advanced model fitting ",
                "interface produce the same result for non-i.i.d. case"),{

  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Vary the measurement error slightly for each observation.
  V <- vector("list",n)
  for (i in 1:n)
    V[[i]] <- dat$V + 0.01 * sim_unconstrained(2)
  
  # Perform 4 EM updates.
  fit <- ud_init(X,V = V)
  control <- list(scaled.update = "none",rank1.update = "none",
                  unconstrained.update = "none")
  capture.output(fit <- ud_fit(fit,control = c(control,list(maxiter = 4))))

  # Update the responsibilities matrix.
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

  # Update the prior covariance matrices, U, using the ED updates.
  control <- list(scaled.update = "none",rank1.update = "none",
                  unconstrained.update = "ed")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  updates1 <- assign_prior_covariance_updates(fit1,control1)$covupdates
  updates2 <- assign_prior_covariance_updates(fit2,control2)$covupdates
  fit1 <- update_prior_covariances(fit1,updates1)
  fit2 <- update_prior_covariances(fit2,updates2)
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-8)
  
  # Attempting to perform the TED updates should result in an error.
  control <- list(scaled.update = "none",rank1.update = "ted",
                  unconstrained.update = "ted")
  control1 <- control
  control2 <- control
  control1$version <- "R"
  control2$version <- "Rcpp"
  updates1 <- assign_prior_covariance_updates(fit1,control1)$covupdates
  updates2 <- assign_prior_covariance_updates(fit2,control2)$covupdates
  expect_error(update_prior_covariances(fit1,updates1))
  expect_error(update_prior_covariances(fit2,updates2))
})
