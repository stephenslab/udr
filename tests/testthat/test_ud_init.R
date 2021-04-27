context("ud_init")

test_that(paste("Error is thrown when X has less than 2 rows or columns, or",
                "when less than 2 prior covariance matrices are requested"),{
  set.seed(1)
  X <- matrix(rnorm(20))
  expect_error(ud_init(X))
  expect_error(ud_init(t(X)))
  X <- matrix(rnorm(20),10,2)
  expect_error(ud_init(X,U_scaled = NULL,n_rank1 = 0,n_unconstrained = 1))
})

test_that("ud_init produces same result when V is a matrix or list",{

  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Run ud_fit when V is a matrix and a list.
  set.seed(1); fit1 <- ud_init(X,V = dat$V)
  set.seed(1); fit2 <- ud_init(X,V = rep(list(dat$V),n))

  # The two ud_fit objects should be the same, except for V.
  fit1["V"] <- NULL
  fit2["V"] <- NULL
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-15)
})

test_that("R and C++ implementations of ud_init produce same result",{

  # Simulate data.
  set.seed(1)
  n   <- 100
  dat <- simulate_ud_data_2d(n)
  X   <- dat$X

  # Run ud_fit when V is a matrix. The two outputs should be the same.
  set.seed(1); fit1 <- ud_init(X,V = dat$V,control = list(version = "R"))
  set.seed(1); fit2 <- ud_init(X,V = dat$V,control = list(version = "Rcpp"))
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-12)

  # Run ud_fit when V is a list. The two outputs should be the same.
  V <- vector("list",n)
  for (i in 1:n)
    V[[i]] <- sim_unconstrained(2) + diag(2)
  set.seed(1); fit1 <- ud_init(X,V,control = list(version = "R"))
  set.seed(1); fit2 <- ud_init(X,V,control = list(version = "Rcpp"))
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-12)
})
