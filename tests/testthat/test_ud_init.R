context("ud_init")

test_that("Error is thrown when X has less than 2 rows or columns",{
  set.seed(1)
  X <- matrix(rnorm(20))
  expect_error(ud_init(X))
  expect_error(ud_init(t(X)))
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
  expect_equal(fit1,fit2)
})
