
test_that("unconstrained FA increases loglikelihood", {
  
  # Simulate data
  seed <- 1
  n <- 500
  U.true <- diag(2)
  V <- diag(2)
  X <- simulate_one_component(seed, n, U.true, V)

  # Initialization
  p <- runif(n)
  maxiter <- 10
  U <- sim_unconstrained(2)
  logliks <- rep(NA, maxiter)
  
  for (i in 1:maxiter){
    Unew <- update_prior_covariance_unconstrained_fa(X, U, p)
    U <- Unew
    logliks[i] <- loglik_weighted_single_component(X, U, V, p)
  }
  expect_nondecreasing(logliks)
})

test_that("unconstrained ED_iid increases loglikelihood", {
  # Simulate data
  seed <- 1
  n <- 500
  U.true <- diag(2)
  V <- diag(2)
  X <- simulate_one_component(seed, n, U.true, V)

  # Initialization
  p <- runif(n)
  maxiter <- 10
  U <- sim_unconstrained(2)
  logliks <- rep(NA, maxiter)
  
  for (i in 1:maxiter){
    Unew <-  update_prior_covariance_ed_iid(X, U, V, p)
    U <- Unew
    logliks[i] <- loglik_weighted_single_component(X, U, V, p)
  }
  expect_nondecreasing(logliks)
})

test_that("scaled FA_general increases loglikelihood", {
  
  # Simulate data
  seed <- 1
  n <- 500
  U0 <- matrix(c(1,0,0,0), ncol = 2, nrow = 2)
  V.g <- simplify2array(replicate(n, diag(2), simplify=FALSE))
  X = simulate_one_component(seed, n, U0, V)
  
  # Initialization
  maxiter <- 10
  U <- create_scaled_matrix_struct(X, U0)
  r <- sum(eigen(U$U0)$values > 1e-15)
  minval <- 0
  p <- runif(n)
  logliks <- rep(NA, maxiter)
  
  for (i in 1: maxiter){
    s <- update_prior_covariance_scaled_fa_general(X, U$U0, V.g, p, U$s, r)
    Unew <- update_prior_covariance_scaled_struct(U, s)
    logliks[i] <- loglik_weighted_single_component(X, U$mat, V.g, p)
    U <- Unew
}
  expect_nondecreasing(logliks)
})

test_that("rank1 FA_iid increases loglikelihood", {
  
  # Simulate data
  seed <- 1
  n <- 500
  U.true <- matrix(c(1,0,0,0), ncol = 2, nrow = 2)
  V <- diag(2)
  X <- simulate_one_component(seed, n, U.true, V)
  
  # Initialization
  maxiter <- 10
  U <- sim_rank1(2)
  U <- create_rank1_matrix_struct(X, U)
  p <- runif(n)
  logliks <- rep(NA, maxiter)
  
  for (i in 1: maxiter){
    vec <- update_prior_covariance_rank1_fa_iid(X, U$vec, p)
    Unew <- update_prior_covariance_rank1_struct(U, vec)
    logliks[i] <- loglik_weighted_single_component(X, U$mat, V, p)
    U <- Unew
}
  expect_nondecreasing(logliks)
})


test_that("unconstrained ED_general vs. ED_iid", {
  
  set.seed(1)
  # Simulate data
  n <- 500
  U.true <- diag(2)
  V <- diag(2)
  V.g <- replicate(n, V, simplify = "array")
  X <- simulate_one_component(1, n,  U.true, V)
  
  # Initialization
  p <- runif(n)
  maxiter <- 10
  U.iid <- U.general <- sim_unconstrained(2)
  
  logliks_iid = rep(NA, maxiter)
  logliks_general = rep(NA, maxiter)
  
  for (i in 1:maxiter){
    Unew.iid <- update_prior_covariance_ed_iid(X, U.iid, V, p)
    Unew.general <- update_prior_covariance_ed_general(X, U.general, V.g, p)
    
    U.iid <- Unew.iid
    U.general <- Unew.general
    
    logliks_iid[i] <- loglik_weighted_single_component(X, U.iid, V, p)
    logliks_general[i] <- loglik_weighted_single_component(X, U.general, V, p)
  }
  expect_equal(logliks_general, logliks_iid)
  expect_equal(U.iid, U.general)
})


test_that("rank1 fa_general vs. fa_iid", {
  set.seed(1)
  # Simulate data
  n <- 500
  U.true <- matrix(c(1,0,0,0), ncol = 2, nrow = 2)
  V <- diag(2)
  V.g <- replicate(n, V, simplify = "array")
  X <- simulate_one_component(1, n,  U.true, V)
  
  # Initialization
  p <- runif(n)
  maxiter <- 10
  U.iid <- U.general <- sim_rank1(2)
  logliks_iid <- rep(NA, maxiter)
  logliks_general <- rep(NA, maxiter)
  U.iid <- create_rank1_matrix_struct(X, U.iid)
  U.general <- create_rank1_matrix_struct(X, U.general)
  
  for (i in 1:maxiter){
    
    vec1 <- update_prior_covariance_rank1_fa_iid(X, U.iid$vec, p)
    Unew.iid <- update_prior_covariance_rank1_struct(U.iid, vec1)
    
    vec2 <- update_prior_covariance_rank1_fa_general(X, U.general$vec, V.g, p)
    Unew.general <- update_prior_covariance_rank1_struct(U.general, vec2)
    
    U.iid <- Unew.iid
    U.general <- Unew.general
    
    logliks_iid[i] <- loglik_weighted_single_component(X, U.iid$mat, V, p)
    logliks_general[i] <- loglik_weighted_single_component(X, U.general$mat, V.g, p)
  }
  expect_equal(logliks_general, logliks_iid)
  expect_equal(U.iid, U.general)
})

