context("fastr_fit")

dt <- runif(1, 0.001, 0.01)
n_bin <- sample(10:2000, 1)
n_cell <- sample(2:10, 1)
k <- runif(n_cell, 0.1, 0.4)
alpha <- runif(n_cell, 2, 5)
n_trial <- sample(1:5, 1)
Sig <- cov2cor(rWishart(1, n_cell+1, diag(n_cell))[,,1])
rho <- Sig[lower.tri(Sig)]
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, rho=rho)
Y <- sim$Y
x <- sim$x
n_factor <- ifelse(n_cell >= 5, 2, 1)
init_param <- list(log_k = rep(-2, n_cell),
		   log_a = rep(0, n_cell),
		   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
		   x = prop_paths(Y, dt, rep(-2, n_cell), rep(0, n_cell)))

test_that("`fastr_fit()` works without bugs for method='joint'",{
  fit <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, init=init_param,
		   method="joint", woodbury=F, silent=T)
  success <- TRUE
  expect_true(success)
})

test_that("`fastr_fit()` works without bugs for method='2-step'",{
  fit <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, init=init_param,
		   method="2step", woodbury=F, silent=T)
  success <- TRUE
  expect_true(success)
})

test_that("`fastr_fit()` works without bugs using Woodbury formula",{
  fit <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, init=init_param,
		   method="2step", woodbury=T, silent=T)
  success <- TRUE
  expect_true(success)
})
