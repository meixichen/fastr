context("joint_factor_model")
test_that("Joint optimization of the factor model converges and gives PD Hessian.",{
  dt <- 0.005
  n_bin <- 2000
  k <- c(0.15, 0.3, 0.2, 0.27, 0.25)
  alpha <- c(2.5, 4, 3.5, 5, 3)
  n_cell <- length(alpha)
  n_trial <- 5
  l1 <- c(0.8, 0.08, 0.75, 0.1, 0.8)
  l2 <- c(0, 0.8, 0.05, 0.85, 0.1)
  L <- cbind(l1, l2)
  sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
  Y <- sim$Y
  x <- sim$x
  n_factor <- 2
  init_param <- list(log_k = rep(-2, n_cell),
                     log_a = rep(0, n_cell),
                     Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                     x = prop_paths(Y, dt, rep(-2, n_cell), rep(0, n_cell)))

  # Efficient matrix inversion not applied, joint est
  fit_joint <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, woodbury=FALSE,
                         method="joint", init=init_param, silent=TRUE,
                         simplified=FALSE, lam=1, nu=5.)
  # Expect convergence indicator is 0, i.e., nlminb has converged
  expect_equal(fit_joint$env$nlminb_fit$convergence, 0)
  # Expect the Hessian matrix is positive definite, i.e., no NA standard errors
  expect_true(fit_joint$env$tmb_report$pdHess)

  # Efficient matrix inversion applied, joint est
  fit_eff_joint <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, woodbury=TRUE,
                             method="joint", init=init_param, silent=TRUE,
                             simplified=FALSE, lam=1, nu=5.)
  expect_equal(fit_eff_joint$env$nlminb_fit$convergence, 0)
  expect_true(fit_eff_joint$env$tmb_report$pdHess)
})

