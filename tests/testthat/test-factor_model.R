context("factor_model")

test_that("Factor models in C++ has same nll as the nll computed in R. ",{
  n_test <- 10
  for (ii in 1:n_test){
    #---- Data simulation ---------
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
    data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=1, nu=5.)

    #---- R nll --------------
    nll_r <- compute_Rnll(data, init_param)

    #---- TMB model nll -----------
    # Efficient matrix inversion not applied
    adfun <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, lam=1, nu=5., init=init_param,
		       woodbury=FALSE, method="joint", silent=TRUE, adfun_only=TRUE,
		       ignore_random=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param)) # negative log-likelihood computed by TMB/C++
    expect_equal(nll_tmb, nll_r)

    # Efficient matrix inversion applied
    adfun_eff <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, lam=1, nu=5., init=init_param,
			   woodbury=TRUE, method="joint", silent=TRUE, adfun_only=TRUE,
			   ignore_random=TRUE)
    nll_tmb_eff <- adfun_eff$fn(unlist(init_param))
    expect_equal(nll_tmb_eff, nll_r)
  }
})

test_that("2-step optimization of the factor model converges and gives PD Hessian.",{
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
  fit_eff_2step <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, woodbury=TRUE,
			     method="2step", init=init_param, silent=TRUE,
			     simplified=FALSE, lam=1, nu=5.)
  expect_equal(fit_eff_2step$env$nlminb_fit$convergence, 0)
  expect_true(fit_eff_2step$env$tmb_report$pdHess)
})


test_that("Factor model works for a single neuron.",{
  k <- 0.4
  alpha <- 1.2
  dt <- 0.05
  n_bin <- 2000
  n_trial <- 5
  data <- simdata(dt, n_bin, n_trial, alpha, k)
  Y <- data$Y
  init_param <- list(log_k = -2,
                     log_a = 0,
                     x = prop_paths(Y, dt, -2, 0))
  fit <- fastr_fit(data=Y, dt=dt, simplified=F)
  expect_equal(fit$env$nlminb_fit$convergence, 0)
  expect_true(fit$env$tmb_report$pdHess)
})
