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
    data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=1)
    
    #---- R nll --------------
    nll_r <- compute_Rnll(data, init_param)
    
    #---- TMB model nll -----------
    # Efficient matrix inversion not applied
    adfun <- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                            parameters=init_param,
                            DLL = "mnfa_TMBExports", 
                            silent = TRUE)
    nll_tmb <- adfun$fn(unlist(init_param)) # negative log-likelihood computed by TMB/C++
    expect_equal(nll_tmb, nll_r)
    # Efficient matrix inversion applied
    adfun_big <- TMB::MakeADFun(data=list(model="factor_model_big", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                            parameters=init_param,
                            DLL = "mnfa_TMBExports", 
                            silent = TRUE)
    
    nll_tmb_big <- adfun_big$fn(unlist(init_param))
    expect_equal(nll_tmb_big, nll_r)
  }
})

test_that("Optimization of the factor models converges and gives PD Hessian.",{
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
  
  # Efficient matrix inversion not applied
  adfun <- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                          parameters=init_param,
                          random = "x",
                          DLL = "mnfa_TMBExports", 
                          silent = T)
  mod_fit <- nlminb(adfun$par, adfun$fn, adfun$gr) # optimization
  rep <- TMB::sdreport(adfun) # get Hessian
  expect_equal(mod_fit$convergence, 0) # Expect convergence indicator is 0, i.e., nlminb has converged
  expect_true(rep$pdHess) # Expect the Hessian matrix is positive definite, i.e., no NA standard errors
  
  # Efficient matrix inversion applied
  adfun_big <- TMB::MakeADFun(data=list(model="factor_model_big", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                              parameters=init_param,
                              random = "x",
                              DLL = "mnfa_TMBExports", 
                              silent = T)
  mod_fit_big <- nlminb(adfun_big$par, adfun_big$fn, adfun_big$gr) # optimization
  rep_big <- TMB::sdreport(adfun_big) # get Hessian
  expect_equal(mod_fit_big$convergence, 0) # Expect convergence indicator is 0, i.e., nlminb has converged
  expect_true(rep_big$pdHess)
})
