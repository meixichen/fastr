context("factor_model_big")

test_that("Big factor model in C++ has same nll as the nll computed in R. ",{
  n_test <- 10
  for (ii in 1:n_test){
    #---- Data simulation ---------
    dt <- 1
    n_bin <- 2
    n_cell <- sample(2:10, 1)
    k <- runif(n_cell, 0.1, 0.4)
    alpha <- runif(n_cell, 2, 5)
    n_trial <- 1
    Sig <- cov2cor(rWishart(1, n_cell+1, diag(n_cell))[,,1])
    rho <- Sig[lower.tri(Sig)]
    sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, rho=rho)
    Y <- sim$Y
    x <- sim$x
    
    #---- TMB model nll -----------
    n_factor <- ifelse(n_cell >= 5, 2, 1)
    init_param <- list(log_k = rep(-2, n_cell),
                       log_a = rep(0, n_cell),
                       Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                       x = prop_paths(Y, dt, rep(-2, n_cell), rep(0, n_cell)))
    data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=1)
    adfun_reg <- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                                parameters=init_param,
                                DLL = "mnfa_TMBExports", 
                                silent = TRUE)
    adfun_big <- TMB::MakeADFun(data=list(model="factor_model_big", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                            parameters=init_param,
                            DLL = "mnfa_TMBExports", 
                            silent = TRUE)
    nll_reg <- adfun_reg$fn(unlist(init_param)) # negative log-likelihood computed by TMB/C++
    nll_big <- adfun_big$fn(unlist(init_param))
    
    #---- R nll --------------
    nll_r <- compute_Rnll(data, init_param)
    nll_r-nll_reg
    nll_r-nll_big
  }
})
