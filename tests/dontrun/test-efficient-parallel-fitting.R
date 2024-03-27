# Check the following:
# 1. does the parallel model produce correct nll
# 2. how much the parallel model is faster

set.seed(123)
require(TMB)
#require(fastr)
mod_name1 <- "factor_model_eff_parallel"
compile(paste0(mod_name1, ".cpp"))
dyn.load(dynlib(mod_name1))
mod_name2 <- "factor_model_parallel"
compile(paste0(mod_name2, ".cpp"))
dyn.load(dynlib(mod_name2))

dt <- 0.005
n_bin <- 2000
n_cell <- 10
n_factor <- 2
n_trial <- 10
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l1[1:5] <- runif(5, 0.7, 0.95)
l2[6:10] <- runif(5, 0.7, 0.95)
L <- cbind(l1, l2)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y
x <- sim$x
init_param <- list(log_k = rep(-2, n_cell),
                   log_a = rep(0, n_cell),
                   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                   x = prop_paths(Y, dt, rep(-2, n_cell), rep(0, n_cell)))
data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=1, nu=5.)
cat("Data simulated.\n")
#---- Check nll --------------
nll_r <- compute_Rnll(data, init_param)
adfun_bigparallel <- TMB::MakeADFun(data=data,
                        parameters=init_param,
                        DLL = mod_name1,
                        silent = FALSE)
nll_tmb_bigparallel <- adfun_bigparallel$fn(unlist(init_param))
cat("Difference between r nll and tmb nll is", nll_r - nll_tmb_bigparallel, ". \n")

#---- Compare speed ----------------
adfun_bigparallel <- TMB::MakeADFun(data=list(n_factor=n_factor, dt=dt, Y=Y, lam=1, nu=5.),
                                 parameters=init_param,
                                 random = "x",
                                 DLL = mod_name1,
                                 silent = F)
adfun_eff <- TMB::MakeADFun(data=list(model="factor_model_eff", n_factor=n_factor,
				      dt=dt, Y=Y, lam=1, nu=5.),
                                 parameters=init_param,
                                 random = "x",
                                 DLL = "fastr_TMBExports",
                                 silent = F)
adfun_parallel <- TMB::MakeADFun(data=list(n_factor=n_factor, dt=dt, Y=Y, lam=1, nu=5.),
                                 parameters=init_param,
                                 random = "x",
                                 DLL = mod_name2,
                                 silent = F)
adfun_serial<- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor,
					dt=dt, Y=Y, lam=1, nu=5.),
                                 parameters=init_param,
                                 random = "x",
                                 DLL = "fastr_TMBExports",
                                 silent = F)
cat("All ADFuns constructed. Start fitting Woodbury + parallel. \n")
t_bigparallel <- system.time({
  fit_bigparallel <- nlminb(adfun_bigparallel$par, adfun_bigparallel$fn, adfun_bigparallel$gr)
})
cat("Finished fitting Woodbury + parallel. Time elapsed is", t_bigparallel[3], ". \n")
t_eff <- system.time({
  fit_eff <- nlminb(adfun_eff$par, adfun_eff$fn, adfun_eff$gr)
})
cat("Finished fitting efficient version. Time elapsed is", t_eff[3], ". \n")
t_parallel <- system.time({
  fit_parallel <- nlminb(adfun_parallel$par, adfun_parallel$fn, adfun_parallel$gr)
})
cat("Finished fitting parallel. Time elapsed is", t_parallel[3], ". \n")
t_serial <- system.time({
  fit_serial <- nlminb(adfun_serial$par, adfun_serial$fn, adfun_serial$gr)
})
cat("Finished fitting serial. Time elapsed is", t_serial[3], ". \n")


print(rbind(t_bigparallel, t_eff, t_parallel, t_serial)[,1:3])

