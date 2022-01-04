# Check the following:
# 1. does the parallel model produce correct nll
# 2. how much the parallel model is faster
require(TMB)
#require(mnfa)
mod_name <- "factor_model_parallel"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name)) 

dt <- 0.005
n_bin <- 2000
k <- c(0.15, 0.3, 0.2, 0.27, 0.25)
alpha <- c(2.5, 4, 3.5, 5, 3)
n_cell <- length(alpha)
n_trial <- 4
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
data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=1)

#---- Check nll --------------
nll_r <- compute_Rnll(data, init_param)
adfun_parallel <- TMB::MakeADFun(data=list(n_factor=n_factor, dt=dt, Y=Y, lam=1),
                        parameters=init_param,
                        DLL = mod_name, 
                        silent = FALSE)
nll_tmb_parallel <- adfun_parallel$fn(unlist(init_param))
print(nll_r - nll_tmb_parallel)

#---- Compare speed ----------------
adfun_parallel <- TMB::MakeADFun(data=list(n_factor=n_factor, dt=dt, Y=Y, lam=1),
                                 parameters=init_param,
                                 random = "x",
                                 DLL = mod_name,
                                 silent = T)
adfun_serial<- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                                 parameters=init_param,
                                 random = "x",
                                 DLL = "mnfa_TMBExports", 
                                 silent = T)
adfun_parallelhpp <- TMB::MakeADFun(data=list(model="factor_model_parallel", n_factor=n_factor, dt=dt, Y=Y, lam=1),
                              parameters=init_param,
                              random = "x",
                              DLL = "mnfa_TMBExports", 
                              silent = T)
t_serial <- system.time({
  fit_serial <- nlminb(adfun_serial$par, adfun_serial$fn, adfun_serial$gr)
})

t_parallel <- system.time({
  fit_parallel <- nlminb(adfun_parallel$par, adfun_parallel$fn, adfun_parallel$gr)
})

t_parallelhpp <- system.time({
  fit_parallelhpp <- nlminb(adfun_parallelhpp$par, adfun_parallelhpp$fn, adfun_parallelhpp$gr)
})

print(rbind(t_serial, t_parallel, t_parallelhpp)[,1:3])
