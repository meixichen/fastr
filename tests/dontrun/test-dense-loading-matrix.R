# Test the factor model with dense loading matrix 
require(TMB)
dt <- 0.005
n_bin <- 2000
k <- c(0.15, 0.3, 0.2, 0.27, 0.25, 0.1, 0.32, 0.4)
alpha <- c(2.5, 4, 3.5, 5, 3, 2.7, 6, 5.5)
n_cell <- length(alpha)
n_trial <- 5
l1 <- c(0.17, 0.08, 0.2, 0.1, 0.8, 0.2, 0.3, 0.7)
l2 <- c(0.72, 0.8, 0.05, 0.85, 0.1, 0.5, 0.1, 0.18)
l3 <- c(0.1, 0.1, 0.75,0.1, 0.15, 0.6, 0.7, 0.2)
L <- cbind(l1, l2, l3)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y
x <- sim$x
n_factor <- 3

#---------- Dense Lambda --------------
mod_name <- "factor_model_dense"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name))

init_param <- list(log_k = rep(-2, n_cell),
                   log_a = rep(0, n_cell),
                   Lt = rep(1, n_cell*n_factor),
                   x = prop_paths(Y, dt, rep(-2, n_cell), rep(0, n_cell)))

# Efficient matrix inversion not applied
adfun <- TMB::MakeADFun(data=list(n_factor=n_factor, dt=dt, Y=Y, lam=1, nu=5.),
                        parameters=init_param,
                        random="x",
                        DLL = mod_name, 
                        silent = T)

t_start <- Sys.time()
mod_fit <- nlminb(adfun$par, adfun$fn, adfun$gr) # optimization
rep <- TMB::sdreport(adfun) # get Hessian
t_taken <- Sys.time() - t_start
print("Time taken to fit the dense matrix model is", t_taken)
print(summary(rep, "fixed"))
#---------- Lower tri Lambda -----------
mod_name <- "factor_model_parallel"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name))

init_param <- list(log_k = rep(-2, n_cell),
                   log_a = rep(0, n_cell),
                   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                   x = prop_paths(Y, dt, rep(-2, n_cell), rep(0, n_cell)))

# Efficient matrix inversion not applied
adfun_lt <- TMB::MakeADFun(data=list(n_factor=n_factor, dt=dt, Y=Y, lam=1, nu=5.),
                        parameters=init_param,
                        random="x",
                        DLL = mod_name, 
                        silent = T)

t_start <- Sys.time()
mod_fit_lt <- nlminb(adfun_lt$par, adfun_lt$fn, adfun_lt$gr) # optimization
rep_lt <- TMB::sdreport(adfun_lt) # get Hessian
t_taken <- Sys.time() - t_start
print("Time taken to fit the lower triangular matrix model is", t_taken)
print(summary(rep_lt, "fixed"))
