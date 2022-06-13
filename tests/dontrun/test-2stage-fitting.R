set.seed(123)
require(TMB)
dt <- 0.005
n_bin <- 1000
n_cell <- 8
n_factor <- 2
n_trial <- 3
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l1[1:(n_cell/2)] <- runif(n_cell/2, 0.7, 0.95)
l2[(n_cell/2+1):n_cell] <- runif(n_cell/2, 0.7, 0.95)
L <- cbind(l1, l2)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y
x <- sim$x

all_mle <- get_ig_mle(Y, dt)
log_k_mle <- all_mle$log_k
log_a_mle <- all_mle$log_a
init_param <- list(log_k = log_k_mle,
                   log_a = log_a_mle,
                   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                   x = prop_paths(Y, dt, log_k_mle, log_a_mle))
cat("Data simulated.\n")

adfun<- TMB::MakeADFun(data=list(model="factor_model_eff", n_factor=n_factor, 
				 dt=dt, Y=Y, lam=1, nu=5.),
                                 parameters=init_param,
				 map = list(log_k = rep(factor(NA), n_cell),
					    log_a = rep(factor(NA), n_cell)),
                                 random = "x",
                                 DLL = "fastr_TMBExports",
                                 silent = F)

t_start <- Sys.time()
fit <- nlminb(adfun$par, adfun$fn, adfun$gr)
rep <- sdreport(adfun)
t_taken <- difftime(Sys.time(), t_start, units="mins")
cat("Time taken to fit is", t_taken)
varimax(get_FA_estim(fit, n_cell, n_factor)$L)
