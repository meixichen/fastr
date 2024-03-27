devtools::load_all()
set.seed(123)
require(TMB)

dt <- 0.005
n_bin <- 1000
n_cell <- 8
n_trial <- 3
n_factor <- 2
n_basis <- 20
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

lam <- 1  # regularization
nu <- 5. # sigmoid function scale
T_len <- n_bin*dt  # total time length

# basis function
phi_x <- function(t, i){
  coef <- (i-0.5)*pi
  sqrt(2)*sin(coef*t/T_len)/coef
}

# Initial parameter values
Phi_args <- expand.grid(seq(dt, T_len, by=dt), 1:n_basis)
Phi <- t(matrix(apply(Phi_args, 1, function(row) phi_x(row[1], row[2])),
		nrow=n_bin, ncol=n_basis))  # n_basis x n_bin
Xi <- array(rnorm(n_basis*n_cell*n_trial), c(n_basis, n_cell, n_trial)) # n_basis x n_cell x n_trial 

data <- list(model="factor_model_basis", 
	     n_factor=n_factor,
	     dt=dt,
	     Y=Y,
	     lam=lam,
	     nu=nu,
	     Phi=Phi,
	     T=T_len)
init_param <- list(log_k=log(k),
		   log_a=log(alpha),
		   Lt=rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
		   Xi=Xi)

map <- list(log_k=rep(factor(NA), n_cell),
	    log_a=rep(factor(NA), n_cell))

adfun <- TMB::MakeADFun(data=data,
                        parameters=init_param,
                        map = map,
                        random = "Xi",
                        DLL = "fastr_TMBExports",
                        silent = FALSE)
# compare with two-stage fitting without basis expansion
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

# Without basis expansion it took around 40 secs
start_t <- Sys.time()
fit <- nlminb(adfun$par, adfun$fn, adfun$gr)
time_nlminb <- difftime(Sys.time(), start_t, units="secs")
start_t <- Sys.time()
rep <- sdreport(adfun, getJointPrecision=TRUE)
time_sdrep <- difftime(Sys.time(), start_t, units="secs") - time_nlminb
lmat_hat <- get_FA_estim(fit, n_cell=n_cell, n_factor=n_factor)$L
lmat_varimax <- varimax(lmat_hat)$loadings[1:n_cell,]
lmat_unnorm_hat <- rep$par.fixed[which(names(rep$par.fixed)=="Lt")]
lmat_unnorm_cov <- rep$cov.fixed[which(colnames(rep$cov.fixed)=="Lt"),
                               which(colnames(rep$cov.fixed)=="Lt")]
