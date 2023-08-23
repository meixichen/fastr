devtools::load_all()
set.seed(123)
require(TMB)

n_bin <- 1000
n_cell <- 8
n_trial <- 3
n_factor <- 2
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

dt <- 0.005
Y <- NULL  # q x n x r
lam <- NULL  # regularization
nu <- 20
Phi <- NULL  # m x n
T_len <- 10

log_k <- -1
log_a <- 0
Lt <- NULL  # dq-d(d-1)/2 vector
Xi <- NULL  # m x q x r 

data <- list(model="factor_model_basis", 
	     n_factor=n_factor,
	     dt=dt,
	     Y=Y,
	     lam=lam,
	     nu=nu,
	     Phi=Phi,
	     T=T_len)
init_param <- list(log_k=log_k,
		   log_a=log_a,
		   Lt=Lt,
		   Xi=Xi)

map <- list(log_k=rep(factor(NA), n_cell),
	    log_a=rep(factor(NA), n_cell))

adfun <- TMB::MakeADFun(data=data,
                        parameters=init_param,
                        map = map,
                        random = "Xi",
                        DLL = "fastr_TMBExports",
                        silent = FALSE)

start_t <- Sys.time()
fit <- nlminb(adfun$par, adfun$fn, adfun$gr, control = list(eval.max=1000, iter.max=1000))
time_nlminb <- difftime(Sys.time(), start_t, units="secs")
rep <- sdreport(adfun, getJointPrecision=TRUE)
time_sdrep <- difftime(Sys.time(), start_t, units="secs") - time_nlminb
lmat_hat <- get_FA_estim(fit, n_cell=n_cell, n_factor=n_factor)$L
lmat_varimax <- varimax(lmat_hat)$loadings[1:n_cell,]
lmat_unnorm_hat <- rep$par.fixed[which(names(rep$par.fixed)=="Lt")]
lmat_unnorm_cov <- rep$cov.fixed[which(colnames(rep$cov.fixed)=="Lt"),
                               which(colnames(rep$cov.fixed)=="Lt")]
