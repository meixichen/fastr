devtools::load_all()
require(TMB)
set.seed(123)
dt <- 0.005
n_bin <- 2000
n_cell <- 6
n_factor <- 2
n_trial <- 5
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

# Fix k and a to MLEs
ig_mles <- get_ig_mle(Y, dt=dt)
k_hat <- sqrt(ig_mles$lam)
alpha_hat <- k_hat/ig_mles$mu
init_param2 <- list(log_k = log(k_hat),
                    log_a = log(alpha_hat),
                    Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                    x = prop_paths(Y, dt, log_k=log(k_hat), log_a=log(alpha_hat)))

############### Notations ################ 
# theta1 = (log(k), log(alpha)) 
## where k and alpha are length `n_cell` vectors
# theta2 = Lt 
## which is the computational basis of the upper triangle of the loading matrix
# L2 = joint likelihood (sigmoid) with y2 being the binary spike indicators
# L1 = marginal likelihood (IG) with y1 being the ISIs
#########################################

# Joint estimation (estimate theta1 and theta2 jointly)
adfun1 <- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor,
                                        dt=dt, Y=Y, lam=1, nu=15.),
                                  parameters=init_param2,
                                  random = "x",
                                  DLL = "mnfa_TMBExports",
                                  silent = F)
# Marginal estimation (with theta1 plugged in)
adfun2 <- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor,
                                        dt=dt, Y=Y, lam=1, nu=15.),
                                  parameters=init_param2,
                                  map = list(log_k = as.factor(rep(NA, n_cell)),
                                             log_a = as.factor(rep(NA, n_cell))),
                                  random = "x",
                                  DLL = "mnfa_TMBExports",
                                  silent = F)

t_start <- Sys.time()
fit_tmb2 <- nlminb(adfun2$par, adfun2$fn, adfun2$gr)
rep2 <- sdreport(adfun2)
time_tmb2 <- difftime(Sys.time(), t_start, units="mins")

L_est <- fit_tmb2$par
x_est <- rep2$par.random

marg_params <- list(log_k = log(k_hat),
		    log_a = log(alpha_hat),
		    Lt = L_est)
# \partial L2 / \partial theta1
gr_l2_t1 <- adfun1$gr(unlist(marg_params))[1:(2*n_cell)] 

R2inv <- rep2$cov.fixed #----------------------------------------------------- R2 inverse

# \partial L2 / \partial theta2
# checked it is only off from adfun1$gr(unlist(marg_params))[13:23] by e-8 to e-10 
gr_l2_t2 <- rep2$gradient.fixed 
R3 <- gr_l2_t1 %*% t(gr_l2_t2) #---------------------------------------------- R3

# A function to calculate the gradient: d L2 / d theta1
# @param X A list of length `n_cell` that contains the ISIs of all neurons
# @param log_k Vector of log threshold parameters for all neurons
# @param log_alpha Vector of log drift parameters for all neurons
# @return A vector of the gradient evaluated at the provided log_k and log_alpha
calc_gr_l1_t1 <- function(X, log_k, log_alpha){
  k <- exp(log_k)
  alpha <- exp(log_alpha)
  mu <- k/alpha
  lam <- k^2
  # \partial L1 / \partial (lam,mu)
  deriv_lam_mu <- function(x, lam, mu){
    n <- length(x)
    x_bar <- mean(x)
    # deriv w.r.t. mu and deriv w.r.t. lambda
    c(n*lam/mu^3*(x_bar-mu), 0.5*n/lam - 0.5/mu^2 * sum((x-mu)^2/x)) 
  }
  # \partial (mu, lam) / \partial (log(k), log(alpha))
  deriv_logk <- function(log_k, log_alpha){
    temp <- exp(log_k-log_alpha)
    matrix(c(temp, 2*exp(2*log_k), -temp, 0), ncol=2, nrow=2)
  }
  n_cell <- length(X)
  gr_logk <- rep(0, n_cell)
  gr_loga <- rep(0, n_cell)
  for (i in 1:n_cell){
    deriv1 <- deriv_lam_mu(X[[i]], lam[i], mu[i])
    deriv2 <- deriv_logk(log_k[i], log_alpha[i])
    gr_logk[i] <- t(deriv1) %*% deriv2[,1]
    gr_loga[i] <- t(deriv1) %*% deriv2[,2] 
  }
  c(gr_logk, gr_loga)
}
all_ISI <- apply(Y, 1,
		 function(yy) {
		   unlist(apply(yy, 2, function(y) { diff(which(y==1)*dt) }))
		 })
# We expect below to be zero
calc_gr_l1_t1(all_ISI, log(k_hat), log(alpha_hat))

R4 <- matrix(0, nrow=n_cell*2, ncol=length(L_est)) #---------------------------- R4

# Numerically calculate the Hessian of L1
library(numDeriv)
library(LaplacesDemon) # for the IG density
# Compute the inverse Gaussian log likelihood as a function of log(k) and log(a)
# @param theta Vector: c(log_k, log_alpha)
ig_loglik <- function(theta, ISI){
  #k <- exp(theta[1])
  #alpha <- exp(theta[2])
  #mu <- k/alpha 
  #lam <- k^2
  #n <- length(ISI)
  #n*theta[1] - 0.5 * lam/mu^2 * sum((ISI-mu)^2/ISI) + 0.5*n*log(2*pi)+1.5*sum(log(ISI))
  sum(dinvgaussian(all_ISI[[i]],
		   exp(theta[1]-theta[2]), 
		   exp(2*theta[1]), log=T))
}
R1 <- matrix(0, nrow=2*n_cell, ncol=2*n_cell)
for (i in 1:n_cell){
  hess_comp <- hessian(ig_loglik, c(log(k_hat[i]), log(alpha_hat[i])), 
		       ISI=all_ISI[[i]])
  idx <- (2*(i-1)+1):(2*i)
  R1[idx, idx] <- hess_comp
}
R1 <- -R1

# Another way to calculate the Hessian
fit <- optim(par=c(log(k_hat[i]), log(alpha_hat[i])),
	     # theta=c(log_k, log_alpha)
	     fn=function(theta){-sum(dinvgaussian(all_ISI[[i]],
						  exp(theta[1]-theta[2]), 
	   		                          exp(2*theta[1]), log=T))}, 
	     method="BFGS", hessian=TRUE)

# Actually analytical form of the Hessian in terms of mu and lam
# Hess <- diag(c(-n*lam/(mu^3), -0.5*n*lam^2)) # taken from the STAR package

# Finally, calculate the corrected variance estimator
R1inv <- solve(R1)
comp1 <- t(R3) %*% R1inv %*% R3
#comp2 <- t(R4) %*% R1inv %*% R3
#comp3 <- t(R3) %*% R1inv %*% R4
#Sig <- R2inv + R2inv %*% (comp1 - comp2 - comp3) %*% R2inv
Sig <- R2inv + R2inv %*% comp1 %*% R2inv # terms involving R4 are zero

