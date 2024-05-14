# Context of the test: 
# Make sure the joint MVN approx to the posterior predictive distribution (PPD)
# is doing what we hope it does. 

# Simple test model:
# theta ~ N(0, a^2)
# mu ~ N(0, b^2)
# y ~ N(mu, c^2)

# Sampling procedure:
# 1. Fixed effect sampling: theta ~ N(theta_hat, Sig_theta)
# 2. Conditional posterior: mu | theta, y ~ N(mu_n, sig_n)
# 3. Joint MVN approx sampling: (mu, theta) ~ p(mu, theta | theta_hat, y) 
# Note: step 3 is essentially using 1st Taylor expansion around theta_hat for 
# the mean and a zero-th order expansion around theta_hat for the cov

require(TMB)
mod_name <- "test_ppd_joint_mvn"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name))

log_a <- 2.
log_b <- 1.
log_c <- 0.5

n_obs <- 5000
theta <- 0.5
mu <- rnorm(1, theta, exp(log_b))
y <- rnorm(n_obs, mu, rep(exp(log_c), n_obs))
y_bar <- mean(y)
adf <- MakeADFun(data=list(y=y, log_a=log_a),
                 parameters=list(mu=mu, theta=theta, log_b=log_b, log_c=log_c),
                 random="mu",
                 map=list(log_b=factor(NA), log_c=factor(NA)),
                 DLL=mod_name)
fit <- nlminb(adf$par, adf$fn, adf$gr)

# Theoretical samples from p(mu, theta | theta_hat, y)
p_mu_given_y_theta <- function(n, theta, y){ # p(mu | y, theta)
  prec_n <- 1/exp(2*log_b) + length(y)/exp(2*log_c)
  var_n <- 1/prec_n
  mu_n <- (1/exp(2*log_b)) * var_n * theta + 
    (length(y)/exp(2*log_c)) * var_n * mean(y)
  rnorm(n, mu_n, sqrt(var_n))
}
n_rep <- 1000
mu_sam <- rep(0, n_rep)
pos_theta_mean <- 0.4
pos_theta_sd <- 2
for (ii in 1:n_rep){
  theta_s <- rnorm(1, pos_theta_mean, pos_theta_sd)
  mu_sam[ii] <- p_mu_given_y_theta(1, theta_s, y)
}
print(mean(mu_sam)); print(sd(mu_sam))

# TMB's way of doing the above using sdreport()
tmb_rep <- sdreport(adf, par.fixed=pos_theta_mean, hessian.fixed = 1/(pos_theta_sd^2), getJointPrecision = TRUE)
rmvn_prec <- function(n, mean, prec) {
  d <- ncol(prec) 
  if(!is(prec, "CHMfactor")) {
    prec <- Matrix::Cholesky(prec, super = TRUE)
  }
  u <- matrix(rnorm(d*n),d,n)
  u <- Matrix::solve(prec,u,system="Lt")
  u <- Matrix::solve(prec,u,system="Pt")
  u <- t(as(u, "matrix") + mean)
}
mu_sam2 <- rmvn_prec(n_rep, mean=c(tmb_rep$par.random, tmb_rep$par.fixed), prec=tmb_rep$jointPrecision)
print(mean(mu_sam2[,1])); print(sd(mu_sam2[,1]))
