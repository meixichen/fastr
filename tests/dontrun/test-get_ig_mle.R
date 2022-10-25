devtools::load_all()
require(numDeriv)
require(testthat)
set.seed(123)
####### Simulate data #############
dt <- 0.005
n_bin <- 2000
n_cell <- 20
n_factor <- 4
n_trial <- 1
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l3 <- runif(n_cell, 0.1, 0.25)
l4 <- runif(n_cell, 0.1, 0.25)
l1[1:5] <- runif(5, 0.7, 0.95)
l2[6:10] <- runif(5, 0.7, 0.95)
l3[11:15] <- runif(5, 0.7, 0.95)
l4[16:20] <- runif(5, 0.7, 0.95)
L <- cbind(l1, l2, l3, l4)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y

# Get se estimates analytically using get_ig_mle()
my_est <- get_ig_mle(Y, dt)
my_se <- my_est$se_ig

# Calculate the SEs numerically using numDeriv
negllk <- function(theta, ISI){
  gam <- theta[1]
  lam <- theta[2]
  n <- length(ISI)
  -.5*log(lam)*n + sum(.5* lam * gam^2 * ISI) - n*lam*gam + sum(.5*lam/ISI)
}

gam_se <- rep(NA, n_cell)
lam_se <- rep(NA, n_cell)
for (i in 1:n_cell){
  ISI <- bin2isi(Y,dt,i)
  x <- c(my_est$gamma[i], my_est$lambda[i])
  hess <- numDeriv::hessian(function(theta){negllk(theta, ISI)}, x)
  se <- sqrt(diag(solve(hess)))
  gam_se[i] <- se[1]
  lam_se[i] <- se[2]
}

testthat::expect_equal(gam_se, unname(my_se[1:n_cell]))
testthat::expect_equal(lam_se, unname(my_se[(1+n_cell):(2*n_cell)]))
