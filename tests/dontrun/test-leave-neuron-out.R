set.seed(123)
require(TMB)
mod_name <- "1neuron_model"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name))

#-------- Simulate data -----------------
dt <- 0.005
n_bin <- 3000
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

#-------- Fitting training data -----------------
n_train <- 2000
Y_train <- Y[, 1:n_train, ]
fit_train <- fastr_fit(data=Y_train, dt=dt, n_factor=2)
Lam_train <- fit_train$lmat_hat
k_train <- exp(fit_train$ig_params$log_k_hat)
a_train <- exp(fit_train$ig_params$log_a_hat)

#-------- Fitting test data ------------------
test_cell <- 5
n_test <- n_bin-n_train
Y_test <- Y[, (n_train+1):dim(Y)[2],]
X_test <- array(NA, dim=c(n_cell-1, n_test, n_trial))
for (ii in 1:(n_cell-length(test_cell))){
  cell_ind <- setdiff(1:n_cell, test_cell)[ii]
  cat("Fitting neuron", cell_ind, "\n")
  fit_i <- fastr_fit(data=Y_test[cell_ind,,], dt=dt, silent = T)
  X_test_i <- array(fit_i$paths, dim=c(n_test, n_trial))
  X_test[ii,,] <- X_test_i
}


