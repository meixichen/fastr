set.seed(123)
require(dplyr)
require(dtw)
#-------- Simulate data -----------------
dt <- 0.005
n_bin <- 400
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

#-------- Predict held-out neuron in the fitted time period -----------------
test_cell <- 3
# Fit the model to all data
fit_all <- fastr_fit(data=Y, dt=dt, n_factor=2)
Lam_hat <- fit_all$lmat_hat
k_hat <- exp(fit_all$ig_params$log_k_hat)
a_hat <- exp(fit_all$ig_params$log_a_hat)
cov_hat <- Lam_hat %*% t(Lam_hat)
diag(cov_hat) <- 1

# Extract estimated paths
X_hat <- array(fit_all$paths, dim=c(n_cell, n_bin, n_trial))
dX_others <- apply(X_hat[-test_cell,,], MARGIN=c(1,3), function(x) c(x[1], diff(x)))
dX_others <- aperm(dX_others, c(2, 1, 3))
dXq_recov <- matrix(0, nrow=n_bin, ncol=n_trial)
cond_coeff <- cov_hat[test_cell, -test_cell, drop=F]%*%
  solve(cov_hat[-test_cell, -test_cell]) 
dXq_mean <- a_hat[test_cell]*dt
dX_others_mean <- a_hat[-test_cell]*dt
for (t in 1:n_bin){
  for (r in 1:n_trial){
    dXq_recov[t, r] <- dXq_mean + cond_coeff %*% (dX_others[, t, r] - dX_others_mean)
  }
}
Xq_recov <- apply(dXq_recov, 2, cumsum)
# Recovered spikes per trial for the held-out neuron
recovered_spks <- apply(Xq_recov, 2, 
                        function(x){
                          temp <- floor(x/k_hat[test_cell])
                          temp[which(temp<0)] <- 0
                          c(temp[1], diff(temp))
                        })
apply(recovered_spks, 2, sum) # Recovered total spike counts
apply(Y[test_cell,,], 2, sum) # Observed total spike counts
# Compare spike times
obs_spk_times <- apply(Y[test_cell, , ], 2,
                       bin2num, dt=dt)
recov_spk_times <- apply(recovered_spks, 2,
                         bin2num, dt=dt)
par(mfrow=c(1,3))
for (r in 1:n_trial){
  obs <- obs_spk_times[[r]]
  pred <- recov_spk_times[[r]]
  plot(obs, col="blue", 
       xlim=range(1, length(obs), length(pred)),
       ylim=range(obs, pred), xlab="Spike index", ylab="Spike time") 
  points(pred, col="red")
  cat("Dynamic time warping distance between the two is", 
      dtw(obs, pred)$normalizedDistance, "\n")
}

# Spike counts in bins
get_spk_count <- function(binary_spk, bin_len){
  n_bin <- length(binary_spk)
  aggregate(binary_spk,
            list(cut(1:n_bin, (0:(n_bin/bin_len))*bin_len)),
            sum)[,2]
}
bin_len <- 20
par(mfrow=c(1,3))
for (r in 1:n_trial){
  recov_spk_counts <- get_spk_count(recovered_spks[, r], bin_len)
  obs_spk_counts <- get_spk_count(Y[test_cell, , r], bin_len)
  plot(recov_spk_counts, ylim=range(recov_spk_counts, obs_spk_counts), 
       type="l", col="red", 
       xlab="Bin index", ylab="Spike count", main=paste("trial", r))
  lines(obs_spk_counts, col="blue")
  cat("Mean absolute error is", mean(abs(recov_spk_counts-obs_spk_counts)), "\n")
}

#-------- Predict held-out neuron in held-out time period ------------------
test_cell <- 2
# Fit the training data
n_train <- 400
Y_train <- Y[, 1:n_train, ]
fit_train <- fastr_fit(data=Y_train, dt=dt, n_factor=2)
Lam_train <- fit_train$lmat_hat
k_train <- exp(fit_train$ig_params$log_k_hat)
a_train <- exp(fit_train$ig_params$log_a_hat)
cov_train <- Lam_train %*% t(Lam_train)
diag(cov_train) <- 1

# Fit a single-neuron model to each training neuron in the held-out time period
train_cells <- setdiff(1:n_cell, test_cell)
n_test <- n_bin-n_train
Y_test <- Y[, (n_train+1):dim(Y)[2],]
X_test <- array(NA, dim=c(n_cell-1, n_test, n_trial))
for (ii in 1:(n_cell-length(test_cell))){
  cell_ind <- train_cells[ii]
  cat("Fitting neuron", cell_ind, "\n")
  fit_i <- fastr_fit(data=Y_test[cell_ind,,], dt=dt, silent = T)
  X_test_i <- array(fit_i$paths, dim=c(n_test, n_trial))
  X_test[ii,,] <- X_test_i
}

# Predict activity of the test neuron in the held-out time period
dX_test <- apply(X_test, MARGIN=c(1,3), function(x) c(x[1], diff(x)))
dX_test <- aperm(dX_test, c(2, 1, 3))
dXq_pred <- matrix(0, nrow=n_test, ncol=n_trial)
cond_coeff <- cov_train[test_cell, -test_cell, drop=F]%*%
  solve(cov_train[-test_cell, -test_cell]) 
dXq_mean <- a_train[test_cell]*dt
dX_train_mean <- a_train[-test_cell]*dt
for (t in 1:n_test){
  for (r in 1:n_trial){
    dXq_pred[t, r] <- dXq_mean + cond_coeff %*% (dX_test[, t, r] - dX_train_mean)
  }
}
Xq_pred <- apply(dXq_pred, 2, cumsum)
# Predicted spikes per trial for the held-out neuron
pred_spks_test <- apply(Xq_pred, 2, 
                        function(x){
                          temp <- floor(x/k_train[test_cell])
                          temp[which(temp<0)] <- 0
                          c(temp[1], diff(temp))
                        })
apply(pred_spks_test, 2, sum)
# Observed total spike count per trial for the held-out neuron
apply(Y[test_cell,(n_train+1):n_bin,], 2, sum)

# Not a good idea to compare ISIs because they are determined only by k and alpha,
# so the estimated correlation does not contribute much to prediction

# Compare spike times
obs_spk_times <- apply(Y[test_cell, (n_train+1):n_bin, ], 2,
                       bin2num, dt=dt)
pred_spk_times <- apply(pred_spks_test, 2,
                        bin2num, dt=dt)
par(mfrow=c(1,3))
for (r in 1:n_trial){
  obs <- obs_spk_times[[r]]
  pred <- pred_spk_times[[r]]
  plot(obs, col="blue", 
       xlim=range(1, length(obs), length(pred)),
       ylim=range(obs, pred), xlab="Spike index", ylab="Spike time") 
  points(pred, col="red")
}

# Spike counts in bins
bin_len <- 50
par(mfrow=c(1,3))
for (r in 1:n_trial){
  pred_spk_counts <- get_spk_count(pred_spks_test[, r], bin_len)
  obs_spk_counts <- get_spk_count(Y[test_cell,(n_train+1):n_bin, r], bin_len)
  plot(pred_spk_counts, ylim=range(pred_spk_counts, obs_spk_counts), 
       type="l", col="red", 
       xlab="Bin index", ylab="Spike count", main=paste("trial", r))
  lines(obs_spk_counts, col="blue")
  cat("Mean absolute error is", mean(abs(pred_spk_counts-obs_spk_counts)), "\n")
}

