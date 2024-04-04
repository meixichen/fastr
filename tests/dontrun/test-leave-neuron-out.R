require(tidyverse)
require(dtw)
require(ggpubr)
require(MASS)
set.seed(123)
#-------- Simulate data -----------------
dt <- 0.005
n_bin <- 400
n_cell <- 16
n_factor <- 4
n_trial <- 3
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l3 <- runif(n_cell, 0.1, 0.25)
l4 <- runif(n_cell, 0.1, 0.25)
l1[1:(n_cell/n_factor)] <- runif(n_cell/n_factor, 0.7, 0.95)
l2[(n_cell/n_factor+1):(n_cell/n_factor*2)] <- runif(n_cell/n_factor, 0.7, 0.95)
l3[(n_cell/n_factor*2+1):(n_cell/n_factor*3)] <- runif(n_cell/n_factor, 0.7, 0.95)
l4[(n_cell/n_factor*3+1):(n_cell)] <- runif(n_cell/n_factor, 0.7, 0.95)
L <- cbind(l1, l2, l3, l4)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y
x <- sim$x

#-------- PPD for held-out neuron in the fitted time period -----------------
get_spk_times <- function(mod_fit, test_cell, test_trial,
                          X_sam_mean, X_sam_sd, test_stat_fun=median){
  n_bin <- mod_fit$env$n_bin
  k_sam <- exp(rnorm(n=1, 
                     mean=mod_fit$ig_params$log_k_hat[test_cell], 
                     sd=mod_fit$ig_params$logk_se[test_cell]))
  X_sam <- rnorm(n=n_bin, mean=X_sam_mean[test_cell,,test_trial],
                 sd=X_sam_sd[test_cell,,test_trial])
  temp <- floor(X_sam/k_sam)
  temp[which(temp<0)] <- 0
  spk_bin_ind <- sapply(1:max(temp),
                        function(u){
                          which(temp==u)[1]
                        })
  test_stat_fun(spk_bin_ind[!is.na(spk_bin_ind)])
}
fit_d <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, silent = TRUE)
X_sam_mean <- array(fit_d$paths, c(n_cell, n_bin, n_trial))
X_sam_sd <- array(fit_d$paths_se, c(n_cell, n_bin, n_trial))
test_cells <- sample(1:16, 3)
test_trials <- 1:3
n_sam <- 1e4

par(mfrow=c(3, 3))
for (test_cell in test_cells){
  for (test_trial in test_trials){
    test_T <- rep(0, n_sam)
    for (i in 1:n_sam){
      test_T[i] <- get_spk_times(fit_d, test_cell, test_trial, 
                                 X_sam_mean, X_sam_sd)  
    }
    obs_T <- median(which(Y[test_cell,,test_trial]==1))*dt
    hist(test_T*dt, 
         main=paste("Neuron", test_cell, "Trial", test_trial,
                    "p =", mean(test_T*dt < obs_T)))
    abline(v=obs_T, col="red")
  }
}

#-------- Predict held-out neuron in the fitted time period -----------------
get_spk_count <- function(binary_spk, bin_len){
  n_bin <- length(binary_spk)
  aggregate(binary_spk,
            list(cut(1:n_bin, (0:(n_bin/bin_len))*bin_len)),
            sum)[,2]
}

leave_neuron_out_pred <- function(mod_fit, test_cell, bin_len=20){
  n_factor <- mod_fit$env$n_factor
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
      dXq_recov[t, r] <- dXq_mean + cond_coeff %*% (dX_others[, t, r, drop=FALSE] - dX_others_mean)
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
  # Compare spike times using Dynamic Time Warping (DTW) distance
  obs_spk_times <- apply(Y[test_cell, , ], 2,
                         bin2num, dt=dt)
  recov_spk_times <- apply(recovered_spks, 2,
                           bin2num, dt=dt)
  
  mean_dtw_dist <- mean(sapply(1:n_trial, 
                     function(r){dtw(obs_spk_times[[r]], 
                                     recov_spk_times[[r]])$normalizedDistance}))
  
  recov_spk_counts <- apply(recovered_spks, 2, get_spk_count, bin_len=bin_len)
  obs_spk_counts <- apply(Y[test_cell, , ], 2, get_spk_count, bin_len=bin_len)
  mae_bin_count <- mean(sapply(1:n_trial,
                          function(r){
                            mean(abs(recov_spk_counts-
                                 obs_spk_counts))
                        }))
  list(mean_dtw_dist=mean_dtw_dist, mae_bin_count=mae_bin_count,
       recov_spk_times=recov_spk_times, recovered_spks=recovered_spks,
       recov_spk_counts)
}


factor_test_list <- 2:6
lno_results <- as.data.frame(expand.grid(n_factor=factor_test_list, test_cell=1:n_cell))
lno_results$dtw_dist <- NA
lno_results$mae_bin_count <- NA
for (d in factor_test_list){
  fit_d <- fastr_fit(data=Y, dt=dt, n_factor=d, silent = TRUE)
  for (q in 1:n_cell){
    result_dq <- leave_neuron_out_pred(fit_d, test_cell=q)
    row_selected <- which(lno_results$n_factor==d & lno_results$test_cell==q)
    lno_results$dtw_dist[row_selected] <- result_dq$mean_dtw_dist
    lno_results$mae_bin_count[row_selected] <- result_dq$mae_bin_count
  }
}

avg_lno_results <- lno_results %>% group_by(n_factor) %>% 
  summarise(mean_dtw_dist=mean(dtw_dist),
            mean_error_count=mean(mae_bin_count))
avg_dtw_plot <- ggplot(avg_lno_results) + 
  geom_line(aes(x=n_factor, y=mean_dtw_dist)) +
  labs(x="Number of factors", y="MAE for binned spike counts")
avg_error_plot <- ggplot(avg_lno_results) + 
  geom_line(aes(x=n_factor, y=mean_error_count))+
  labs(x="Number of factors", y="Mean of DTW distances between spike times")
ggarrange(avg_dtw_plot, avg_error_plot)
dtw_plot <- ggplot(lno_results) + 
  geom_line(aes(x=n_factor, y=dtw_dist, color=as.factor(test_cell))) +
  labs(x="Number of factors", y="MAE for binned spike counts")
error_plot <- ggplot(lno_results) + 
  geom_line(aes(x=n_factor, y=mae_bin_count, color=as.factor(test_cell)))+
  labs(x="Number of factors", y="Mean of DTW distances between spike times")
ggarrange(dtw_plot, error_plot)
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
    dXq_pred[t, r] <- dXq_mean + cond_coeff %*% (dX_test[, t, r, drop=FALSE] - dX_train_mean)
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
  obs_spk_counts <- get_spk_count(Y[test_cell,(n_train+1):n_bin, r, drop=FALSE], bin_len)
  plot(pred_spk_counts, ylim=range(pred_spk_counts, obs_spk_counts), 
       type="l", col="red", 
       xlab="Bin index", ylab="Spike count", main=paste("trial", r))
  lines(obs_spk_counts, col="blue")
  cat("Mean absolute error is", mean(abs(pred_spk_counts-obs_spk_counts)), "\n")
}

