devtools::load_all()
require(tidyverse)
require(dtw)
require(ggpubr)
require(MASS)
require(TMB)
#-------- Simulate data -----------------
set.seed(123)
dt <- 0.005
n_bin <- 1000
n_cell <- 12
n_factor <- 3
n_trial <- 3
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l3 <- runif(n_cell, 0.1, 0.25)
#l4 <- runif(n_cell, 0.1, 0.25)
l1[1:(n_cell/n_factor)] <- runif(n_cell/n_factor, 0.7, 0.95)
l2[(n_cell/n_factor+1):(n_cell/n_factor*2)] <- runif(n_cell/n_factor, 0.7, 0.95)
l3[(n_cell/n_factor*2+1):(n_cell/n_factor*3)] <- runif(n_cell/n_factor, 0.7, 0.95)
#l4[(n_cell/n_factor*3+1):(n_cell)] <- runif(n_cell/n_factor, 0.7, 0.95)
L <- cbind(l1, l2, l3)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y
x <- sim$x

#-------- PPD using TMB ----------------------
get_spk_times_from_path <- function(path, k){
  temp <- floor(path/k)
  temp[which(temp<0)] <- 0
  spk_bin_ind <- sapply(1:max(temp),
                        function(u){
                          which(temp==u)[1]
                        })
  spk_bin_ind
}
mod_name <- "factor_model_ppd"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name))

# Check if the cpp template gives the same results as the fastr hpp template
pkg_template_params <- fastr_model(data=Y, dt=dt, n_factor=n_factor)
pkg_template <- fastr_fit(data=Y, dt=dt, n_factor=n_factor,
                          init=pkg_template_params$init_param, 
                          adfun_only=TRUE, integrate_random = FALSE)

ppd_template_test <- MakeADFun(data=list(n_factor=n_factor,
                                    dt=dt,
                                    Y=Y, lam=ifelse(n_cell<10, 1, 0.5),
                                    nu=15, held_out_cell=0),
                          parameters=pkg_template_params$init_param,
                          map=pkg_template_params$map,
                          DLL=mod_name)
expect_equal(pkg_template$fn(pkg_template$par), ppd_template_test$fn(ppd_template_test$par))

# First get parameter estimates
# mod_fit <- fastr_fit(data=Y, dt=dt, n_factor=n_factor,
#                      init=pkg_template_params$init_param, report_sd = TRUE,
#                      simplified = F)
# saveRDS(mod_fit, "mod_fit.rds")
mod_fit <- readRDS("mod_fit.rds")
theta_hat_list <- list(log_k=mod_fit$log_k_hat,
                       log_a=mod_fit$log_a_hat,
                       Lt=mod_fit$lmat_unnorm_hat)
# Construct a new template with a held-out neuron
held_out_cell <- as.integer(6)
Y_train <- Y
Y_train[held_out_cell,,] <- matrix(NA, nrow=n_bin, ncol=n_trial)
init_params_train <- pkg_template_params$init_param
init_params_train$x[held_out_cell,,] <- init_params_train$x[1,,]
ppd_template <- MakeADFun(data=list(n_factor=n_factor,
                                    dt=dt,
                                    Y=Y_train, lam=ifelse(n_cell<10, 1, 0.5),
                                    nu=15, held_out_cell=held_out_cell),
                          parameters=init_params_train,
                          random="x",
                          DLL=mod_name, silent = TRUE)

# Prepare mean and var of fixed effect posterior
fixed_mean <- unlist(theta_hat_list)
fixed_cov <- matrix(0, nrow=length(fixed_mean), ncol=length(fixed_mean))
fixed_cov[1:(2*n_cell), 1:(2*n_cell)] <- mod_fit$marg_cov
fixed_cov[(2*n_cell+1):nrow(fixed_cov), (2*n_cell+1):nrow(fixed_cov)] <- mod_fit$lmat_unnorm_cov
# PPD sampling
calc_psth <- function(binary_spks, bin_size=10){
  n_bin <- length(binary_spks)
  aggregate(binary_spks, 
            list(cut(1:n_bin, (0:(n_bin/bin_size))*bin_size, right = TRUE)), 
            sum)[,2]
}
n_iter <- 1000
bi_spk_samples <- matrix(0, nrow=n_bin, ncol=n_trial*n_iter)
for (iter in 1:n_iter){
  cat("Sampling iteration:", iter, "\n")
  fixed_sample <- mvtnorm::rmvnorm(
    n = 1,
    mean = fixed_mean,
    sigma = fixed_cov)
  nll <- ppd_template$fn(drop(fixed_sample))
  random_sample <- rlang::with_env(ppd_template$env, last.par[lrandom()])
  bi_spk_samples[, (iter*n_trial-n_trial+1):(iter*n_trial)] <- ppd_template$simulate()$y_pred
}

bin_size <- 10
psth_rep_bounds <- matrix(0, nrow=n_bin/bin_size, ncol=n_trial*2)
colnames(psth_rep_bounds) <- apply(expand.grid(c("lower", "upper"), 1:n_trial), 1, paste, collapse="_")
for (u in 1:n_trial){
  psth_u <- sapply(1:n_iter, 
         function(iter) calc_psth(bi_spk_samples[,iter*n_trial-(n_trial-u)], bin_size=bin_size))
  psth_rep_bounds[, paste(c("lower", "upper"), u, sep = "_")] <- t(apply(psth_u, 1, range))
}
par(mfrow=c(3,1))
for (u in 1:n_trial){
  obs_psth <- calc_psth(Y[held_out_cell,,u], bin_size = bin_size)
  rep_psth <- psth_rep_bounds[, paste(c("lower", "upper"), u, sep = "_")]
  cover_percent <- mean(sapply(1:length(obs_psth), 
                               function(ii) {between(obs_psth[ii], rep_psth[ii,1], rep_psth[ii,2])}))
  plot(obs_psth, ylim=range(obs_psth, rep_psth), 
       main=paste("Trial", u, "Cover%:", cover_percent),
       ylab="Spike count", xlab="Bin index", col="red", type="l", lwd=2)
  polygon(c(1:length(obs_psth), length(obs_psth):1), c(rep_psth[,1], rev(rep_psth[,2])),
          col = adjustcolor("gray", alpha.f=0.5), border=FALSE)
}

# Cumulative spike count graphical GOF test
cum_spk_sample <- apply(bi_spk_samples, 2, cumsum)
cum_spk_obs <- apply(Y[held_out_cell,,], 2, cumsum)
par(mfrow=c(1,3))
for (u in 1:n_trial){
  plim <- range(cum_spk_sample[,u],cum_spk_obs[,u])
  plot(cum_spk_sample[,u], cum_spk_obs[,u], 
       xlim=plim, ylim=plim, type="l", 
       col=adjustcolor("gray", alpha.f=0.5),
       xlab="Replicated cumulative spike count",
       ylab="Observed cumulative spike count", main=paste("Trial", u))
  for (iter in seq(u, n_trial*100, by=3)){
    lines(cum_spk_sample[,iter], cum_spk_obs[,u], 
          col=adjustcolor("gray", alpha.f=0.5))
  }
  abline(0, 1, col="green", lty="dashed", lwd=3)
}

# Examples of bad alignment
par(mfrow=c(2,3))
match_trial_df <- cbind(c(1,2,3,2,3,1), c(2,3,1,1,2,3))
for (u in 1:nrow(match_trial_df)){
  plim <- range(cum_spk_sample[,match_trial_df[u,1]], cum_spk_obs[,match_trial_df[u,2]])
  plot(cum_spk_sample[,match_trial_df[u,1]], cum_spk_obs[,match_trial_df[u,2]], 
       xlim=plim, ylim=plim, type="l", 
       col=adjustcolor("gray", alpha.f=0.5),
       xlab="Replicated cumulative spike count",
       ylab="Observed cumulative spike count", 
       main=paste("Obs from trial", match_trial_df[u,2], "\n",
                  "Rep from trial", match_trial_df[u,1]))
  for (iter in seq(match_trial_df[u,1], n_trial*n_iter, by=3)){
    lines(cum_spk_sample[,iter], cum_spk_obs[,match_trial_df[u,2]], 
          col=adjustcolor("gray", alpha.f=0.5))
  }
  abline(0, 1, col="green", lty="dashed", lwd=3)
}


# Test sdreport
sdrep_pred <- sdreport(ppd_template, par.fixed=unlist(theta_hat_list))
x_ppd <- array(sdrep_pred$par.random, c(n_cell,n_bin,n_trial))
x_exceeds <- matrix(sdrep_pred$value[1:1200], nrow=n_bin, ncol=n_trial)
spk_probs <- matrix(sdrep_pred$value[1201:2400], nrow=n_bin, ncol=n_trial)
# Plot predicted path
par(mfrow=c(2, n_trial))
for (u in 1:n_trial){
  plot(x_ppd[held_out_cell,,u], type="l", col="red", 
       main="Red: predicted\n Green: true")
  lines(x[held_out_cell,,u], col="green")
}
for (u in 1:n_trial){
  obs_spks_ind <- which(Y[held_out_cell,,u]==1)
  plot(x_exceeds[,u], type="l")
  points(obs_spks_ind, rep(0, length(obs_spks_ind)), pch="|", col="red")
}

# Reset the predicted path using the estimated k_held_out
x_pred_reset <- matrix(0, nrow=n_bin, ncol=n_trial)
k_pred  <- exp(theta_hat_list$log_k[held_out_cell])
for (u in 1:n_trial){
  N_pred <- 1
  for (j in 1:n_bin){
    x_pred_reset[j, u] <- x_ppd[held_out_cell, j, u] - N_pred * k_pred
    N_pred <- N_pred + (x_pred_reset[j, u] >= 0)
  }
}
par(mfrow=c(1,n_trial))
for (u in 1:n_trial){
  obs_spks_ind <- which(Y[held_out_cell,,u]==1)
  plot(x_pred_reset[,u], type="l")
  points(obs_spks_ind, rep(0, length(obs_spks_ind)), pch="|", col="red")
}

# Predicted versus observed spike times (in terms of bin index)
# Output by resetting the estimated x path for the held-out neuron
get_spk_times_from_path(x_ppd[held_out_cell,,1], exp(theta_hat_list$log_k[held_out_cell]))
# Output from c++
which(diff(x_ppd[held_out_cell,,1]-x_exceeds[,1])!=0) 
# Truth
which(Y[held_out_cell,,1]==1)

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

