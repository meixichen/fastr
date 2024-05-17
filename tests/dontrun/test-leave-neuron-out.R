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

#-------- PPD GOF test using TMB ----------------------
mod_fit <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, report_sd = TRUE,
                     simplified = F)
mod_fit <- readRDS("mod_fit.rds")
# Construct a new template with a held-out neuron
left_out_neuron <- as.integer(2)
#marg_param_conditions <- c(TRUE, FALSE)
#pred_methods <- c("mvn_joint", "mvn_2step_marg", "mvn_2step_both")
marg_param_conditions <- c(TRUE)
pred_methods <- c("mvn_joint")
all_test_conditions <- expand.grid(marg_param_conditions=marg_param_conditions, 
                                   pred_methods=pred_methods)
all_test_conditions$pred_methods <- as.character(all_test_conditions$pred_methods)
n_test <- nrow(all_test_conditions)
time_taken <- rep(0, n_test)
pred_spk_list <- vector("list", n_test)
n_iter <- 1000
for (ii in 1:n_test) {
  cat("=========== test", ii, "out of", n_test, "==========\n")
  time_start <- Sys.time()
  lno_pp_fun <- leave_one_out_predict(mod_fit, left_out_neuron = left_out_neuron, 
                                      data = Y, 
                                      fix_marg_param = all_test_conditions[ii,1], 
                                      predict_method = all_test_conditions[ii,2])
  bin_spks <- matrix(0, nrow=n_bin, ncol=n_trial*n_iter)
  for (iter in 1:n_iter){
    bin_spks[,(iter*n_trial-n_trial+1):(iter*n_trial)] <- lno_pp_fun()
  }
  pred_spk_list[[ii]] <- bin_spks
  time_taken[ii] <- as.numeric(difftime(Sys.time(), time_start, units="secs"))
}

# Cumulative spike count graphical GOF test
par(mfrow=c(1,n_trial))
for (ii in 1:n_test){
  for (u in 1:n_trial){
    pred <- pred_spk_list[[ii]][,seq(u, n_trial*n_iter, by=n_trial)]
    obs <- Y[left_out_neuron,,u]
    plot_gof_spk_csum(obs, pred)
  }
}

# Spike time rep vs obs histograms
pred <- pred_spk_list[[ii]][,seq(u, n_trial*n_iter, by=n_trial)]
obs <- Y[left_out_neuron, , u]
plot_gof_spk_times(obs, pred, dt, plot_type = "hist")

# Example of bad alignment
par(mfrow=c(1,2))
pred <- bin_spks_bad[,seq(1, n_trial*n_iter, by=n_trial)]
obs <- Y[left_out_neuron, , 2]
plot_gof_spk_csum(obs, pred)
pred <- bin_spks_bad[,seq(2, n_trial*n_iter, by=n_trial)]
obs <- Y[5, , 1]
plot_gof_spk_csum(obs, pred)
