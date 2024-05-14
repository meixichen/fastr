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
#' Sample from a multivariate normal with sparse precision matrix.
#'
#' @param n Number of random draws.
#' @param mean Mean vector.
#' @param prec Sparse precision matrix, i.e., inheriting from [Matrix::sparseMatrix] or its Cholesky factor, i.e., inheriting from [Matrix::CHMfactor].
#'
#' @return A matrix with `n` rows, each of which is a draw from the corresponding normal distribution.
#'
#' @details If the matrix is provided in precision form, it is converted to Cholesky form using `Matrix::Cholesky(prec, super = TRUE)`.  Once it is of form [Matrix::CHMfactor], this function is essentially copied from local function `rmvnorm()` in function `MC()` defined in [TMB::MakeADFun()].
rmvn_prec <- function(n, mean, prec) {
  d <- ncol(prec) # number of mvn dimensions
  if(!is(prec, "CHMfactor")) {
    prec <- Matrix::Cholesky(prec, super = TRUE)
  }
  u <- matrix(rnorm(d*n),d,n)
  u <- Matrix::solve(prec,u,system="Lt")
  u <- Matrix::solve(prec,u,system="Pt")
  u <- t(as(u, "matrix") + mean)
}
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
rm("pkg_template"); rm("ppd_template_test")
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
                          map=pkg_template_params$map,
                          random="x",
                          DLL=mod_name, silent = TRUE)

# Prepare mean and var of fixed effect posterior
fixed_mean <- unlist(theta_hat_list)
fixed_cov <- matrix(0, nrow=length(fixed_mean), ncol=length(fixed_mean))
fixed_cov[1:(2*n_cell), 1:(2*n_cell)] <- mod_fit$marg_cov
fixed_cov[(2*n_cell+1):nrow(fixed_cov), (2*n_cell+1):nrow(fixed_cov)] <- mod_fit$lmat_unnorm_cov
# 2-step PPD sampling
n_iter <- 1000
bi_spk_samples <- matrix(0, nrow=n_bin, ncol=n_trial*n_iter)
bi_spk_samples2 <- matrix(0, nrow=n_bin, ncol=n_trial*n_iter)
system.time({
  for (iter in 1:n_iter){
    cat("Sampling iteration:", iter, "\n")
    fixed_sample <- mvtnorm::rmvnorm(
      n = 1,
      mean = fixed_mean,
      sigma = fixed_cov)
    nll <- ppd_template$fn(drop(fixed_sample))
    
    # Method 1: sample X_i ~ N(hat(X_i), Sigma)
    X_mean <- rlang::with_env(ppd_template$env, last.par[lrandom()])
    X_prec <- rlang::with_env(ppd_template$env, L.created.by.newton)
    X_sam <- rmvn_prec(1, X_mean, X_prec)
    yi_sam <- matrix(0, nrow=n_bin, ncol=n_trial)
    for (tr in 1:n_trial){
      Xi_sam <- X_sam[,
                      seq(held_out_cell+(n_cell*n_bin*(tr-1)), n_cell*n_bin*tr, by=n_cell)]
      temp <- rep(0, n_bin)
      temp[get_spk_times_from_path(Xi_sam, exp(fixed_sample[held_out_cell]))] <- 1
      yi_sam[,tr] <- temp
    }
    bi_spk_samples[, (iter*n_trial-n_trial+1):(iter*n_trial)] <- yi_sam
    
    # # Method 2: let X_i = hat{X_i} (more efficient)
    # bi_spk_samples2[, (iter*n_trial-n_trial+1):(iter*n_trial)] <- ppd_template$simulate()$y_pred 
  }
})

# Joint PPD sampling
system.time({
  sdrep_pred <- sdreport(ppd_template, par.fixed=unlist(theta_hat_list$Lt), 
                         hessian.fixed = solve(mod_fit$env$tmb_report$cov.fixed),
                         getJointPrecision = TRUE, getReportCovariance = FALSE)
})
sdrep_pred$pdHess
n_iter <- 1e4
bi_spk_samples_j <- matrix(0, nrow=n_bin, ncol=n_trial*n_iter)
joint_prec <- sdrep_pred$jointPrecision
joint_mean <- c(sdrep_pred$par.fixed, sdrep_pred$par.random)
for (iter in 1:n_iter){
  cat("Sampling iteration:", iter, "\n")
  joint_sam <- rmvn_prec(1, joint_mean, joint_prec)
  X_sam <- joint_sam[, which(colnames(joint_sam)=="x")]
  yi_sam <- matrix(0, nrow=n_bin, ncol=n_trial)
  for (tr in 1:n_trial){
    Xi_sam <- X_sam[seq(held_out_cell+(n_cell*n_bin*(tr-1)), n_cell*n_bin*tr, by=n_cell)]
    temp <- rep(0, n_bin)
    temp[get_spk_times_from_path(Xi_sam, exp(theta_hat_list$log_k[held_out_cell]))] <- 1
    yi_sam[,tr] <- temp
  }
  bi_spk_samples_j[, (iter*n_trial-n_trial+1):(iter*n_trial)] <- yi_sam
}

# Cumulative spike count graphical GOF test
cum_spk_obs <- apply(Y[held_out_cell,,], 2, cumsum)
par(mfrow=c(1,3))
cum_spk_sample <- apply(bi_spk_samples_j, 2, cumsum)
for (u in 1:n_trial){
  plim <- range(cum_spk_sample[,u],cum_spk_obs[,u])
  plot(cum_spk_sample[,u], cum_spk_obs[,u], 
       xlim=plim, ylim=plim, type="l", 
       col=adjustcolor("gray", alpha.f=0.5),
       xlab="Replicated cumulative spike count",
       ylab="Observed cumulative spike count", main=paste("Trial", u))
  trial_iters <- seq(u, n_trial*n_iter, by=n_trial)
  avg_rep <- apply(cum_spk_sample[,trial_iters], 1, mean)
  for (iter in trial_iters){
    lines(cum_spk_sample[,iter], cum_spk_obs[,u], 
          col=adjustcolor("gray", alpha.f=0.5))
  }
  abline(0, 1, col="green", lty="dashed", lwd=3)
  lines(avg_rep, cum_spk_obs[,u], col="red", lwd=2)
}

## Change Y-axis to also replication
par(mfrow=c(1,3))
iter_sam <- sample(1:n_iter, 1)
for (u in 1:n_trial){
  trial_iters <- seq(u, n_trial*n_iter, by=n_trial)
  plim <- range(cum_spk_sample[,u])
  plot(cum_spk_sample[,u], cum_spk_sample[,trial_iters[iter_sam]], 
       xlim=plim, ylim=plim, type="l", 
       col=adjustcolor("gray", alpha.f=0.5),
       xlab="Replicated cumulative spike count",
       ylab="Replicated cumulative spike count", main=paste("Trial", u))
  avg_rep <- apply(cum_spk_sample[,trial_iters], 1, mean)
  for (iter in trial_iters){
    lines(cum_spk_sample[,iter], cum_spk_sample[,trial_iters[iter_sam]], 
          col=adjustcolor("gray", alpha.f=0.5))
  }
  abline(0, 1, col="green", lty="dashed", lwd=3)
  lines(avg_rep, cum_spk_sample[,trial_iters[iter_sam]], col="red", lwd=2)
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

# Check PSTH from replicated data
calc_psth <- function(binary_spks, bin_size=10){
  n_bin <- length(binary_spks)
  aggregate(binary_spks,
            list(cut(1:n_bin, (0:(n_bin/bin_size))*bin_size, right = TRUE)),
            sum)[,2]
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


# Spike time rep vs obs histograms
n_sam <- 5
obs_spk_times <- apply(Y[held_out_cell,,], 2, function(x) which(x==1)*dt)
for (tr in 1:n_trial){
  rep_spk_times <- apply(bi_spk_samples_j[,seq(tr, n_trial*n_iter,by=3)], 2, function(x) which(x==1)*dt)
  spk_inds <- sort(sample(1:min(c(length(obs_spk_times[[tr]]), sapply(rep_spk_times, length))), n_sam))
  spk_df <- data.frame(spike_time=obs_spk_times[[tr]][spk_inds], 
                       spike_index=as.factor(spk_inds),
                       type=rep("Observed", n_sam))
  rep_spk_times_mat <- sapply(spk_inds, function(spk_ind) sapply(rep_spk_times, function(rep) rep[spk_ind]))
  spk_df <- rbind(spk_df,
                  data.frame(spike_time=as.vector(rep_spk_times_mat),
                             spike_index=as.factor(rep(spk_inds, each=n_iter)),
                             type=rep("Replicated", n_sam*n_iter)))
  spk_df %>% filter(type=="Replicated") %>%
  ggplot(aes(x=spike_time, color=spike_index)) +
    geom_density() +
    geom_vline(data=spk_df%>%filter(type=="Observed"), 
               aes(xintercept=spike_time, color=spike_index), linetype="dotdash")+
    facet_grid(spike_index~.) +
    xlab("Spike time (in sec)")
}

