#' Leave-neuron-out posterior prediction.
#' @param mod_fit The fitted model of class `fastr_fit`.
#' @param left_out_neuron Nonzero integer index for the neuron to be predicted.
#' @return A function that can be run iteratively to sample from the posterior
#' predictive distribution.
#' @param data  A `n_cell x n_bin x n_trial` (multiple neurons) array of binary 
#' spike trains.
#' @param predict_method Character string "mvn_joint", "mvn_2step_marg", or 
#' "mvn_2step_both". See details.
#' @param fix_marg_param Are the marginal parameters (k and alpha) treated as 
#' fixed? Default is TRUE. When k and alpha are fixed, only the loading matrix
#' is sampled. See details.
#' @param woodbury Should the Woodbury matrix identity be used for covariance 
#' matrix operations? Default is TRUE. 
#' @param lam Scalar. Penalization parameter for the loading matrix elements. 
#' Default is `lam=ifelse(n_cell<10, 1, 0.5)`.
#' @param nu Scalar. Parameter that controls the "steepness" around 0 of the sigmoid
#' function applied at the spike data likelihood layer. Default is 15.
#' @param silent Suppress model fitting messages?
#' @param init_param Optional. A list of initial parameters for constructing the
#' model template (ADfun).
#' @param integrate_random Integrate random effects? If TRUE, latent paths are 
#' integrated out in the model. Setting to FALSE can be helpful for checking the 
#' marginal likelihood.
#' @param adfun_only Only outputs of ADFun created by TMB? This is for debugging
#' purposes only.
#' @details
#' Three method are available for leave-neuron-out posterior prediction. They
#' are listed below in terms of preference related to computational cost.
#' 1. `mvn_joint` samples parameters and latent paths jontly from a MVN:
#' \deqn{(\theta, X) \sim MVN\Big((\hat{\theta}, \hat{X}(\hat{\theta})),
#' \hat{\Sigma}_{\mathrm{joint}}\Big)}
#' where \eqn{\hat{X}(\hat{\theta})} is a first order Taylor expansion of 
#' \eqn{X(\theta)} around \eqn{\hat{\theta}}, and 
#' \eqn{\hat{\Sigma}_{\mathrm{joint}}} is a zero-th order Taylor expansion of 
#' the joint covariance at \eqn{\hat{\theta}}.
#' 
#' 2. `mvn_2step_marg` samples the paths in two steps. First, sample 
#' \eqn{\theta^s \sim q(\theta \mid y)}, where \eqn{q()} 
#' indicates an approximate posterior. 
#' Then conditional on \eqn{\theta^s} value, set \eqn{X^s = \hat{X}(\hat{\theta})}. 
#' I.e., only the parameters are sampled randomly whereas the paths are obtained
#' deterministically.
#' 
#' 3. `mvn_2step_both` samples the paths in two steps but with both steps 
#' performing random sampling. The first step is the same as that for 
#' `mvn_2step_marg`. In the second step, sample the latent paths
#' \eqn{X^s \sim MVN\Big(\hat{X}(\hat{\theta}), \hat{\Sigma}(\hat{\theta})\Big)}.
#' I.e., both parameters and paths are sampled randomly.
#' 
#' In terms of computational costs, methods 2 and 3 are similar in speed since 
#' the additional cost from randomly sampling step 2 is not massive. That said,
#' method 1 is preferred since it is much faster. 
#' @export
leave_one_out_predict <- function(mod_fit, left_out_neuron, data, 
                              predict_method=c("mvn_joint", "mvn_2step_marg", 
                                               "mvn_2step_both"), 
                              fix_marg_param=TRUE, woodbury=TRUE, lam=NULL, nu=15, 
                              silent=TRUE, init_param=NULL, integrate_random=TRUE,
                              adfun_only=FALSE){
  
  n_factor <- mod_fit$env$n_factor
  n_cell <- mod_fit$env$n_cell
  n_bin <- mod_fit$env$n_bin
  n_trial <- mod_fit$env$n_trial
  dt <- mod_fit$env$dt
  
  # Build cov.fixed and par.fixed
  cov1 <- mod_fit$marg_cov
  cov2 <- mod_fit$lmat_unnorm_cov
  if (fix_marg_param){ # the model template only takes in Lambda estimate
    par_fixed <- mod_fit$lmat_unnorm_hat
    cov_fixed <- cov2
  }
  else{ # the model template takes in both Lambda and marginal parameters
    par_fixed <- c(mod_fit$log_k_hat, mod_fit$log_a_hat, mod_fit$lmat_unnorm_hat)
    if (is.null(cov1)){ 
      # if marg_cov is not in the output, get SEs to cov 
      # but corr between log_a and log_k is lost
      cov1 <- diag(mod_fit$ig_params$loga_se, mod_fit$ig_params$logk_se)^2
    }
    #cov_fixed <- rbind(cbind(cov1,matrix(0,nrow=nrow(cov1),ncol=ncol(cov2))),
    #                   cbind(matrix(0,nrow=nrow(cov2),ncol=ncol(cov1)),cov2))
    cov_fixed <- matrix(0, nrow=length(par_fixed), ncol=length(par_fixed))
    cov_fixed[1:(2*n_cell), 1:(2*n_cell)] <- cov1
    cov_fixed[(2*n_cell+1):nrow(cov_fixed), 
              (2*n_cell+1):nrow(cov_fixed)] <- cov2
  }
  
  # Build a new adfun obj for prediction without the left-out neuron data
  if (is.null(lam)) lam <- ifelse(n_cell<10, 1, 0.5)
  ppd_model <- fastr_model(data, dt=dt, 
                           n_factor=n_factor, init=init_param, 
                           method=ifelse(fix_marg_param, "2step", "joint"), 
                           lam=lam, nu=nu, woodbury=woodbury, 
                           left_out_neuron=left_out_neuron,
                           integrate_random=integrate_random)
  ppd_template <- TMB::MakeADFun(data=ppd_model$data,
                                 parameters=ppd_model$init_param,
                                 map=ppd_model$map,
                                 random=ppd_model$random,
                                 DLL="fastr_TMBExports", silent = silent)
  if (adfun_only){
    return(ppd_template)
  } else{
    # Posterior predictive function to simulate from
    predict_method <- match.arg(predict_method)
    
    if (predict_method == "mvn_2step_both"){
      # Method 1: both marg param and X are randomly sampled
      # Step1: sample theta ~ N(theta_hat, Sigma_theta)
      # Step2: sample X_i ~ N(hat(X_i; theta), Sigma_X)
      out <- function(){
        fixed_sample <- mvtnorm::rmvnorm(
          n = 1,
          mean = par_fixed,
          sigma = cov_fixed)
        fixed_sample <- drop(fixed_sample)
        nll <- ppd_template$fn(fixed_sample)
        X_mean <- rlang::with_env(ppd_template$env, last.par[lrandom()])
        X_prec <- rlang::with_env(ppd_template$env, L.created.by.newton)
        X_sam <- rmvn_prec(1, X_mean, X_prec)
        yi_sam <- matrix(0, nrow=n_bin, ncol=n_trial)
        k_est <- ifelse(fix_marg_param, 
                        exp(mod_fit$log_k_hat[left_out_neuron]), 
                        exp(fixed_sample[left_out_neuron]))#FIX: NOT get k by position
        for (tr in 1:n_trial){
          Xi_sam <- X_sam[, seq(left_out_neuron+(n_cell*n_bin*(tr-1)), 
                                n_cell*n_bin*tr, by=n_cell)]
          temp <- rep(0, n_bin)
          temp[get_spk_times_from_path(Xi_sam, k_est)] <- 1
          yi_sam[,tr] <- temp
        }
        yi_sam
      }
    } else if (predict_method == "mvn_2step_marg"){
      # Method 2: only marg is randomly sampled and X is determined from marg
      # Step1: sample theta ~ N(theta_hat, Sigma_theta)
      # Step2: set X_i = argmax_x p(theta, X | Y) (laplace approximated)
      out <- function(){ 
        fixed_sample <- mvtnorm::rmvnorm(
          n = 1,
          mean = par_fixed,
          sigma = cov_fixed)
        nll <- ppd_template$fn(drop(fixed_sample))
        ppd_template$simulate()$y_pred 
      }
    } else{
      # Method 3: jointly sample marg and X
      # Step1: sample [theta, X_i] ~ N([theta, X_i(theta)], joint_Sig)
      sdrep_pred <- TMB::sdreport(ppd_template, 
                             par.fixed=par_fixed, 
                             hessian.fixed = solve(cov_fixed),
                             getJointPrecision = TRUE, 
                             getReportCovariance = FALSE)
      if (!sdrep_pred$pdHess) stop("Hess not PD! Prediction halted.")
      joint_prec <- sdrep_pred$jointPrecision
      joint_mean <- c(sdrep_pred$par.fixed, sdrep_pred$par.random)
      out <- function(){
        joint_sam <- rmvn_prec(1, joint_mean, joint_prec)
        X_sam <- joint_sam[, which(colnames(joint_sam)=="x")]
        yi_sam <- matrix(0, nrow=n_bin, ncol=n_trial)
        k_est <- ifelse(fix_marg_param, 
                        exp(mod_fit$log_k_hat[left_out_neuron]), 
                        exp(joint_sam[, left_out_neuron])) #FIX: NOT get k by position
        for (tr in 1:n_trial){
          Xi_sam <- X_sam[seq(left_out_neuron+(n_cell*n_bin*(tr-1)), 
                              n_cell*n_bin*tr, by=n_cell)]
          temp <- rep(0, n_bin)
          temp[get_spk_times_from_path(Xi_sam, k_est)] <- 1
          yi_sam[,tr] <- temp
        }
        yi_sam
      }
    }
    return(out)
  }
}