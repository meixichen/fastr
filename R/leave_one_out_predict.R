#' Leave-neuron-out posterior prediction.
#' @param mod_fit The fitted model of class `fastr_fit`.
#' @param left_out_neuron Nonzero integer index for the neuron to be predicted.
#' @return A function that can be run iteratively to sample from the posterior
#' predictive distribution.
#' @details
#' Additional details...
#' @export
predict.fastr_fit <- function(mod_fit, left_out_neuron, data, init_param,
                              method="2step", lam=NULL, nu=15, woodbury=TRUE,
                              silent=TRUE){
  if ("tmb_report" %in% names(mod_fit$env)){
    cov_fixed <- mod_fit$env$tmb_report$cov.fixed
  } else{
    cov_fixed <- mod_fit$lmat_unnorm_cov
  }
  # Check cov dimension to see if need to include log_a and log_k vars
  # ---- goes there
  
  # Build a new adfun obj for prediction without the left-out neuron data
  n_factor <- mod_fit$env$n_factor
  n_cell <- mod_fit$env$n_cell
  n_bin <- mod_fit$env$n_bin
  n_trial <- mod_fit$env$n_trial
  dt <- mod_fit$env$dt
  if (is.null(lam)) lam <- ifelse(n_cell<10, 1, 0.5)
  ppd_model <- fastr_model(data, dt=dt, 
                           n_factor=n_factor, init=init_param, 
                           method=method, lam=lam, nu=nu, woodbury=woodbury, 
                           left_out_neuron=left_out_neuron,
                           integrate_random=TRUE)
  ppd_template <- MakeADFun(data=ppd_model$data,
                            parameters=ppd_model$init_param,
                            map=ppd_model$map,
                            random=ppd_model$random,
                            DLL="fastr_TMBExports", silent = silent)
}