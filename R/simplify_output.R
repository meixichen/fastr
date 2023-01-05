#' Makes simplified output from an object of class `fastr_fit`
#'
#' @param fobj An object of class `fastr_fit`.
#' @param ... Additional arguments.
#' @return A list of class `fastr_fit`.
#' @details This is both a helper function used within the function `fastr_fit()` and
#' can also be used standalone.
#' @export 

simplify_output <- function(fobj, ...){
  if (!inherits(fobj, "fastr_fit")) stop("`fobj` must be an object of class `fastr_fit`.")
  n_cell <- fobj$env$n_cell
  ig_params <- list(lam_hat = fobj$lam_hat, 
		    lam_se = fobj$se_ig[(n_cell+1):(2*n_cell)],
                    log_k_hat = fobj$log_k_hat,
		    log_a_hat = fobj$log_a_hat,
		    logk_se = fobj$logk_se,
		    loga_se = fobj$loga_se
                    )
  env <- list(start_time = fobj$env$start_time,
              time_nlminb = fobj$env$time_nlminb,
	      time_sdrep = fobj$env$time_sdrep,
	      n_factor = fobj$env$n_factor,
	      n_cell = fobj$env$n_cell,
	      n_bin = fobj$env$n_bin,
	      n_trial = fobj$env$n_trial,
	      dt = fobj$env$dt
              )
  out <- list(time = fobj$time, 
	      rate_hat = fobj$rate_hat,
              rate_se = fobj$se_ig[1:n_cell],
              lmat_hat = fobj$lmat_hat,
	      lmat_varimax = fobj$lmat_varimax,
	      lmat_unnorm_hat = fobj$lmat_unnorm_hat,
	      lmat_unnorm_cov = fobj$lmat_unnorm_cov,
	      paths = fobj$paths,
	      paths_se = fobj$paths_se,
	      rho_se = fobj$env$tmb_report$sd,
	      ig_params = ig_params,
	      env = env
              )
  class(out) <- "fastr_fit"
  return(out)
}
