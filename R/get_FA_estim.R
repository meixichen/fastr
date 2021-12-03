#' Obtain L and a/k from model fit
#' 
#' @param mod_fit Fitted factor model output by nlminb(adfun)
#' @param n_cell Number of neurons
#' @param n_factor Number of factors
#' @return A list of normalized L matrix and the firing rate ratio (a/k)
#' @export

get_FA_estim <- function(mod_fit, n_cell, n_factor){
  estimates <- mod_fit$par
  off_diag_L_est <- estimates[names(estimates)=="Lt"]
  L_est <- matrix(0, nrow=n_cell, ncol=n_factor)
  L_est[lower.tri(L_est, diag=TRUE)] <- off_diag_L_est
  psi <- rep(1, n_cell)
  res <- norm_mat(L_est, psi, 3)
  
  fi <- unname(exp(estimates[names(estimates)=="log_a"]-
                     estimates[names(estimates)=="log_k"]))
  return(list(
    L=res$L,
    fi=round(fi,3)
  ))
}