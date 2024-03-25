#' Obtain L and a/k from model fit
#' 
#' @param mod_fit Fitted factor model output by nlminb(adfun), or a vector of off
#' diagonal elements of the loading matrix.
#' @param n_cell Number of neurons
#' @param n_factor Number of factors
#' @return A list of normalized L matrix and, if the model fit is provided, the firing 
#' rate ratio (a/k)
#' @export

get_FA_estim <- function(mod_fit, n_cell, n_factor){
  if (!("par" %in% names(mod_fit))){
    if (length(mod_fit) == n_cell*n_factor-n_factor*(n_factor-1)/2){ 
      off_diag_L_est <- mod_fit
      fi <- NULL
    } else{
      stop("mod_fit should be either an output of nlminb/optim, or a vector of Lt estimates.")
    }
  } else{
    estimates <- mod_fit$par
    off_diag_L_est <- estimates[names(estimates)=="Lt"]
    fi <- unname(exp(estimates[names(estimates)=="log_a"]-
                     estimates[names(estimates)=="log_k"]))
    fi <- round(fi, 3)
  }
  L_est <- matrix(0, nrow=n_cell, ncol=n_factor)
  L_est[lower.tri(L_est, diag=TRUE)] <- off_diag_L_est
  psi <- rep(1, n_cell)
  res <- norm_mat(L_est, psi, 3)
  
  return(list(
    L=res$L,
    fi=fi
  ))
}
