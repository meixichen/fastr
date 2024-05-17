context("leave_one_out_predict")

test_that("Leave neuron out posterior prediction works for different settings. ",{
  set.seed(123)
  dt <- 0.005
  n_bin <- 200
  n_cell <- 8
  n_factor <- 2
  n_trial <- 1
  k <- runif(n_cell, 0.1, 0.6)
  alpha <- runif(n_cell, 2, 5.5)
  l1 <- runif(n_cell, 0.1, 0.25)
  l2 <- runif(n_cell, 0.1, 0.25)
  l1[1:(n_cell/n_factor)] <- runif(n_cell/n_factor, 0.7, 0.95)
  l2[(n_cell/n_factor+1):(n_cell/n_factor*2)] <- runif(n_cell/n_factor, 0.7, 0.95)
  L <- cbind(l1, l2)
  sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
  Y <- sim$Y
  x <- sim$x
  
  # First get parameter estimates
  mod_fit <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, report_sd = TRUE,
                       simplified = F)
  
  # LNO-PP test
  held_out_cell <- left_out_neuron <- as.integer(2)
  marg_param_conditions <- c(TRUE, FALSE)
  pred_methods <- c("mvn_joint", "mvn_2step_marg", "mvn_2step_both")
  all_test_conditions <- expand.grid(marg_param_conditions=marg_param_conditions, 
                                     pred_methods=pred_methods)
  n_iter <- 5
  time_taken <- rep(0, nrow(all_test_conditions))
  for (ii in 1:nrow(all_test_conditions)) {
    time_start <- Sys.time()
    lno_pp_fun <- leave_one_out_predict(mod_fit, left_out_neuron = 6, data = Y, 
                                        fix_marg_param = all_test_conditions[ii,1], 
                                        predict_method = as.character(all_test_conditions[ii,2]))
    for (iter in 1:n_iter){
      out <- lno_pp_fun()
    }
    time_taken[ii] <- as.numeric(difftime(Sys.time(), time_start, units="secs"))
  }
})
