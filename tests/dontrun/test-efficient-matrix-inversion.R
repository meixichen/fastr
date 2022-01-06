set.seed(123)
dt <- 0.01
n_bin <- 50
n_cell <- 300
n_trial <- 5
n_factor <- 3
if (is.integer(n_cell/n_factor)) {stop("n_factor must divide n_cell.")}
cluster_size <- n_cell/n_factor
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 3, 5.5)
l1 <- runif(n_cell, 0.05, 0.2)
l2 <- runif(n_cell, 0.05, 0.2)
l3 <- runif(n_cell, 0.05, 0.2)
cluster1 <- 1:cluster_size
cluster2 <- (cluster_size+1):(2*cluster_size)
cluster3 <- (2*cluster_size+1):(3*cluster_size)
l1[cluster1] <- runif(cluster_size, 0.8, 0.95)
l2[cluster2] <- runif(cluster_size, 0.8, 0.95)
l3[cluster3] <- runif(cluster_size, 0.9, 0.95)
L <- cbind(l1, l2, l3)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y
x <- sim$x
init_param <- list(log_k = rep(-1, n_cell),
                   log_a = rep(1, n_cell),
                   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                   x = prop_paths(Y, dt, rep(-1, n_cell), rep(1, n_cell)))
cat("Data simulated.\n")
adfun <- TMB::MakeADFun(data=list(model="factor_model", n_bin=n_bin, n_cell=n_cell, n_trial=n_trial, n_factor=n_factor, dt=dt, Y=Y, lam=0.1),
                        parameters=init_param,
                        DLL = "mnfa_TMBExports", silent = TRUE)
adfun_eff <- TMB::MakeADFun(data=list(model="factor_model_big", n_bin=n_bin, n_cell=n_cell, n_trial=n_trial, n_factor=n_factor, dt=dt, Y=Y, lam=0.1),
                            parameters=init_param,
                            DLL = "mnfa_TMBExports", silent = TRUE)
cat("Adfun constructed.\n")
n_rep <- 10
n_test <- 100
t_reg_mat <- NULL
t_eff_mat <- NULL
for (i in 1:n_rep){
  t_reg <- system.time({
    for (ii in 1:n_test){
      adfun$fn(unlist(init_param))
    }
  })
  t_reg_mat <- rbind(t_reg_mat, t_reg[1:3])
  t_eff <- system.time({
    for (ii in 1:n_test){
      adfun_eff$fn(unlist(init_param))
    }
  })
  t_eff_mat <- rbind(t_eff_mat, t_eff[1:3])
}

colnames(t_reg_mat) <- paste("reg", colnames(t_reg_mat))
colnames(t_eff_mat) <- paste("eff", colnames(t_eff_mat))
print(cbind(t_reg_mat, t_eff_mat))
