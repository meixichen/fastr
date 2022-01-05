# Test whether the model can perform efficiently in the large neuron ensemble scenario
set.seed(123)
require(TMB)
mod_name <- "factor_model_parallel"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name))
openmp(4) # use 4 threads 
dt <- 0.01
n_bin <- 1000
n_cell <- 16
n_factor <- 4
cluster_size <- n_cell/n_factor
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
n_trial <- 10
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l3 <- runif(n_cell, 0.1, 0.25)
l4 <- runif(n_cell, 0.1, 0.25)
cluster1 <- 1:cluster_size
cluster2 <- (cluster_size+1):(2*cluster_size)
cluster3 <- (2*cluster_size+1):(3*cluster_size)
cluster4 <- (3*cluster_size+1):(4*cluster_size)
l1[cluster1] <- runif(cluster_size, 0.7, 0.95)
l2[cluster2] <- runif(cluster_size, 0.7, 0.95)
l3[cluster3] <- runif(cluster_size, 0.7, 0.95)
l4[cluster4] <- runif(cluster_size, 0.7, 0.95)
L <- cbind(l1, l2, l3, l4)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
cat("Data simulated.\n")
Y <- sim$Y
x <- sim$x
init_param <- list(log_k = rep(-1, n_cell),
                   log_a = rep(1, n_cell),
                   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                   x = prop_paths(Y, dt, rep(-1, n_cell), rep(1, n_cell)))
adfun <- TMB::MakeADFun(data=list(n_factor=n_factor, dt=dt, Y=Y, lam=0.1),
                                 parameters=init_param,
                                 random = "x",
                                 DLL = mod_name,
                                 silent = FALSE)
t_start <- Sys.time()
mod_fit <- nlminb(adfun$par, adfun$fn, adfun$gr, control=list(eval.max=400, iter.max=400)) 
rep <- TMB::sdreport(adfun) 
t_taken <- Sys.time() - t_start
print("Time taken to fit the model is", t_taken)
estim <- get_FA_estim(mod_fit, n_cell, n_factor)
print("The normalized loading matrix is")
print(estim$L)
print("The true loading matrix is")
print(L)
cat("The estimated firing rates are", estim$fi, "\n")
cat("The true firing rates are", alpha/k, "\n")
