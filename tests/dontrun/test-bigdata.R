# Test whether the model can perform efficiently in the large neuron ensemble scenario
set.seed(123)

#---------- Inputs from user -------------
model_flag <- readline(prompt = "Choose model flag from [hpp/hpp_woodbury/cpp_parallel]: ")

if (!(model_flag %in% c("hpp", "hpp_woodbury", "cpp_parallel"))){
  stop("model_flag must be one of 'hpp', 'hpp_woodbury', or 'cpp_parallel'.")
}

dt <- readline(prompt = "Enter dt: ")
n_bin <- readline(prompt = "Enter n_bin: ")
n_cell <- readline(prompt = "Enter n_cell: ")
if (n_cell%%n_factor != 0){
  stop("n_factor(4) must divide n_cell.")
}
n_trial <- readline(prompt = "Enter n_trial: ")
lam <- readline(prompt = "Enter penalty term lam: ")
neval <- readline(prompt = "Enter max number of function evaluations for nlminb: ")
niter <- readline(prompt = "Enter max number of iterations for nlminb: ")
n_factor <- 4 
#-------- Transform entered values to numeric/integer ----------------
dt <- as.numeric(dt)
n_bin <- as.integer(n_bin)
n_cell <- as.integer(n_cell)
n_trial <- as.integer(n_trial)
lam <- as.numeric(lam)
neval <- as.integer(neval)
niter <- as.integer(niter)

#-------- Data simulation -----------------------
cluster_size <- n_cell/n_factor
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
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
Y <- sim$Y
x <- sim$x
init_param <- list(log_k = rep(-1, n_cell),
                   log_a = rep(1, n_cell),
                   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                   x = prop_paths(Y, dt, rep(-1, n_cell), rep(1, n_cell)))
cat("Data simulated.\n")

#-------- Construct ADfun ---------------------
if (model_flag == "cpp_parallel"){
  require(TMB)
  mod_name <- "factor_model_parallel"
  compile(paste0(mod_name, ".cpp"))
  dyn.load(dynlib(mod_name))
  openmp(4) # use 4 threads
  data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=lam)
  DLL <- mod_name
}else if (model_flag == "hpp"){
  data <- list(model="factor_model", n_factor=n_factor, dt=dt, Y=Y, lam=lam)
  DLL <- "mnfa_TMBExports"
}else {
  data <- list(model="factor_model_big", n_factor=n_factor, dt=dt, Y=Y, lam=lam)
  DLL <- "mnfa_TMBExports"
}

adfun <- TMB::MakeADFun(data=data,
                        parameters=init_param,
                        random = "x",
                        DLL = DLL,
                        silent = FALSE)

cat("ADFun constructed. Ready for optimization. \n")

#--------- Optimization--------------------
t_start <- Sys.time()
mod_fit <- nlminb(adfun$par, adfun$fn, adfun$gr, control=list(eval.max=neval, iter.max=niter)) 
rep <- TMB::sdreport(adfun) 
t_taken <- Sys.time() - t_start
estim <- get_FA_estim(mod_fit, n_cell, n_factor)
obj <- list(true_alpha = alpha, true_k = k, true_L = L, 
            TMB_rep = rep, estimates = estim)
saveRDS(obj, file = "test-bigdata-results.rds") # save data

#--------- Display results ----------------
print(t_taken)
cat("The normalized loading matrix is \n")
print(estim$L)
cat("The true loading matrix is \n")
print(L)
cat("The estimated firing rates are", estim$fi, "\n")
cat("The true firing rates are", alpha/k, "\n")
