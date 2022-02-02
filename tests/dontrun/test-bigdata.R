# Test whether the model can perform efficiently in the large neuron ensemble scenario
require(TMB)
set.seed(123)
#---------- Inputs from user -------------
model_flag <- readline(prompt = "Choose model flag from [hpp/hpp_woodbury/cpp_parallel/cpp_bigparallel]: ")

if (!(model_flag %in% c("hpp", "hpp_woodbury", "cpp_parallel", "cpp_bigparallel"))){
  stop("model_flag must be one of 'hpp', 'hpp_woodbury', 'cpp_parallel', or 'cpp_bigparallel'.")
}

dt <- readline(prompt = "Enter dt: ")
dt <- as.numeric(dt)
n_bin <- readline(prompt = "Enter n_bin: ")
n_bin <- as.integer(n_bin)
n_cell <- readline(prompt = "Enter n_cell: ")
n_cell <- as.integer(n_cell)
n_factor <- readline(prompt = "Enter n_factor (currently only supporting n_factor=2, 4): ")
n_factor <- as.integer(n_factor)
if (n_cell%%n_factor != 0){
  stop("n_factor must divide n_cell.")
}
n_trial <- readline(prompt = "Enter n_trial: ")
n_trial <- as.integer(n_trial)
lam <- readline(prompt = "Enter penalty term lam: ")
lam <- as.numeric(lam)
neval <- readline(prompt = "Enter max number of function evaluations for nlminb: ")
neval <- as.integer(neval)
niter <- readline(prompt = "Enter max number of iterations for nlminb: ")
niter <- as.integer(niter)
filename <- readline(prompt = "Enter file name for saving R data: ")
if (tolower(tail(unlist(strsplit(filename, "[.]")), n=1)) != "rds"){
  filename <- paste0(filename, ".rds")
}

#-------- Data simulation -----------------------
cluster_size <- n_cell/n_factor
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
if (n_factor == 2){
  l1 <- runif(n_cell, 0.1, 0.25)
  l2 <- runif(n_cell, 0.1, 0.25)
  cluster1 <- 1:cluster_size
  cluster2 <- (cluster_size+1):(2*cluster_size)
  l1[cluster1] <- runif(cluster_size, 0.7, 0.95)
  l2[cluster2] <- runif(cluster_size, 0.7, 0.95)
  L <- cbind(l1, l2)
}else if (n_factor == 4){
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
}else{
 stop("Currently only supporting n_factor=2, 4.")
}
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
  mod_name <- "factor_model_parallel"
  compile(paste0(mod_name, ".cpp"))
  dyn.load(dynlib(mod_name))
  openmp(4) # use 4 threads
  data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=lam)
  DLL <- mod_name
}else if (model_flag == "cpp_bigparallel"){
  mod_name <- "factor_model_eff_parallel"
  compile(paste0(mod_name, ".cpp"))
  dyn.load(dynlib(mod_name))
  openmp(4) # use 4 threads
  data <- list(n_factor=n_factor, dt=dt, Y=Y, lam=lam)
  DLL <- mod_name
}else if (model_flag == "hpp"){
  data <- list(model="factor_model", n_factor=n_factor, dt=dt, Y=Y, lam=lam)
  DLL <- "mnfa_TMBExports"
}else {
  data <- list(model="factor_model_eff", n_factor=n_factor, dt=dt, Y=Y, lam=lam)
  DLL <- "mnfa_TMBExports"
}

adfun <- MakeADFun(data=data,
                        parameters=init_param,
                        random = "x",
                        DLL = DLL,
                        silent = FALSE)

cat("ADFun constructed. Ready for optimization. \n")

#--------- Optimization--------------------
t_start <- Sys.time()
mod_fit <- nlminb(adfun$par, adfun$fn, adfun$gr, control=list(eval.max=neval, iter.max=niter)) 
cat("nlminb optimization finised...running sdreport... \n")
rep <- sdreport(adfun) 
t_taken <- Sys.time() - t_start
estim <- get_FA_estim(mod_fit, n_cell, n_factor)
obj <- list(t_taken = t_taken, true_alpha = alpha, true_k = k, true_L = L, 
            fit = mod_fit, report = rep, estimates = estim)
saveRDS(obj, file = filename) # save data

#--------- Display results ----------------
print(t_taken)
cat("The normalized loading matrix is \n")
print(estim$L)
cat("The true loading matrix is \n")
print(L)
cat("The estimated firing rates are", estim$fi, "\n")
cat("The true firing rates are", alpha/k, "\n")
