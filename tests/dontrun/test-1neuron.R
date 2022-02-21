set.seed(123)
require(TMB)
#require(mnfa)
mod_name <- "1neuron_model"
compile(paste0(mod_name, ".cpp"))
dyn.load(dynlib(mod_name))

k <- 0.4
alpha <- 1.2
dt <- 0.05
n_bin <- 2000
n_trial <- 5

data <- simdata(dt, n_bin, n_trial, alpha, k)
Y <- data$Y
init_param <- list(log_k = -2,
                   log_a = 0,
                   x = prop_paths(Y, dt, -2, 0))
data <- list(dt=dt, Y=Y, nu=5)

adfun <- TMB::MakeADFun(data=data,
                        parameters=init_param,
                        random = "x",
                        DLL = mod_name,
                        silent = F)
fit <- nlminb(adfun$par, adfun$fn, adfun$gr)
rep <- sdreport(adfun)
exp(summary(rep, "fixed"))

