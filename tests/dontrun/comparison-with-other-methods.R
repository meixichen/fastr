require(TMB)
set.seed(123)
dt <- 0.005
n_bin <- 2000
n_cell <- 6
n_factor <- 2
n_trial <- 5
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l1[1:(n_cell/2)] <- runif(n_cell/2, 0.7, 0.95)
l2[(n_cell/2+1):n_cell] <- runif(n_cell/2, 0.7, 0.95)
L <- cbind(l1, l2)
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y
x <- sim$x
init_param <- list(log_k = rep(-2, n_cell),
                   log_a = rep(0, n_cell),
                   Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
                   x = prop_paths(Y, dt, rep(-2, n_cell), rep(0, n_cell)))
adfun <- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor,
                                        dt=dt, Y=Y, lam=1, nu=15.),
                                  parameters=init_param,
                                  random = "x",
                                  DLL = "mnfa_TMBExports",
                                  silent = F)
t_start <- Sys.time()
fit_tmb <- nlminb(adfun$par, adfun$fn, adfun$gr)
rep <- sdreport(adfun)
time_tmb <- difftime(Sys.time(), t_start, units="mins")



############# Competing method 1: simple k-means #########################
# Added Gaussian noise
Y.noisy <- Y + array(rnorm(length(Y), 0, 0.1), dim(Y))
km.res <- kmeans(Y.noisy[,,1], center=2)$cluster
# Convolved with exponentially decaying kernel
expdec.conv <- function(y, tau){
  T <- length(y)
  out <- rep(0, T)
  out[1] <- y[1]
  for (t in 2:T){
    out[t] <- ifelse(y[t]==1, 
                     out[t-1] + 1,
		     out[t-1]*exp(-1/tau)
                    )
  }
  out
}
Y.conv <- apply(Y, c(1,3), expdec.conv, tau=10)
km.res.conv <- kmeans(t(Y.conv[,,1]), center=2)$cluster
############# Competing method 2: fuzzy k-means ##########################
library(fclust)
fkm.res <- FKM(Y.noisy[,,1], k=2, m=1.2)$clus
fkm.res.conv <- FKM(t(Y.conv[,,1]), k=2, m=1.2)$clus

