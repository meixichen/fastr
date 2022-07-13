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
                                  DLL = "fastr_TMBExports",
                                  silent = F)
t_start <- Sys.time()
fit_tmb <- nlminb(adfun$par, adfun$fn, adfun$gr)
rep <- sdreport(adfun)
time_tmb <- difftime(Sys.time(), t_start, units="mins")

# Fix k and a to MLEs
ig_mles <- get_ig_mle(Y, dt=dt)
logk_hat <- ig_mles$log_k
loga_hat <- ig_mles$log_a
init_param2 <- list(log_k = logk_hat, 
		    log_a = loga_hat, 
		    Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
		    x = prop_paths(Y, dt, log_k=logk_hat, log_a=loga_hat))
adfun2 <- TMB::MakeADFun(data=list(model="factor_model", n_factor=n_factor,
                                        dt=dt, Y=Y, lam=1, nu=15.),
                                  parameters=init_param2,
				  map = list(log_k = as.factor(rep(NA, n_cell)),
					     log_a = as.factor(rep(NA, n_cell))),
                                  random = "x",
                                  DLL = "fastr_TMBExports",
                                  silent = F)

t_start <- Sys.time()
fit_tmb2 <- nlminb(adfun2$par, adfun2$fn, adfun2$gr)
rep2 <- sdreport(adfun2)
time_tmb2 <- difftime(Sys.time(), t_start, units="mins")

############# Competing method 1: simple k-means #########################
# Added Gaussian noise
Y.noisy <- Y + array(rnorm(length(Y), 0, 0.1), dim(Y))
km.res <- kmeans(Y.noisy[,,1], center=5 )$cluster
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


############# Factor model: applying k-means on the loading matrix ############
res <- readRDS("TH2015results_cluster.rds")
fit <- res$fit
n_cell <- 43
n_factor <- 6
devtools::load_all("~/projects/NeuronModel/mnfa")
L <- varimax(get_FA_estim(fit, n_cell=n_cell, n_factor=n_factor)$L)$loadings[1:n_cell,]
fa.res <- kmeans(L, centers=5)$cluster

############# Plotting #############################
n_cell <- 43
n_factor <- 6
true.group <- c(rep(1, 2), rep(2,8), rep(3,15), rep(4,8), rep(5,10))
par(mfrow=c(1,4))
image(x=1, y=1:n_cell, matrix(true.group, ncol=n_cell), ylab="Neuron index", 
      xlab="", xaxt="n", yaxt="n",
      main="True group memberships")
axis(side=2, at=seq(1,n_cell,by=1))
image(x=1, y=1:n_cell, matrix(km.res, ncol=n_cell), ylab="", 
      xlab="", xaxt="n", yaxt="n",
      main="K-means assigned group memberships")
axis(side=2, at=seq(1,n_cell,by=1))
image(x=1, y=1:n_cell, matrix(fa.res, ncol=n_cell), ylab="", 
      xlab="", xaxt="n", yaxt="n",
      main="K-means applied on factor loadings")
axis(side=2, at=seq(1,n_cell,by=1))



