# Study the sensitivity of the model to the number of factors specified
require(TMB)
set.seed(123)

#------------ Simulate data ---------------
dt <- 0.005
n_bin <- 2000
n_trial <- 5
n_cell <- 12

alpha <- runif(n_cell, 2, 5.5)
k <- runif(n_cell, 0.1, 0.6)
# Suppose the "true" number of factors is 4
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l3 <- runif(n_cell, 0.1, 0.25)
l4 <- runif(n_cell, 0.1, 0.25)
cluster_size <- n_cell/4
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

#------------ Model fitting ---------------
n_factor_list <- c(2, 3, 4, 5, 6)
#all_mle <- get_ig_mle(Y, dt)
#log_k_mle <- log(all_mle$lam)/2
#log_a_mle <- log_k_mle - log(all_mle$mu)
out <- list()
for (i in 1:length(n_factor_list)){
  n_factor <- n_factor_list[i]
  init_param <- list(log_k = rep(-1,n_cell),
		     log_a = rep(0,n_cell),
		     Lt = rep(1, n_cell*n_factor-n_factor*(n_factor-1)/2),
		     x = prop_paths(Y, dt, rep(-1,n_cell), rep(0,n_cell)))
  data <- list(model="factor_model_eff", n_factor=n_factor, 
	       dt=dt, Y=Y, lam=0.5, nu=15.)
  adfun<- TMB::MakeADFun(data=data,
			 parameters=init_param,
			 random = "x",
			 DLL = "mnfa_TMBExports",
			 silent = F)

  t_start <- Sys.time()
  fit <- nlminb(adfun$par, adfun$fn, adfun$gr, 
                control = list(eval.max=1000, iter.max=1000))
  rep <- sdreport(adfun)
  t_taken <- difftime(Sys.time(), t_start, units="mins")
  out_i <- list(fit=fit, rep=rep, time=t_taken, 
		 L=varimax(get_FA_estim(fit, n_cell, n_factor)$L))
  saveRDS(out_i, file = paste0("sensitivity-analysis-factor", n_factor, ".rds"))
  out[[i]] <- out_i
}
saveRDS(out, file="sensitivity-analysis-results.rds")

#-------------- Plots ------------------
fields::image.plot(1:4, 1:n_cell,t(L), zlim=c(-1,1),
                   xlab="Column (factor)", ylab="Row (Neuron)",
                   main=expression(Visualization~of~true~loading~matrix~Lambda))
fields::image.plot(1:2, 1:n_cell,t(out[[1]]$L$loadings[1:n_cell,]), zlim=c(-1,1),
                   xlab="Column (factor)", ylab="Row (Neuron)",
                   main=expression(Visualization~of~Varimax-transformed~widehat(Lambda)~with~2~factors))
fields::image.plot(1:3, 1:n_cell,t(out[[2]]$L$loadings[1:n_cell,]), zlim=c(-1,1),
                   xlab="Column (factor)", ylab="Row (Neuron)",
                   main=expression(Visualization~of~Varimax-transformed~widehat(Lambda)~with~3~factors))
fields::image.plot(1:4, 1:n_cell,t(out[[3]]$L$loadings[1:n_cell,]), zlim=c(-1,1),
                   xlab="Column (factor)", ylab="Row (Neuron)",
                   main=expression(Visualization~of~Varimax-transformed~widehat(Lambda)~with~4~factors))
fields::image.plot(1:5, 1:n_cell,t(out[[4]]$L$loadings[1:n_cell,]), zlim=c(-1,1),
                   xlab="Column (factor)", ylab="Row (Neuron)",
                   main=expression(Visualization~of~Varimax-transformed~widehat(Lambda)~with~5~factors))
fields::image.plot(1:6, 1:n_cell,t(out[[5]]$L$loadings[1:n_cell,]), zlim=c(-1,1),
                   xlab="Column (factor)", ylab="Row (Neuron)",
                   main=expression(Visualization~of~Varimax-transformed~widehat(Lambda)~with~6~factors))


par(mfrow=c(2,3))
for (i in 1:5){
  plot(k, exp(out[[i]]$fit$par[1:n_cell]), xlab="k", ylab=expression(widehat(k)))
  abline(0, 1, lty="dashed", col="blue")
  plot(alpha, exp(out[[i]]$fit$par[(1+n_cell):(2*n_cell)]), 
       xlab=expression(alpha), ylab=expression(widehat(alpha)))
  abline(0, 1, lty="dashed", col="blue")
}
