library(mvtnorm) # for rmvnorm
library(fiels) # for image.plot
res <- readRDS("test-bigdata-results-40neurons.rds")
n_cell <- 40
alpha <- res$true_a
k <- res$true_k
fi <- alpha/k
L <- res$true_L
est <- res$estimates
est_L <- est$L
varimax_L <- matrix(as.numeric(varimax(est_L)$loadings), ncol=4, nrow=40)
est_fi <- est$fi
rep <- res$TMB_rep
cov <- rep$cov.fixed
mean <- rep$par.fixed
n_sam <- 10000
pos_sam <- rmvnorm(n_sam, mean, cov)
pos_mean <- apply(pos_sam, 2, mean)
pos_ci <- apply(pos_sam, 2, quantile, probs=c(0.025, 0.975)) 

pos_sam_fi <- apply(pos_sam, 1, function(x){exp(x[(n_cell+1):(2*n_cell)]-x[1:n_cell])})
pos_mean_fi <- apply(pos_sam_fi, 1, mean)
pos_ci_fi <- apply(pos_sam_fi, 1, quantile, probs=c(0.025, 0.975))

pdf("sim-firing.pdf",width=6, height=6)
par(mar=c(5.1, 4.8, 1, 2.1))
plot(fi, est_fi, xlab=expression(True~alpha/k), ylab=expression(widehat(alpha)/widehat(k)), pch=16)
abline(0, 1, lty="dashed", col="grey")
lines(fi[order(fi)], pos_ci_fi[1,order(fi)], lty="dashed", col="red", cex=0.8)
lines(fi[order(fi)], pos_ci_fi[2,order(fi)], lty="dashed", col="red", cex=0.8)
legend("topleft", pch=c(16, NA, NA), lty=c(NA, "dashed", "dashed"),
       col=c("black", "red", "grey"),
       legend=c("Posterior mean for one neuron", "95% credible interval", "45 degree line"))
dev.off()

pdf("sim-loading-mat.pdf", width=9, height=6)
par(mfrow=c(1,2), mar=c(5, 5, 5, 5))
fields::image.plot(1:4, 1:40,t(true_L), zlim=c(-1,1), 
		   xlab="Column (factor)", ylab="Row (Neuron)", 
		   main=expression(Visualization~of~true~loading~matrix~Lambda))
fields::image.plot(1:4, 1:40,t(varimax_L), zlim=c(-1,1),
		   xlab="Column (factor)", ylab="Row (Neuron)", 
		   main=expression(Visualization~of~Varimax-transformed~widehat(Lambda)))
dev.off()
