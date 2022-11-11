library(ggplot2)
library(fields)
library(heatmaply)
devtools::load_all()
set.seed(123)
save_path <- "~/Thesis_proposal/"
####### Simulate data #############
dt <- 0.005
n_bin <- 2000
n_cell <- 20
n_factor <- 4
n_trial <- 1
k <- runif(n_cell, 0.1, 0.6)
alpha <- runif(n_cell, 2, 5.5)
l1 <- runif(n_cell, 0.1, 0.25)
l2 <- runif(n_cell, 0.1, 0.25)
l3 <- runif(n_cell, 0.1, 0.25)
l4 <- runif(n_cell, 0.1, 0.25)
l1[1:5] <- runif(5, 0.7, 0.95)
l2[6:10] <- runif(5, 0.7, 0.95)
l3[11:15] <- runif(5, 0.7, 0.95)
l4[16:20] <- runif(5, 0.7, 0.95)
L <- cbind(l1, l2, l3, l4)
true_cor <- L%*%t(L)
Psi <- 1-diag(true_cor) 
diag(true_cor) <- 1
sim <- simdata(dt=dt, n_bin=n_bin, n_trial=n_trial, alpha=alpha, k=k, L=L)
Y <- sim$Y

cat("The simulation parameters are:\n")
cat("alpha:\n")
cat(paste(round(alpha,2), collapse=", "), "\n")
cat("k:\n")
cat(paste(round(k,2), collapse=", "), "\n")
cat("Lambda:\n")
cat(paste(round(t(L)[1,],2), collapse=" & "), "\n")
cat(paste(round(t(L)[2,],2), collapse=" & "), "\n")
cat(paste(round(t(L)[3,],2), collapse=" & "), "\n")
cat(paste(round(t(L)[4,],2), collapse=" & "), "\n")
cat("Psi:\n")
cat(paste(round(Psi,2), collapse=", "), "\n")

####### Choose n_factor #########
# Fit a single neuron model to each of the neurons
d_max <- 10
d_seq <- 1:d_max
d_res <- choose_n_factor(Y, dt, d_seq, plot=FALSE)
cumvars <- d_res$cumvars
pdf(paste0(save_path, "app-choose-d.pdf"), width=12, height=4)
par(mfrow=c(1,3))
plot(d_seq, cumvars, type="l", ylab="Cumvar(d)",
     xlab="Number of factors (d)", 
     cex.axis=1.4, cex.lab=1.5)
abline(v=n_factor, lty="dashed", col="green", lwd=3)
plot(d_seq[1]:tail(d_seq,2)[1], diff(cumvars), type="l",
     ylab="Cumvar(d+1) - Cumvar(d)",
     xlab="Number of factors (d)", 
     cex.axis=1.4, cex.lab=1.5)
abline(v=n_factor, lty="dashed", col="green", lwd=3)
plot(d_seq, d_res$minmaxload, type="l",
     ylab="Min of max loadings per factor",
     xlab="Number of factors (d)", 
     cex.axis=1.4, cex.lab=1.5)
abline(v=n_factor, lty="dashed", col="green", lwd=3)
dev.off()

####### Model fitting #############
fit_2step <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, method="2step")
loga_hat <- fit_2step$log_a_hat
logk_hat <- fit_2step$log_k_hat
loga_se <- fit_2step$loga_se
logk_se <- fit_2step$logk_se
rate_hat <- fit_2step$rate_hat
rate_se <- fit_2step$se_ig[1:n_cell]
cor4 <- (fit_2step$lmat_hat) %*% t(fit_2step$lmat_hat)
diag(cor4) <- 1
err4 <- cor4[lower.tri(cor4)] - true_cor[lower.tri(true_cor)]

run_joint_fit <- FALSE
if (run_joint_fit){
	logk_hat <- fit_2step$log_k_hat
	loga_hat <- fit_2step$log_a_hat
	fit_joint <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, 
			       init=list(log_k=logk_hat, 
					 log_a=loga_hat, 
					 Lt=rep(1,
					 n_cell*n_factor-n_factor*(n_factor-1)/2)),
			       method="joint")
}
# With 6 neurons, 2 factors, 2000 bins, and 5 trials, 
# the joint method is 6.8 times slower than the 2-step method even with 
# initial values of a and k set to MLEs.

####### Accuracy check #############
loga_lim <- range(log(alpha), loga_hat+2*loga_se, loga_hat-2*loga_se)
logk_lim <- range(log(k), logk_hat+2*logk_se, logk_hat-2*logk_se)
rate_lim <- range(alpha/k, rate_hat+2*rate_se, 
                  rate_hat-2*rate_se)
ggplot(mapping=aes(x=log(alpha), y=loga_hat)) + 
  xlim(loga_lim)+ylim(loga_lim)+
  geom_point()+
  geom_errorbar(aes(ymin=loga_hat-2*loga_se, ymax=loga_hat+2*loga_se), width=.02)+
  geom_abline(slope=1, linetype="dashed", color="red")+
  labs(x=expression(log(alpha)), y=expression(widehat(log(alpha))))+
  theme_bw()
ggsave(paste0(save_path, "true-vs-est-loga.pdf"), width=4, height=4)

ggplot(mapping=aes(x=log(k), y=logk_hat)) + 
  xlim(logk_lim)+ylim(logk_lim)+
  geom_point()+
  geom_errorbar(aes(ymin=logk_hat-2*logk_se, ymax=logk_hat+2*logk_se), width=.02)+
  geom_abline(slope=1, linetype="dashed", color="red")+
  labs(x=expression(log(k)), y=expression(widehat(log(k))))+
  theme_bw()
ggsave(paste0(save_path, "true-vs-est-logk.pdf"), width=4, height=4)

ggplot(mapping=aes(x=alpha/k, y=rate_hat)) + 
  xlim(rate_lim)+ylim(rate_lim)+
  geom_point()+
  geom_errorbar(aes(ymin=rate_hat-2*rate_se, ymax=rate_hat+2*rate_se), width=.02)+
  geom_abline(slope=1, linetype="dashed", color="red")+
  labs(x=expression(alpha/k), y=expression(widehat(alpha/k)))+
  theme_bw()
ggsave(paste0(save_path, "true-vs-est-rate.pdf"), width=4, height=4)


pdf(paste0(save_path, "True-Lambda.pdf"), width=3.4, height=7)
par(mar=c(2,4,0.5,1))
image(x=1:n_factor, y=1:n_cell,
       z=t(L), zlim=c(-1,1), col=viridis::viridis(20),
       ylab="Neuron index", xlab="", axes=FALSE)
axis(1, at=1:n_factor, tick=FALSE)
axis(2, at=1:n_cell, labels=1:n_cell, las=2, tick=FALSE)
dev.off()
pdf(paste0(save_path, "Est-Lambda.pdf"), width=3.4, height=7)
par(mar=c(2,4,0.5,1))
plot(fit_2step, legend=F)
dev.off()

####### Check fitting with more or less factors #########
fit2 <- fastr_fit(data=Y, dt=dt, n_factor=2, method="2step")
fit6 <- fastr_fit(data=Y, dt=dt, n_factor=6, method="2step")
# Uniqueness of neurons 1-10 in 2-factor model
Lam2 <- fit2$lmat_varimax
uniq2 <- round((1-diag(Lam2 %*% t(Lam2)))[1:10], 2)
# Estimated corrlation matrix
cor2 <- (fit2$lmat_hat) %*% t(fit2$lmat_hat)
diag(cor2) <- 1
err2 <- cor2[lower.tri(cor2)] - true_cor[lower.tri(true_cor)]
cor6 <- (fit6$lmat_hat) %*% t(fit6$lmat_hat)
diag(cor6) <- 1
err6 <- cor6[lower.tri(cor6)] - true_cor[lower.tri(true_cor)]
pdf(paste0(save_path, "sim-boxplot-err.pdf"), width=5, height=3)
par(mar=c(4, 3, 0.5, 1))
boxplot(err2, err4, err6, xlab="Number of factors", 
	names=c("2", "4", "6"),
        ylim=c(-0.6,0.4), outcex=0.5)
abline(h=0, lty="dashed", col="red")
dev.off()

cat(paste(uniq2, collapse=", "), "\n")
pdf(paste0(save_path, "Est-Lambda-2factors.pdf"), width=2.5, height=7)
par(mar=c(2,4,0.5,1))
plot(fit2, legend=F)
dev.off()
pdf(paste0(save_path, "Est-Lambda-6factors.pdf"), width=5, height=7)
par(mar=c(2,4,0.5,1))
plot(fit6)
dev.off()

############# Compare with simple correlation plot on Y ##################
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
Y.conv <- apply(Y, c(1,3), expdec.conv, tau=5)
get_lag_cc <- function(x, y, lag=0, nonsig=0){
  n <- length(x)
  cor <- ccf(x,y,lag.max=lag,plot=F)$acf[[1]]
  sig <- 2/sqrt(n)
  if (abs(cor) > abs(sig)) return(cor)
  else return(nonsig)
}
Y_cor <- matrix(1, nrow=n_cell, ncol=n_cell)
for (i in 1:n_cell){
  for (j in 1:(n_cell-1)){
    Y_cor[i,j] <- get_lag_cc(Y.conv[,i,1], Y.conv[,j,1])
    Y_cor[j,i] <- Y_cor[i,j]
  }
}
colnames(Y_cor) <- rownames(Y_cor) <- 1:n_cell
pdf(paste0(save_path, "convolved-Y.pdf"), height=3, width=9)
par(mar=c(4,4.5,0.1,0.1))
plot(Y.conv[1:1000,1,1], type="l", ylab="Convolved Y", xlab="Time (bin index)", 
     cex.lab=1.5, cex.axis=1.5)
points(which(Y[1,1:1000,1]==1), rep(0,sum(Y[1,1:1000,1])), pch="l", col="red")
dev.off()

pdf(paste0(save_path, "correlation-plot.pdf"), height=4, width=4.5)
par(mar=c(0.1,3,3, 5))
corrplot(Y_cor)
dev.off()

heatmaply_cor(Y_cor)

