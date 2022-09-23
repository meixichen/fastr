library(ggplot2)
library(fields)
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


####### Model fitting #############
fit_2step <- fastr_fit(data=Y, dt=dt, n_factor=n_factor, method="2step")
loga_hat <- fit_2step$log_a_hat
logk_hat <- fit_2step$log_k_hat
loga_se <- fit_2step$loga_se
logk_se <- fit_2step$logk_se
lograte_hat <- loga_hat-logk_hat
lograte_se <- sqrt(loga_se^2+logk_se^2) # <------- need to change

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

####### Check fitting with more or less factors #########
fit2 <- fastr_fit(data=Y, dt=dt, n_factor=2, method="2step")
fit6 <- fastr_fit(data=Y, dt=dt, n_factor=6, method="2step")
pdf(paste0(save_path, "Est-Lambda-2factors.pdf"), width=3, height=7)
par(mar=c(2,4,0.5,1))
plot(fit2)
dev.off()
pdf(paste0(save_path, "Est-Lambda-6factors.pdf"), width=5, height=7)
par(mar=c(2,4,0.5,1))
plot(fit6)
dev.off()

####### Accuracy check #############
loga_lim <- range(log(alpha), loga_hat+2*loga_se, loga_hat-2*loga_se)
logk_lim <- range(log(k), logk_hat+2*logk_se, logk_hat-2*logk_se)
lograte_lim <- range(log(alpha)-log(k), lograte_hat+2*lograte_se, lograte_hat-2*lograte_se)
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

ggplot(mapping=aes(x=log(alpha)-log(k), y=lograte_hat)) + 
  xlim(lograte_lim)+ylim(lograte_lim)+
  geom_point()+
  geom_errorbar(aes(ymin=lograte_hat-2*lograte_se, ymax=lograte_hat+2*lograte_se), width=.02)+
  geom_abline(slope=1, linetype="dashed", color="red")+
  labs(x=expression(log(alpha/k)), y=expression(widehat(log(alpha/k))))+
  theme_bw()
ggsave(paste0(save_path, "true-vs-est-lograte.pdf"), width=4, height=4)

pdf(paste0(save_path, "True-Lambda.pdf"), width=4, height=7)
par(mar=c(2,4,0.5,1))
fields::image.plot(x=1:n_factor, y=1:n_cell,
       z=t(L), zlim=c(-1,1),
       ylab="Neuron index", xlab="", axes=FALSE)
axis(1, at=1:n_factor, tick=FALSE)
axis(2, at=1:n_cell, labels=1:n_cell, las=2, tick=FALSE)
dev.off()
pdf(paste0(save_path, "Est-Lambda.pdf"), width=4, height=7)
par(mar=c(2,4,0.5,1))
plot(fit_2step)
dev.off()

############# Compare with simple correlation plot on Y ##################
neuron_index <- sample(1:n_cell)
Y_cor <- cor(t(Y[neuron_index,,1]))
colnames(Y_cor) <- rownames(Y_cor) <- neuron_index
corrplot(Y_cor, order="hclust")
# Comment: the correlation at the data level is already pretty strong,
# so clusters can be identified by simple correlation.
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



