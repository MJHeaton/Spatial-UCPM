#library(tidyverse)
library(rmutil)
#library(invgamma)
library(MASS)
library(parallel)

crash <- read.csv("UCPM-UCSMinput_New.csv", header = TRUE)
all.roads <- unique(crash$ROUTE_ID)
rm.roads <- NULL
for(i in 1:length(all.roads)){
  data <- crash[crash$Route_Name == all.roads[i],] 
  if(nrow(data) == 1) rm.roads <- c(rm.roads,data$Route_Name)
}
crash <- crash[!crash$Route_Name %in% rm.roads,]
centroid <- rowMeans(cbind(crash$BEG_MILEPOINT,crash$END_MILEPOINT))
crash$seg_length <- crash$END_MILEPOINT - crash$BEG_MILEPOINT

#######################
## Combine P/N Roads ##
#######################
dup <- which(duplicated(crash[,c("BEG_MILEPOINT", "END_MILEPOINT", "ROUTE_ID")]))
for(d in dup){
  same.seg <- which(crash[d,"BEG_MILEPOINT"]==crash[-d,"BEG_MILEPOINT"] &
                      crash[d,"END_MILEPOINT"]==crash[-d,"END_MILEPOINT"] &
                      crash[d,"ROUTE_ID"]==crash[-d,"ROUTE_ID"])
  crash[same.seg,39:ncol(crash)] <- colSums(crash[c(same.seg, d), 39:ncol(crash)])
}
crash <- crash[-dup,]

########################################
## parallel processing implementation ##
########################################

# setting up priors and starting info
new.crash <- crash
new.crash$VMT_cat <- cut(new.crash$VMT, breaks=quantile(new.crash$VMT, seq(0,1,by=0.1)),
                         include.lowest=TRUE)
new.crash$Lanes_cat <- cut(new.crash$Num_Lanes, breaks=c(-1,2.1,5.1,8.1))
# new.crash %>% group_by(Lanes_cat, VMT_cat) %>% dplyr::summarize(mean=mean(Total_Crashes))
# ggplot(data=new.crash, aes(x=Lanes_cat:VMT_cat, y=Total_Crashes)) + geom_boxplot()
n.cores <- 10#detectCores()
n.draws <- 2500#3500
burn <- 1000000
thin <- 2000
tot.it <- burn+n.draws*thin
kp.seq <- seq(burn+thin,tot.it, by=thin)
x.mat <- model.matrix(Total_Crashes ~ Total_Percent_Trucks + Lanes_cat + SPEED_LIMIT + VMT_cat, data = new.crash)
#tst <- glm.nb(new.crash$Total_Crashes~x.mat-1)
#plot(fitted(tst), new.crash$Total_Crashes, pch=19)
#x.mat <- matrix(1, nrow=nrow(new.crash))

cur.beta <- coef(glm(Total_Crashes~x.mat-1, data=new.crash, family=poisson))
if("Crash_New_Results.RData"%in%list.files()){
  load('Crash_New_Results.RData')
  cur.beta <- rowMeans(all.beta)
}
roads <- unique(new.crash$ROUTE_ID)
crashes <- vector("list",length(roads))

#cur.beta <- coef(glm.nb(Total_Crashes~x.mat-1, data=new.crash))

for(i in 1:length(roads)){
  data <- new.crash[new.crash$Route_Name == roads[i],] 
  if(nrow(data) == 1) next
  X <- matrix(x.mat[new.crash$Route_Name == roads[i],], ncol=ncol(x.mat))
  n <- nrow(data)
  A <- matrix(NA, nrow = n, ncol = n)
  for(j in 1:n){
    A[j,] <- ifelse((data$BEG_MILEPOINT == data$END_MILEPOINT[j] | data$END_MILEPOINT == data$BEG_MILEPOINT[j]),1,0)
  }
  
  ## Hughes and Haran spatial bases
  ones <- matrix(1, nrow=nrow(X), ncol=1)
  # P <- X%*%chol2inv(chol((t(X)%*%X+.001*diag(ncol(X)))))%*%t(X)
  P <- ones%*%chol2inv(chol((t(ones)%*%ones)))%*%t(ones)
  p.orth <- diag(nrow(P)) - P
  E <- eigen(p.orth%*%A%*%p.orth, symmetric=TRUE)
  #E$vectors <- E$vectors[,order(abs(E$values), decreasing=TRUE)]
  #E$values <- sort(abs(E$values), decreasing=TRUE)
  n90 <- sum(cumsum(E$values[E$values>0])/sum(E$values[E$values>0]) < 1)
  if(n90 == 0){
    M <- ones
    y.adj <- log(data$Total_Crashes+0.5)-X%*%cur.beta
    theta <- mean(y.adj)
  } else {
     M <- cbind(1,matrix(E$vectors[,1:n90], ncol=n90))
     y.adj <- log(data$Total_Crashes+0.5)-X%*%cur.beta
     theta <- coef(lm(y.adj~M-1))
  }
  n.coef <- n90 + 1 + 1
  #plot(X%*%cur.beta+M%*%theta, log(data$Total_Crashes+0.5)) 
  
  ## Hodges and Reich RSR
  # P <- matrix(1, nrow=nrow(X), ncol=1)
  # P <- P%*%t(P)/n
  # p.orth <- diag(n)-P
  # E <- eigen(p.orth, symmetric=TRUE)
  # n90 <- n-1
  # n.coef <- n90 + 1 + 1
  
  ## filling in crash data for all roads
  crashes[[i]]$y <- data$Total_Crashes
  crashes[[i]]$X <- X
  crashes[[i]]$M <- M
  crashes[[i]]$amcmc <- list(mn=matrix(0,ncol=1,nrow=n.coef),
                             var=matrix(0,nrow=n.coef,ncol=n.coef))
  crashes[[i]]$theta.draws <- matrix(0,nrow=ncol(crashes[[i]]$M),ncol=n.draws)
  crashes[[i]]$theta <- theta
  crashes[[i]]$logr.draws <- rep(0,n.draws)
  crashes[[i]]$logr <- log(100)
  crashes[[i]]$kp <- 0
}
# r <- sample(length(crashes), 1)
# plot(exp(crashes[[r]]$X%*%cur.beta+crashes[[r]]$M%*%crashes[[r]]$theta), 
#      crashes[[r]]$y, pch=19)

## update theta and log r function
update.theta.logr <- function(x){
  if(i > 500){
    Sigma <- (2.4^2/(ncol(x$M + 1)))*(eps*diag(ncol(x$M) + 1)+x$amcmc$var)
  } else {
    Sigma <- eps*diag(ncol(x$M )+ 1)
  }
  ## use amcmc to propose new theta and log r
  params.star <- c(x$theta,x$logr) + t(chol(Sigma))%*%rnorm(length(x$theta) + 1)
  n.params <- length(params.star)
  # MH.ratio.num <- sum(dnbinom(x$y,mu = c(exp(x$X%*%cur.beta + x$M%*%params.star[-n.params])), size = exp(params.star[n.params]), log = TRUE)) + 
  #   sum(dnorm(t(params.star), prior.theta.logr, sd= prior.theta.logr.sd,log=TRUE)) 
  # MH.ratio.denom <- sum(dnbinom(x$y,mu = c(exp(x$X%*%cur.beta + x$M%*%x$theta)), size = exp(x$logr), log = TRUE)) + 
  #   sum(dnorm(t(c(x$theta,x$logr)), prior.theta.logr, sd= prior.theta.logr.sd,log=TRUE))
  MH.ratio.num <- sum(dpois(x$y,lambda = c(exp(x$X%*%cur.beta + x$M%*%params.star[-n.params])), log = TRUE)) +
    sum(dnorm(t(params.star), prior.theta.logr, sd= prior.theta.logr.sd,log=TRUE))
  MH.ratio.denom <- sum(dpois(x$y,lambda = c(exp(x$X%*%cur.beta + x$M%*%x$theta)), log = TRUE)) +
    sum(dnorm(t(c(x$theta,x$logr)), prior.theta.logr, sd= prior.theta.logr.sd,log=TRUE))
  
  if(log(runif(1)) < (MH.ratio.num-MH.ratio.denom)){
    x$theta <- params.star[-n.params]
    x$logr <- params.star[n.params]
  }
  
  if(i%in%kp.seq){
    x$kp <- x$kp+1
    x$theta.draws[,x$kp] <- x$theta
    x$logr.draws[x$kp] <- x$logr
  }
  
  x$amcmc <- AMCMC.update(matrix(c(x$theta,x$logr),ncol=1),x$amcmc$mn,x$amcmc$var,i)
  return(x)
}
#lapply(crashes,update.theta.logr)

## likelihood for beta across all M and theta
get.likelihood <- function(x,beta.coef){
  # llike <- sum(dnbinom(x$y,mu=exp(x$X%*%beta.coef+x$M%*%x$theta),size=exp(x$logr),log=TRUE))
  llike <- sum(dpois(x$y,lambda=exp(x$X%*%beta.coef+x$M%*%x$theta),log=TRUE))
  llike
}
#sapply(crashes,FUN=get.likelihood,beta.coef=rep(5,ncol(x)))


## parallel AMCMC
## Prior information ##
prior.beta <- 0 
prior.beta.sd <- 3 

prior.theta.logr <- 1 
prior.theta.logr.sd <- 3


# initializing values
#n.draws <- 150000
P <- ncol(x.mat)
all.beta <- matrix(0,nrow=P,ncol=n.draws)
#cur.beta <- c(1.23,-1.35,0.68,0.01,0.38) #
#cur.beta <- c(log(mean(sapply(crashes,function(x){mean(x$y)}))),rep(0,P-1)) 
tau.draws <- rep(0,n.draws)
tau <- 2
a <- 2.01
b <- 1
n.beta <- nrow(all.beta)
eps <- 0.00001^2#0.001^2 #for NB

source("AMCMCUpdate.R")
kp <- 0
amcmc <- list(mn=matrix(0,ncol=1,nrow=ncol(x.mat)),
              var=matrix(0,nrow=ncol(x.mat),ncol=ncol(x.mat)))
pb <- txtProgressBar(min = 0, max = tot.it, style = 3)
for(i in 2:tot.it){
  
  # step 1 - update beta
  if(i > 500){
    beta.Sigma <- (2.4^2/(P))*(eps*diag(ncol(x.mat))+amcmc$var)
  } else {
    beta.Sigma <- eps*diag(ncol(x.mat))
  }
  beta.star <- cur.beta + t(chol(beta.Sigma))%*%rnorm(length(cur.beta))
  MH.ratio.num <- sum(sapply(crashes,FUN=get.likelihood,beta.coef=beta.star)) + 
    sum(dlaplace(t(beta.star), prior.beta, tau, log = TRUE))# sum(dnorm(t(beta.star), prior.beta, sd= prior.beta.sd,log=TRUE)) 
  MH.ratio.denom <- sum(sapply(crashes,FUN=get.likelihood,beta.coef=cur.beta)) + 
    sum(dlaplace(t(cur.beta), prior.beta, tau, log = TRUE))# sum(dnorm(t(beta.cur), prior.beta, sd= prior.beta.sd,log=TRUE))
  #(MH.ratio.num-MH.ratio.denom)
  if(log(runif(1)) < (MH.ratio.num-MH.ratio.denom)){
    #all.beta[,i] <- beta.star
    cur.beta <- beta.star
  }
  amcmc <- AMCMC.update(matrix(cur.beta,ncol=1),amcmc$mn,amcmc$var,i)
  
  # step 2 - updating theta for all roads
  crashes <- mclapply(crashes, update.theta.logr, mc.cores = n.cores)
  
  #step 3 - updating tau
  tau <- 1/rgamma(1, a + n.beta, rate = b + sum(abs(cur.beta - prior.beta)))
  
  ## Keep beta, tau if necessary
  if(i%in%kp.seq){
    kp <- kp+1
    all.beta[,kp] <- cur.beta
    tau.draws[kp] <- tau
  }
  
  #cat(i,"\r")
  setTxtProgressBar(pb, i)
  
  # if(i == 25000) save(crashes,all.beta,tau.draws,file = "q1results.RData")
  # if(i == 50000) save(crashes,all.beta,tau.draws,file = "mid_results.RData")
  # if(i == 75000) save(crashes,all.beta,tau.draws,file = "q3results.RData")
}
close(pb)

## Trace plots
# plot(all.beta[sample(nrow(all.beta),1),-ncol(all.beta)], type="l")
# r <- sample(length(crashes), 1)
# s <- sample(nrow(crashes[[r]]$theta.draws), 1)
# plot(crashes[[r]]$theta.draws[s,], type="l")
# 
# ## Check X%*%beta results
# qplot(rowMeans(exp(x.mat%*%all.beta)), geom="histogram")
# qplot(rowMeans(exp(x.mat%*%all.beta)), new.crash$Total_Crashes, geom="point")
# fitVObs <- sapply(crashes, function(x){
#   cbind(exp(x$X%*%cur.beta+x$M%*%x$theta), x$y, x$X)
# })
# fit <- sapply(crashes, function(x){
#   c(exp(x$X%*%cur.beta+x$M%*%x$theta))
# })
# obs <- sapply(crashes, function(x){
#   x$y
# })
# plot(unlist(fit), unlist(obs))
# abline(a=0,b=1)
# cor(unlist(fit), unlist(obs))^2
# maxs <- sapply(crashes, function(x){
#   max(exp(x$X%*%cur.beta+x$M%*%x$theta))
# })

save(crashes,all.beta,tau.draws,file = "Crash_New_Results_No_Dup.RData")
