library(rvalues)
library(scales)
library(lsei)
library(jmuOutlier)
library(parallel)
library(Rmosek)
library(REBayes)

N <- 100 #orginal sample size
set.seed(1111) #cts
#set.seed(1) #discrete
M=100
iter=10000

#mus <- rnorm(N)
#mus <- rgamma(N,5,2)
#sds <- 0.1*rep(1,N)
#samples <- rnorm(n=N,mean=mus,sd=sds)
#ran = c(-3,3)
ran = c(0.5,6)


#mus <- rnorm(N)
mus <- rgamma(N,5,2)
sds <- 0.1*rep(1,N)
samples <- rnorm(n=N,mean=mus,sd=sds)

np = GLmix(samples, sigma = sds[1])



SSS=cbind(samples,rep(0.1,N))
kapprox1= npmle(SSS, family = gaussian, maxiter =100)
kapprox1$support = np$x
kapprox1$mix.prop = np$y

plot(kapprox1)

modi_npmle <- function(kapprox,eps=10^(-4)){
  kapprox$support=kapprox$support[kapprox1$mix.prop>eps]
  kapprox$mix.prop=kapprox$mix.prop[kapprox1$mix.prop>eps]
  kapprox$mix.prop=kapprox$mix.prop/sum(kapprox$mix.prop)
  
  kapprox$Fhat=stepfun(kapprox$support, c(0, cumsum(kapprox$mix.prop)))
  kapprox$fhat <- density(kapprox$support, weights = kapprox$mix.prop)
  kapprox$fhat <- approxfun(kapprox$fhat$x, kapprox$fhat$y)
  return(kapprox)
}

kapprox = modi_npmle(kapprox1)

plot(kapprox,col="red",main="",xlim = ran)





sgd_boot_x_svrg = function(support,w, sd, bsize){
  components <- sample(1:length(w), prob=w, size=bsize, replace=TRUE)
  #print(components)
  x = rnorm(n=bsize,mean=support[components],sd=sd)
  return(x)
}

sgd_boot_w_svrg = function(support,w, sd, x_new,eta){
  W = w*0
  for (i in 1:length(x_new)){
    new_w = w*0
    comp = dnorm(x_new[i], mean=support, sd=sd)
    new_w = (1*comp)*(100*w)
    new_w = new_w/sum(new_w)
    W = W+new_w
  }
  return(eta*w/sum(w)+(1-eta)*1/length(x_new)*W)
}


sgd_boot_theta_svrg = function(support,w, sd, x_new,eta){
  W = w*0
  for (i in 1:length(x_new)){
    new_w = w*0
    
    
    comp = dnorm(x_new[i], mean=support, sd=sd)
    grad = dnorm(x_new[i], mean=support, sd=sd)*(x_new-support)/(sd^2)
    
    new_w = comp*(100*w)
    new_w = (100*w)/sum(new_w)
    
    W = new_w*grad
  }
  return(support+eta*W)
}

modi_npmle_boot <- function(kapprox,support,w){
  kapprox$support=support
  kapprox$mix.prop=w
  kapprox$mix.prop=w/sum(w)
  
  kapprox$Fhat=stepfun(kapprox$support, c(0, cumsum(kapprox$mix.prop)))
  kapprox$fhat <- density(kapprox$support, weights = kapprox$mix.prop)
  kapprox$fhat <- approxfun(kapprox$fhat$x, kapprox$fhat$y)
  return(kapprox)
}


sd = 0.1
bsize=1

para_boot <- function(M,kapprox,sd,bsiz){
  support = kapprox$support
  w = kapprox$mix.prop
  eta=0.9
  
  for (i in 1:iter){
    eta = 1-1/(N+i)
    x = sgd_boot_x_svrg(support,w,sd,bsize)
    ww=w
    w = sgd_boot_w_svrg(support,w,sd,x,eta)
    support = sgd_boot_theta_svrg(support,ww, sd, x,1-eta)
  }
  support=sort(support)
  w=w[order(support)]
  return( modi_npmle_boot(kapprox,support,w))
  
}




clnum<-detectCores() 

cl <- makeCluster(getOption("cl.cores", clnum-1));

clusterExport(cl, c('para_boot',"iter","N","sgd_boot_x_svrg","bsize","sgd_boot_w_svrg",
                    "modi_npmle_boot","sgd_boot_theta_svrg"), 
              envir=environment())


list_r <- parLapply(cl, 1:M,para_boot, kapprox=kapprox, sd=sd, bsiz=bsize)


stopCluster(cl);



for (i in 1:M){
  par(new=T)
  plot(list_r[[i]],col = alpha("pink",0.3),main="",xlim = ran)
}



#########################################BEGIN###################
### compare with true mixture cdf
#tt <- seq(-3,5, by = .0001)
tt <- seq(-3,20, by = .0001)
#lines(tt, pnorm(tt), lwd = 2, lty = 2, col="blue")
lines(tt, pgamma(tt,5,2), lwd = 2, lty = 2, col="blue",xlim = ran)
par(new=T)
#plot(kapprox,col="red",main="Standard Normal, sample size = 500, bootstrap size=100",xlim = ran)
plot(kapprox,col="red",main="Gamma(5,2), sample size = 500, bootstrap size=100",xlim = ran)
##############END#################################################


