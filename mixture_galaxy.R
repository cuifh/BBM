library(rvalues)
library(scales)
library(lsei)
library(jmuOutlier)
library(parallel)
library(MASS)
library(mclust, quietly = T)

set.seed(1111)


data(galaxies)
samples = galaxies /1000

hist(samples,nclass = 20, xlab = "velocity (Mm/sec)", main = "Histogram of velocities")



GMM = Mclust(samples, model = "E")
k=GMM$G

w_0 = GMM$parameters$pro
theta_0 = GMM$parameters$mean
sigmasq_0 = rep(GMM$parameters$variance$sigmasq,k)



N <- length(samples) #orginal sample size
M=500
iter=10000
s=c(0,40)


# create a npmle list by npmle function
# will be modified later
SSS=cbind(samples,rep(1,N))
kapprox= npmle(SSS, family = gaussian, maxiter =1000)

Theta_0 = cbind(theta_0,sigmasq_0)

kapprox$support = Theta_0
kapprox$mix.prop = w_0



sgd_boot_x_svrg = function(support,w, sigmasq, bsize){
  sd <- sqrt(sigmasq)
  components <- sample(1:length(w), prob=w, size=bsize, replace=TRUE)
  x = rnorm(n=bsize,mean=support[components],sd=sd[components])
  return(x)
}

sgd_boot_w_svrg = function(support,w, sigmasq, x_new,eta){
  eta = eta
  sd = sqrt(sigmasq)
  eta = eta
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

sgd_boot_theta_svrg = function(support,w, sigmasq, x_new,eta){
  sd = sqrt(sigmasq)
  W = w*0
  for (i in 1:length(x_new)){
    new_w = w*0
    
    comp = dnorm(x_new[i], mean=support, sd=sd)
    
    grad1 = dnorm(x_new[i], mean=support, sd=sd)*(x_new-support)/(sd^2) #not change for codes_cts and discrete yet.
    grad2 = dnorm(x_new[i], mean=support, sd=sd)*(-0.5/sigmasq+0.5/(sigmasq)^2*(x_new-support)^2)
    
    l= length(w)
    grad2 = rep(sum(grad2*w),l)
    
    new_w = comp*(100*w)
    new_w = (100*w)/sum(new_w)
    
    W1 = new_w*grad1
    W2 = new_w*grad2
  }
  return(cbind(support+eta*W1,sigmasq+eta*W2))
}


bsize=1


para_boot <- function(M,kapprox,bsiz){
  support = kapprox$support
  w = kapprox$mix.prop
  eta=0.9
  
  
  for (i in 1:iter){
    eta = 1-1/(N+i)
    x = sgd_boot_x_svrg(support[,1],w,support[,2],bsize)
    ww=w
    w = sgd_boot_w_svrg(support[,1],w,support[,2],x,eta)
    support = sgd_boot_theta_svrg(support[,1],ww, support[,2], x,(1-eta)) #20000
  }
  kapprox$support = support
  kapprox$mix.prop = w
  return(kapprox)
  
}




clnum<-detectCores() 

cl <- makeCluster(getOption("cl.cores", clnum-1));

clusterExport(cl, c('para_boot',"iter","N","sgd_boot_x_svrg","bsize","sgd_boot_w_svrg",
                    "sgd_boot_theta_svrg"), 
              envir=environment())


list_r <- parLapply(cl, 1:M,para_boot, kapprox=kapprox, bsiz=bsize)


stopCluster(cl);

dnormix <- function(x, pi, mu, sigmasq){
  sigma = sqrt(sigmasq)
  K <- length(pi)
  
  ans = 0
  for (i in 1:K){
    ans = ans+ pi[i]*dnorm(x, mu[i], sigma[i])
  }
  return(ans)
}


curve(dnormix(x,kapprox$mix.prop, kapprox$support[,1], kapprox$support[,2]),
      from=s[1], to=s[2], xlab="velocity (Mm/sec)", ylab="density", 
      main = "BBM for GMM on Velocities of Galaxies, M=500",ylim = c(0,0.2))


for (i in 1:length(list_r)){
  curve(dnormix(x,list_r[[i]]$mix.prop, list_r[[i]]$support[,1], list_r[[i]]$support[,2]), 
        col = alpha("pink",0.3),add = T)
}

curve(dnormix(x,kapprox$mix.prop, kapprox$support[,1], kapprox$support[,2]), add = T)

rug(samples)

