
library(mclust)
library(MASS)
setwd("/Users/wzhang/Project 2/PeerJ")


##### Part 1: Compare RCV GSI variation between large and small samples ##### 

### Choose whether to repeat CV or not
computeGSI = function(dat,label,rep){

  if(rep==TRUE){
    out = numeric(10)
    for(i in 1:10){
      mod = MclustDA(dat,class=label,G=1:2,verbose=FALSE)
      cv = cvMclustDA(mod,verbose=F,nfold=10)
      out[i] = 1-cv$error
    }
    return(out)} else{
      mod = MclustDA(dat,class=label,G=1:2,verbose=FALSE)
      pred = predict(mod)
      return(mean(pred$classification == label))
    }
}



### Simulate data and calculate GSI with and without CV
rep_norep = function(k, m){

  ## k = feature pair index
  ## m = sample size multiplier
  
  ## Group labels and data
  path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
  y.fil = readRDS("Data/leaf_path12_nonunique.rds")
  temp_2g = combn(23,2)

  dat = y.fil[,temp_2g[,k]]

  out = MclustDA(dat,class=path.fil,G=1:2,verbose=FALSE)

  mod1 = out$models[[1]]
  mod2 = out$models[[2]]

  par1 = mod1$parameters
  par2 = mod2$parameters

  prop1 = par1$pro
  prop2 = par2$pro

  if(mod1$G == 2){
    # Mean vectors for group 1
    mean1_1 = par1$mean[,1]
    mean1_2 = par1$mean[,2]
    cov1_1 = par1$variance$sigma[,,1]
    cov1_2 = par1$variance$sigma[,,2]
  } else{
    mean1 = par1$mean
    cov1 = par1$variance$Sigma
  }

  if(mod2$G == 2){
    # Mean vectors for group 2
    mean2_1 = par2$mean[,1]
    mean2_2 = par2$mean[,2]
    cov2_1 = par2$variance$sigma[,,1]
    cov2_2 = par2$variance$sigma[,,2]
  } else{
    mean2 = par2$mean
    cov2 = par2$variance$Sigma
  }

  # Sample sizes
  n1 = (mod1$n)*m
  n2 = (mod2$n)*m


  ## Generate sample from mixture distribution:
  if(mod1$G == 2){
    U1 = runif(n1)
    dat.group1 = matrix(0,n1,2)
    for(i in 1:n1){
      if(U1[i]<prop1[1]){
        dat.group1[i,] = mvrnorm(1,mean1_1,cov1_1)
      } else{
        dat.group1[i,] = mvrnorm(1,mean1_2,cov1_2)
      }
    }} else{
      dat.group1 = mvrnorm(n1,mean1,cov1)
  }

  if(mod2$G == 2){
    U2 = runif(n2)
    dat.group2 = matrix(0,n2,2)
    for(i in 1:n2){
      if(U2[i]<prop2[1]){
        dat.group2[i,] = mvrnorm(1,mean2_1,cov2_1)
      } else{
        dat.group2[i,] = mvrnorm(1,mean2_2,cov2_2)
      }
    }} else{
      dat.group2 = mvrnorm(n2,mean2,cov2)
  }


  ## Combine data from 2 groups
  dat.sim = rbind(dat.group1, dat.group2)

  ## Generate group labels
  label.sim = rep(1:2, c(n1, n2))

  err.norep = computeGSI(dat.sim, label.sim, rep=FALSE) # no CV used
  err.repVec = computeGSI(dat.sim, label.sim, rep=TRUE) # CV without repetition
  err.rep = mean(err.repVec) # CV w/ repetition

  res = list()
  res$err.norep = err.norep
  res$err.repVec = err.repVec
  res$err.rep = err.rep

  saveRDS(res, paste("Sample Size/Results/Simulated Data/",k,"_",m,"x.rds",sep=""))
}


## GSI from simulated sample 10x
ptm = proc.time()
for(k in 1:choose(23,2)){
  print(paste("Cycle = ",k,sep=""))
  rep_norep(k, m=10)
}
runTime = proc.time() - ptm # 5240.142s

## GSI from simulated sample 1x
ptm = proc.time()
for(k in 1:choose(23,2)){
  print(paste("Cycle = ",k,sep=""))
  rep_norep(k, m=1)
}
runTime = proc.time() - ptm # 504.611s



## GSI from simulated sample 2x
ptm = proc.time()
for(k in 1:choose(23,2)){
  print(paste("Cycle = ",k,sep=""))
  rep_norep(k, m=2)
}
runTime = proc.time() - ptm # 722.054s

## GSI from simulated sample 5x
ptm = proc.time()
for(k in 1:choose(23,2)){
  print(paste("Cycle = ",k,sep=""))
  rep_norep(k, m=5)
}
runTime = proc.time() - ptm # 2040.070s

## GSI from simulated sample 8x
ptm = proc.time()
for(k in 1:choose(23,2)){
  print(paste("Cycle = ",k,sep=""))
  rep_norep(k, m=8)
}
runTime = proc.time() - ptm # 3841.284 



### Use original data to calculate GSI with and without RCV
rx = function(k){

  ## Group labels and data
  path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
  y.fil = readRDS("Data/leaf_path12_nonunique.rds")
  temp_2g = combn(23,2)

  dat = y.fil[,temp_2g[,k]]

  err.norep = computeGSI(dat, path.fil, rep=FALSE) # no CV used
  err.repVec = computeGSI(dat, path.fil, rep=TRUE) # CV without repetition
  err.rep = mean(err.repVec) # CV w/ repetition

  res = list()
  res$err.norep = err.norep
  res$err.repVec = err.repVec
  res$err.rep = err.rep

  saveRDS(res, paste("Sample Size/Results/Original Data/",k,".rds",sep=""))

}

### GSI from original sample
ptm = proc.time()
for(k in 1:choose(23,2)){
  print(paste("Cycle = ",k,sep=""))
  rx(k)
}
runTime = proc.time() - ptm # 546.143s



### Order RCV GSI values ### 
## Simulated sample 10x
gsiMat.10x = matrix(0, nrow=10, ncol=choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_10x.rds",sep="")
  res = readRDS(fname)
  gsiMat.10x[,k] = res$err.repVec
}

md.10x = apply(gsiMat.10x, 2, median)
ordered.10x = gsiMat.10x[,order(md.10x)]

## Simulated sample 8x
gsiMat.8x = matrix(0, nrow=10, ncol=choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_8x.rds",sep="")
  res = readRDS(fname)
  gsiMat.8x[,k] = res$err.repVec
}

md.8x = apply(gsiMat.8x, 2, median)
ordered.8x = gsiMat.8x[,order(md.8x)]

## Simulated sample 5x
gsiMat.5x = matrix(0, nrow=10, ncol=choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_5x.rds",sep="")
  res = readRDS(fname)
  gsiMat.5x[,k] = res$err.repVec
}

md.5x = apply(gsiMat.5x, 2, median)
ordered.5x = gsiMat.5x[,order(md.5x)]

## Simulated sample 2x
gsiMat.2x = matrix(0, nrow=10, ncol=choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_2x.rds",sep="")
  res = readRDS(fname)
  gsiMat.2x[,k] = res$err.repVec
}

md.2x = apply(gsiMat.2x, 2, median)
ordered.2x = gsiMat.2x[,order(md.2x)]

## Simulated sample 1x
gsiMat.1x = matrix(0, nrow=10, ncol=choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_1x.rds",sep="")
  res = readRDS(fname)
  gsiMat.1x[,k] = res$err.repVec
}

md.1x = apply(gsiMat.1x, 2, median)
ordered.1x = gsiMat.1x[,order(md.1x)]

## Original sample
gsiMat.original = matrix(0, nrow=10, ncol=choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Original Data/",k,".rds",sep="")
  res = readRDS(fname)
  gsiMat.original[,k] = res$err.repVec
}

md.original = apply(gsiMat.original, 2, median)
ordered.original = gsiMat.original[,order(md.original)]




## Set uniform lower/upper bounds for graphs
y.low = min(c(gsiMat.10x, gsiMat.8x, gsiMat.5x, gsiMat.2x, gsiMat.1x, gsiMat.original)) - 0.01
y.up = max(c(gsiMat.10x, gsiMat.8x, gsiMat.5x, gsiMat.2x, gsiMat.1x, gsiMat.original)) + 0.01

png("Sample Size/errbox_2.png",width=14,height=10,units="in",res=300)
par(mfrow=c(2,3))
boxplot(ordered.10x, ylim=c(y.low,y.up), main="Simulated Sample, 10x", cex=0.7,
  xlab="Rank", ylab="10 CV GSIs")
boxplot(ordered.8x, ylim=c(y.low,y.up), main="Simulated Sample, 8x", cex=0.7,
  xlab="Rank", ylab="10 CV GSIs")
boxplot(ordered.5x, ylim=c(y.low,y.up), main="Simulated Sample, 5x", cex=0.7,
  xlab="Rank", ylab="10 CV GSIs")
boxplot(ordered.2x, ylim=c(y.low,y.up), main="Simulated Sample, 2x", cex=0.7,
  xlab="Rank", ylab="10 CV GSIs")
boxplot(ordered.1x, ylim=c(y.low,y.up), main="Simulated Sample, 1x", cex=0.7,
  xlab="Rank", ylab="10 CV GSIs")
boxplot(ordered.original, ylim=c(y.low,y.up), main="Original Sample", cex=0.7,
  xlab="Rank", ylab="10 CV GSIs")
par(mfrow=c(1,1))
dev.off()






#####    Part 2: Compare RCV GSI and no-CV GSI   #####

## Simulated sample 10x
gsi.rcv.10x = gsi.nocv.10x = numeric(choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_10x.rds",sep="")
  res = readRDS(fname)
  gsi.rcv.10x[k] = res$err.rep
  gsi.nocv.10x[k] = res$err.norep
}

## Simulated sample 8x
gsi.rcv.8x = gsi.nocv.8x = numeric(choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_8x.rds",sep="")
  res = readRDS(fname)
  gsi.rcv.8x[k] = res$err.rep
  gsi.nocv.8x[k] = res$err.norep
}

## Simulated sample 5x
gsi.rcv.5x = gsi.nocv.5x = numeric(choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_5x.rds",sep="")
  res = readRDS(fname)
  gsi.rcv.5x[k] = res$err.rep
  gsi.nocv.5x[k] = res$err.norep
}

## Simulated sample 2x
gsi.rcv.2x = gsi.nocv.2x = numeric(choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_2x.rds",sep="")
  res = readRDS(fname)
  gsi.rcv.2x[k] = res$err.rep
  gsi.nocv.2x[k] = res$err.norep
}

## Simulated sample 1x
gsi.rcv.1x = gsi.nocv.1x = numeric(choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Simulated Data/",k,"_1x.rds",sep="")
  res = readRDS(fname)
  gsi.rcv.1x[k] = res$err.rep
  gsi.nocv.1x[k] = res$err.norep
}

## Original sample
gsi.rcv.org = gsi.nocv.org = numeric(choose(23,2))
for(k in 1:choose(23,2)){
  fname = paste("Sample Size/Results/Original Data/",k,".rds",sep="")
  res = readRDS(fname)
  gsi.rcv.org[k] = res$err.rep
  gsi.nocv.org[k] = res$err.norep
}

x.low = y.low = min(c(gsi.rcv.10x, gsi.rcv.org, gsi.nocv.10x, gsi.nocv.org,
  gsi.rcv.1x, gsi.nocv.1x, gsi.rcv.5x, gsi.nocv.5x, gsi.rcv.2x, gsi.nocv.2x,
  gsi.rcv.8x, gsi.nocv.8x)) - 0.01
x.up = y.up = max(c(gsi.rcv.10x, gsi.rcv.org, gsi.nocv.10x, gsi.nocv.org,
  gsi.rcv.1x, gsi.nocv.1x, gsi.rcv.5x, gsi.nocv.5x, gsi.rcv.2x, gsi.nocv.2x,
  gsi.rcv.8x, gsi.nocv.8x)) + 0.01

png("Sample Size/rcv_nocv.png",width=15,height=10,units="in",res=300)
par(mfrow=c(2,3))
plot(gsi.rcv.10x, gsi.nocv.10x, xlab="RCV GSI", ylab="No-CV GSI", 
  main="Simulated Sample, 10x", xlim=c(x.low,x.up), ylim=c(y.low,y.up))
abline(0, 1, lty="dashed")

plot(gsi.rcv.8x, gsi.nocv.8x, xlab="RCV GSI", ylab="No-CV GSI", 
  main="Simulated Sample, 8x", xlim=c(x.low,x.up), ylim=c(y.low,y.up))
abline(0, 1, lty="dashed")

plot(gsi.rcv.5x, gsi.nocv.5x, xlab="RCV GSI", ylab="No-CV GSI", 
  main="Simulated Sample, 5x", xlim=c(x.low,x.up), ylim=c(y.low,y.up))
abline(0, 1, lty="dashed")

plot(gsi.rcv.2x, gsi.nocv.2x, xlab="RCV GSI", ylab="No-CV GSI", 
  main="Simulated Sample, 2x", xlim=c(x.low,x.up), ylim=c(y.low,y.up))
abline(0, 1, lty="dashed")

plot(gsi.rcv.1x, gsi.nocv.1x, xlab="RCV GSI", ylab="No-CV GSI", 
  main="Simulated Sample, 1x", xlim=c(x.low,x.up), ylim=c(y.low,y.up))
abline(0, 1, lty="dashed")

plot(gsi.rcv.org, gsi.nocv.org, xlab="RCV GSI", ylab="No-CV GSI", 
  main="Original Sample", xlim=c(x.low,x.up), ylim=c(y.low,y.up))
abline(0, 1, lty="dashed")
par(mfrow=c(1,1))
dev.off()



