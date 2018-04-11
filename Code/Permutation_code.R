
library(mclust)
setwd("/Users/wzhang/Project 2/PeerJ")
  
### 10 repetitions of 10-fold CV
computeGSI = function(dat,label){

  err = numeric(10)
  for(i in 1:10){
    out = MclustDA(dat,class=label,G=1:2,verbose=FALSE)
    #set.seed(seeds[i])
    cv = cvMclustDA(out,verbose=F,nfold=10)
    err[i] = cv$error
  }

  out = 1 - mean(err)
  return(out)
}


## Find the subset of data corresponding to query feature pair
findData = function(fpair){
  y.fil = readRDS("Data/leaf_path12_nonunique.rds")
  expLabel = colnames(y.fil)
  temp_2g = combn(23,2)
  combNames_2g = character(length=ncol(temp_2g))
  for(i in 1:length(combNames_2g)){
    combNames_2g[i] = paste("[",paste(expLabel[temp_2g[,i]],collapse=", "),"]",sep="")
  }
  index = match(fpair, combNames_2g)
  n2 = combn(23,2)
  dat = y.fil[,n2[,index]]/log(2)
  return(dat)
}


## Group labels and data
path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
dat = findData("[4-7, 4-8]")

## Record GSI value
gsi = numeric(2000)

ptm = proc.time()
for(i in 1:2000){
  label = sample(path.fil, length(path.fil)) # permute group labels
  gsi[i] = computeGSI(dat,label) # compute GSI
  print(paste("Cycle = ",i,sep="")) # keep track of iteration
}
runTime = proc.time() - ptm # 3568.586s 

saveRDS(gsi, "Permutation/Results/gsi_empirical.rds")

## Visualize GSI distribution 
png("Permutation/gsi_dist.png",width=10,height=8,units="in",res=300)
hist(gsi, xlab="GSI", main="GSI distribution")
dev.off()

## Observed GSI value
gsi.obs = computeGSI(dat,path.fil)
saveRDS(gsi.obs, "Permutation/Results/gsi_observed.rds")

## Calculate empirical p-value
p-val = mean(gsi>=gsi.obs); p-val 
# p-value = 0.005


