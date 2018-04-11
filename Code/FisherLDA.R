

#### Using Fisher's LDA as classifier to calculate Group Separation Index ####

### Choose 2 groups ###
sim.c2.fisher = function(seed){
  
  ## Choose 2 features to perform mclustDA
  mclustDA.choose2 = function(i,j,seed){
    # Read data: leaf group pathways 1 and 2
    y.fil = readRDS("Data/leaf_path12_nonunique.rds")
    path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
    N = length(path.fil)
    group.names = colnames(y.fil)[c(i,j)]
    dat = y.fil[,c(i,j)] # subset data
    set.seed(seed)
    out = MclustDA(dat,class=path.fil,verbose=FALSE,modelType="EDDA") # run mclustDA
    cv = cvMclustDA(out,verbose=F,nfold=10)
    cv.error = cv$error
    return(cv.error)
  }

  n = 23
  cn = combn(n,2) # all combinations of (23 choose 2)
  K2 = ncol(cn)
  res.choose2 = matrix(nrow=K2,ncol=3) # records fitting results
  colnames(res.choose2) = c("G1","G2","CV_Error")
  res.choose2[,1:2] = t(cn)
  
  ptm = proc.time()
  res.choose2[,3] = mapply(mclustDA.choose2,cn[1,],cn[2,],seed=seed)
  runTime = proc.time() - ptm

  res = list()
  res$time = runTime[3]; #print(runTime[3])
  res$seed = seed
  res$error = res.choose2[!is.nan(res.choose2[,3]),]
  fileName = paste("FisherLDA/fisher_s",seed,".rds",sep="")
  saveRDS(res,fileName)
}



library(mclust)

setwd("/Users/wzhang/Project 2/PeerJ")
Seeds = readRDS("Data/seedList.rds")

for(i in 1:50){
      s = Seeds[i]
      sim.c2.fisher(s)
      print(c(i,Seeds[i]))
}


### Tabulate and save simulation results
rank_tab = function(Seeds){
  errors_tab = matrix(nrow=length(Seeds),ncol=ncol(combn(23,2)))
  for(i in 1:length(Seeds)){
    seed = Seeds[i]
    fileName = paste("FisherLDA/fisher_s",seed,".rds",sep="")
    res = readRDS(fileName)
    errors_tab[i,] = res$error[,3]
  }
  fileName.1 = paste("FisherLDA/err.rds",sep="")
  fileName.2 = paste("FisherLDA/rank.rds",sep="")
  saveRDS(errors_tab, fileName.1)
  rank_tab = apply(errors_tab,1,rank)
  saveRDS(rank_tab, fileName.2)
}

# Save Fisher LDA results
rank_tab(Seeds)



y.fil = readRDS("Data/leaf_path12_nonunique.rds")
expLabel = colnames(y.fil)
temp_2g = combn(23,2)
combNames_2g = character(length=ncol(temp_2g))
for(i in 1:length(combNames_2g)){
  combNames_2g[i] = paste("[",paste(expLabel[temp_2g[,i]],collapse=", "),"]",sep="")
}


err.fisher = readRDS("FisherLDA/Results/Summary/err.rds")
run = list(5)
top5 = list(5)

for(i in 1:5){
  l = (i-1)*10+1
  u = i*10
  err = err.fisher[l:u,]
  err_mean = apply(err,2,mean)
  names(err_mean) = combNames_2g
  list.sort = sort(err_mean,decreasing=F)
  run[[i]] = list.sort
  gsi = 1-list.sort[1:5]
  index = match(names(gsi), combNames_2g)
  t5 = list(gsi=gsi, index=index)
  top5[[i]] = t5
}

saveRDS(run,"FisherLDA/Results/Summary/Fisher_fulllist.rds")
saveRDS(top5,"FisherLDA/Results/Summary/Fisher_top5.rds")


source("Code/functions.R")
library(scales)
library(mclust)

### Plot scatterplots with specified feature subset "sset"
plot_anypair = function(sset, fname){
  y.fil = readRDS("Data/leaf_path12_nonunique.rds")
  path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
  expLabel = colnames(y.fil)
  temp_2g = combn(23,2)
  combNames_2g = character(length=ncol(temp_2g))
  for(i in 1:length(combNames_2g)){
    combNames_2g[i] = paste("[",paste(expLabel[temp_2g[,i]],collapse=", "),"]",sep="")
  }
  index = match(sset, combNames_2g)
  n2 = combn(23,2)
  dat = y.fil[,n2[,index]]/log(2)
  out = MclustDA(dat,class=path.fil,modelType = "EDDA")

  colors = c("magenta","darkcyan")
  symbols = c(3,2)
  dataNames = colnames(dat)
  mod = out$models

  png(fname,height=10,width=10,units="in",res=200)
  plotellipse(dat, 1:2, 2, symbols, colors, models=mod, width=2, size=1.3, scale=log(2))
  legend("bottomright",legend=c("ET","JA"),col=c("magenta", "darkcyan"),
          pch=c(3,2),cex=1.2)
  dev.off()
}



plot_anypair("[3-1, 5-9]", "31_59.png")
plot_anypair("[3-1, 5-5]", "31_55.png")
plot_anypair("[4-4, 4-5]", "44_45.png")





