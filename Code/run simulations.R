
### Load required packages and source functions
library(caret)
library(mclust)
library(foreach)
library(doParallel)
library(R.utils)
source("Code/sim_functions.R")

### Create clusters for parallel computing
registerDoParallel(makeCluster(11)) 

### Import data
setwd("/home/stats/zhangwa")
y.fil = readRDS("Data/leaf_path12_nonunique.rds")
path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
Seeds = readRDS("Data/seedList.rds")
comb = combn(23,4)
I = choose(23,4)

### 
foreach(i=1:50,.combine="c") %dopar% {
      require(mclust)
      s = Seeds[i]
      sim.c2(s,10)
}


foreach(i=1:50,.combine="c") %dopar% {
      require(mclust)
      s = Seeds[i]
      sim.c3(s,10)
}


foreach(i=1:I) %dopar% {
      require(mclust)
      require(R.utils)
      sim.c4(index=i, data=y.fil, class=path.fil, fold=10,
        comb=comb, seedlist=Seeds)
}

