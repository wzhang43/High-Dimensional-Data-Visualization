
###### Run 10x10 CV on PCs obtained from naively combined
###### pathway 1 and 2 data

setwd("/Users/wzhang/Project 2")
Seeds = readRDS("mclustDA_BF/Stability_Sim/seedList.rds")
Seeds.pc = Seeds[1:10]

y.fil = readRDS("data/leaf_path12_nonunique.rds")
path.fil = readRDS("data/leaf_path12_label_nonunique.rds")

### Run PCA on combined data
out = prcomp(y.fil,center=T)
pcvar = (out$sdev)^2

## Percent of variation explained
sum(pcvar[1:2])/sum(pcvar) # 0.6649475
sum(pcvar[1:3])/sum(pcvar) # 0.7819958
sum(pcvar[1:4])/sum(pcvar) # 0.8560355

## PCs used for classification
pc.top2 = out$x[,1:2]
pc.top3 = out$x[,1:3]
pc.top4 = out$x[,1:4]


### 10-fold CV wrapper
cv.pca = function(data,class,seed){
  
  set.seed(seed)
  require(caret)
  cv = createFolds(1:nrow(data), k=10, list=T) 
  
  require(foreach)
  results <- foreach(fold = cv,.combine="c") %dopar% {
      require(mclust)
      data.train <- data[-fold]; class.train = class[-fold]
      data.test <- data[fold]; class.test = class[fold]
      fit = MclustDA(data.train, class=class.train, G=1:2, verbose=F)
      test.pred = predict(fit, newdata=data.test)
      error = 1-mean(test.pred$classification==class.test)
  }
  return(mean(results)) # output 
}


### CV on first 2 PCs
ptm = proc.time()
cv.error2 = sapply(Seeds.pc,cv.pca,data=pc.top2,class=path.fil)
runTime = proc.time() - ptm
print(runTime) # 1.552s
print(mean(cv.error2)) # first 2 PCs: error = 0.361803, GSI = 0.638197

### CV on first 3 PCs
ptm = proc.time()
cv.error3 = sapply(Seeds.pc,cv.pca,data=pc.top3,class=path.fil)
runTime = proc.time() - ptm
print(runTime) # 2.036s
print(mean(cv.error3)) # first 2 PCs: error = 0.3579141, GSI = 0.6420859

### CV on first 4 PCs
ptm = proc.time()
cv.error4 = sapply(Seeds.pc,cv.pca,data=pc.top4,class=path.fil)
runTime = proc.time() - ptm
print(runTime) # 2.357s
print(mean(cv.error4)) # first 2 PCs: error = 0.3605808, GSI = 0.6394192


