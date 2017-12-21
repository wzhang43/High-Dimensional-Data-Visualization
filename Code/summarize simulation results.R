

y.fil = readRDS("Data/leaf_path12_nonunique.rds")
Seeds = readRDS("Data/seedList.rds")

### Tabulate and save simulation results
rank_tab = function(Seeds,fold,k){

  errors_tab = matrix(nrow=length(Seeds),ncol=ncol(combn(23,k)))
  for(i in 1:length(Seeds)){
    seed = Seeds[i]
    fold = fold
    fileName = paste("Results/",k,"groups/",k,"g_s",seed,"_f",fold,".rds",sep="")
    res = readRDS(fileName)
    errors_tab[i,] = res$error[,(k+1)]
  }
  fileName.1 = paste("Results/Summary/",k,"g_",fold,"fold_err.rds",sep="")
  fileName.2 = paste("Results/Summary/",k,"g_",fold,"fold_rank.rds",sep="")
  saveRDS(errors_tab, fileName.1)
  rank_tab = apply(errors_tab,1,rank)
  saveRDS(rank_tab, fileName.2)
}

# Save 2-subset results
rank_tab(Seeds,10,2)

# Save 3-subset results
a1 = which(Seeds==3462)
a2 = which(Seeds==9879)
Seeds_3g = c(Seeds[-c(a1,a2)], 12, 31)
rank_tab(Seeds_3g,10,3)    
    
    




################ 2-subsets ####################

expLabel = colnames(y.fil)
temp_2g = combn(23,2)
combNames_2g = character(length=ncol(temp_2g))
for(i in 1:length(combNames_2g)){
  combNames_2g[i] = paste("[",paste(expLabel[temp_2g[,i]],collapse=", "),"]",sep="")
}


err2_f10 = readRDS("Results/Summary/2g_10fold_err.rds")
run = list(5)
top5 = list(5)

for(i in 1:5){
  l = (i-1)*10+1
  u = i*10
  err = err2_f10[l:u,]
  err_mean = apply(err,2,mean)
  names(err_mean) = combNames_2g
  list.sort = sort(err_mean,decreasing=F)
  run[[i]] = list.sort
  gsi = 1-list.sort[1:5]
  index = match(names(gsi), combNames_2g)
  t5 = list(gsi=gsi, index=index)
  top5[[i]] = t5
}

saveRDS(Seeds,"Results/Summary/Seeds_2g.rds")
saveRDS(run,"Results/Summary/2g_RCV_fulllist.rds")
saveRDS(top5,"Results/Summary/2g_RCV_top5.rds")





################ 3-subsets ####################

expLabel = colnames(y.fil)
temp_3g = combn(23,3)
combNames_3g = character(length=ncol(temp_3g))
for(i in 1:length(combNames_3g)){
  combNames_3g[i] = paste("[",paste(expLabel[temp_3g[,i]],collapse=", "),"]",sep="")
}


err3_f10 = readRDS("Results/Summary/3g_10fold_err.rds")
run = list(5)
top5 = list(5)

for(i in 1:5){
  l = (i-1)*10+1
  u = i*10
  err = err3_f10[l:u,]
  err_mean = apply(err,2,mean)
  names(err_mean) = combNames_3g
  list.sort = sort(err_mean,decreasing=F)
  run[[i]] = list.sort
  gsi = 1-list.sort[1:5]
  index = match(names(gsi), combNames_3g)
  t5 = list(gsi=gsi, index=index)
  top5[[i]] = t5
}

saveRDS(Seeds,"Results/Summary/Seeds_3g.rds")
saveRDS(run,"Results/Summary/3g_RCV_fulllist.rds")
saveRDS(top5,"Results/Summary/3g_RCV_top5.rds")





################ 4-subsets ####################

expLabel = colnames(y.fil)
temp_4g = combn(23,4)
combNames_4g = character(length=ncol(temp_4g))
for(i in 1:length(combNames_4g)){
  combNames_4g[i] = paste("[",paste(expLabel[temp_4g[,i]],collapse=", "),"]",sep="")
}


setwd("/Results/4groups/")
allFiles = list.files(path='.', full.names = TRUE)
temp = regmatches(allFiles, regexpr("[[:digit:]]+", allFiles))
fileNo = sort(as.numeric(temp))
resMat = matrix(0,nrow=length(fileNo),ncol=50)

for(i in 1:length(fileNo)){
  fileName = paste("comb",fileNo[i],".rds",sep="")
  resMat[i,] = readRDS(fileName)
}

row.names(resMat) = combNames_4g[fileNo]
res.comp = resMat[complete.cases(resMat),]

run = list(5)
top5 = list(5)
for(i in 1:5){
  l = (i-1)*10+1
  u = i*10
  err = resMat[,l:u]
  err_mean = apply(err,1,mean,na.rm=TRUE)
  list.sort = sort(err_mean,decreasing=F)
  run[[i]] = list.sort
  gsi = 1-list.sort[1:5]
  index = match(names(gsi), combNames_4g)
  t5 = list(gsi=gsi, index=index)
  top5.1[[i]] = t5
}

saveRDS(run,"Results/Summary/4g_RCV_fulllist.rds")
saveRDS(top5,"Results/Summary/4g_RCV_top5.rds")



