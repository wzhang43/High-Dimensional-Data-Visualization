### Choose 2 groups ###
sim.c2 = function(seed, fold){
  
  ## Choose 2 features to perform mclustDA
  mclustDA.choose2 = function(i,j,fold,seed){
    # Read data: leaf group pathways 1 and 2
    y.fil = readRDS("Data/leaf_path12_nonunique.rds")
    path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
    N = length(path.fil)
    group.names = colnames(y.fil)[c(i,j)]
    dat = y.fil[,c(i,j)] # subset data
    set.seed(seed)
    out = MclustDA(dat,class=path.fil,G=1:2,verbose=FALSE) # run mclustDA
    cv = cvMclustDA(out,verbose=F,nfold=fold)
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
  res.choose2[,3] = mapply(mclustDA.choose2,cn[1,],cn[2,],
    fold=fold,seed=seed)
  runTime = proc.time() - ptm

  res = list()
  res$time = runTime[3];
  res$seed = seed
  res$fold = fold
  res$error = res.choose2[!is.nan(res.choose2[,3]),]
  fileName = paste("Results/2groups/2g_s",seed,"_f",fold,".rds",sep="")
  saveRDS(res,fileName)
}



### Choose 3 groups ###
sim.c3 = function(seed, fold){
  
  ## Choose 3 features to perform mclustDA
  mclustDA.choose3 = function(i,j,k,fold,seed){
    print(paste("i=",i,", j=",j,", k=",k,sep=""))
    # Read data: leaf group pathways 1 and 2
    y.fil = readRDS("Data/leaf_path12_nonunique.rds")
    path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")
    N = length(path.fil)
    group.names = colnames(y.fil)[c(i,j,k)]
    dat = y.fil[,c(i,j,k)] # subset data
    set.seed(seed)
    out = MclustDA(dat,class=path.fil,G=1:2,verbose=FALSE) # run mclustDA
    cv = cvMclustDA(out,verbose=F,nfold=fold)
    cv.error = cv$error
    print("end")
    return(cv.error)
  }

  n = 23
  cn = combn(n,3) # all combinations of (23 choose 3)
  K3 = ncol(cn)
  res.choose3 = matrix(nrow=K3,ncol=4) # records fitting results
  colnames(res.choose3) = c("G1","G2","G3","CV_Error")
  res.choose3[,1:3] = t(cn)
  
  ptm = proc.time()
  res.choose3[,4] = mapply(mclustDA.choose3,cn[1,],cn[2,],cn[3,],
    fold=fold,seed=seed)
  runTime = proc.time() - ptm

  res = list()
  res$time = runTime[3];
  res$seed = seed
  res$fold = fold
  res$error = res.choose3[!is.nan(res.choose3[,4]),]
  fileName = paste("Results/3groups/3g_s",seed,"_f",fold,".rds",sep="")
  saveRDS(res,fileName)
}



### Choose 4 groups ###
sim.c4.new = function(index,data,class,fold,comb,seedlist){
  
  mclustDA.choose4 = function(seed,fold,index,comb,data,class){
    ind = comb[,index]
    print(paste("Seed=",seed,sep=""))
    group.names = colnames(data)[ind]
    dat = data[,ind]
    set.seed(seed)
    out = MclustDA(dat,class=class,G=1:2,verbose=FALSE) # run mclustDA
    
    cv.error = tryCatch(
      expr = {
            evalWithTimeout({cvMclustDA(out,verbose=F,nfold=fold)$error},timeout = 1)
      },
      TimeoutException = function(ex) return(NA)
    )

    return(cv.error)
  }

  ptm = proc.time()
  res = sapply(seedlist,mclustDA.choose4,fold=10,index=index,comb=comb,
    data=data,class=class)    
  runTime = proc.time() - ptm   

  fileName = paste("Results/4groups_new/comb",index,".rds",sep="")
  saveRDS(res,fileName)
}


