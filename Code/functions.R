

### Plot MclustDA results and highlight user-specified outliers (for dim >2)
plotellipse.outlier <- function(outlier, data, dimens, nclass, symbols, colors, models, scale, ...) {
  m <- lapply(models, function(m) {
        m$parameters$mean <- array(m$parameters$mean[dimens, ], c(2, m$G))
        m$parameters$variance$sigma <- array(m$parameters$variance$sigma[dimens, 
            dimens, ], c(2, 2, m$G))
        m
  })
  xup = yup = ceiling(max(data)) # set height and width same for all scatterplots
  xlow = ylow = floor(min(data))
  outlier.name = row.names(data)[outlier]
  plot(data[, dimens], type = "n", xlim=c(xlow,xup), ylim=c(ylow,yup),...)
  points(data[outlier,dimens], pch=1, col="red", cex=2.3)
  text(data[outlier,dimens], labels=outlier.name, cex=1.2, pos=3, col="red", offset=0.7)
  b = log(2)/scale
  lines(x=c(-b,b),y=c(b,b),lty="dashed")
  lines(x=c(-b,-b),y=c(b,-b),lty="dashed")
  lines(x=c(-b,b),y=c(-b,-b),lty="dashed")
  lines(x=c(b,b),y=c(-b,b),lty="dashed")
  abline(v=0,lty="dashed")
  abline(h=0,lty="dashed")
  for (l in 1:nclass) {
    I <- m[[l]]$observations
    group.col = colors[l]
    require(scales)
    pcol = alpha(group.col, 0.9) # make points more transparent
    gcol = alpha(group.col, 0.9) # make ellipses more transparent
    meancol = alpha("white",0)
    points(data[I, dimens[1]], data[I, dimens[2]], 
        pch = symbols[l], col = pcol, cex=1.3, lwd=2)
    
    for (k in 1:(m[[l]]$G)) {
      mvn2plot(mu = m[[l]]$parameters$mean[, k], 
        sigma = m[[l]]$parameters$variance$sigma[, , k], k = 15,
        col = c(gcol, "grey60", meancol), lwd=c(1.2, 1.2))
    }
  }
}



### Plot MclustDA result and highlight outliers (for dim=2)
### Note: requires function "plotellipse.outlier()"
plot.outlier = function(data, mod, outlier, scale){
  colors = c("magenta","darkcyan")
  symbols = c(3,2)
  dataNames = colnames(data)
  mod = mod$models
  plotellipse.outlier(outlier, data, 1:2, 2, symbols, colors, models=mod, scale=scale)
  legend("bottomright",legend=c("ET","JA"),col=c("magenta", "darkcyan"),
          pch=c(3,2),cex=1.2)
}




### Generate scatterplot matrices showing MclustDA results (for dim>2)
splom.mclustda = function(data, out, scale=1, width, size){

  gap = 0.2
  p = ncol(data)
  colors = c("magenta","darkcyan")
  symbols = c(3,2)
  dataNames = colnames(data)
  mod = out$models
  
  low = floor(min(data))
  up = ceiling(max(data))
  ticks = (low:up)[which((low:up)%%2 == 0)]

  par(mfrow = c(p, p), mar = rep(c(gap, gap/2), each = 2), 
                oma = c(4, 4, 4, 4))
  for (i in seq(p)) {
    for (j in seq(p)) {
      if (i == j) {
        plot(0, 0, type = "n", xlab = "", ylab = "", 
            axes = FALSE, xlim=c(low,up), ylim=c(low,up))
        text((low+up)/2, (low+up)/2, dataNames[i], cex = 1.5, adj = 0.5)
        box()
      }
      else {
        plotellipse(data, c(j, i), 2, symbols, 
                  colors, width=width, size=size, models=mod, scale=scale, xaxt = "n", yaxt = "n")
        legend("bottomright",legend=c("ET","JA"),col=c("magenta", "darkcyan"),
          pch=c(3,2),cex=0.8)
      }
      if (i == 1 && (!(j%%2))) 
        axis(3, at = ticks)
      if (i == p && (j%%2)) 
        axis(1, at = ticks)
      if (j == 1 && (!(i%%2))) 
        axis(2, at = ticks)
      if (j == p && (i%%2)) 
        axis(4, at = ticks)
    }
  }
}



plotellipse <- function(data, dimens, nclass, symbols, colors, models, scale, width, size, ...) {
  m <- lapply(models, function(m) {
        m$parameters$mean <- array(m$parameters$mean[dimens, ], c(2, m$G))
        m$parameters$variance$sigma <- array(m$parameters$variance$sigma[dimens, 
            dimens, ], c(2, 2, m$G))
        m
  })
  xup = yup = ceiling(max(data)) # set height and width same for all scatterplots
  xlow = ylow = floor(min(data))
  plot(data[, dimens], type = "n", xlim=c(xlow,xup), ylim=c(ylow,yup),...)
  b = log(2)/scale
  lines(x=c(-b,b),y=c(b,b),lty="dashed")
  lines(x=c(-b,-b),y=c(b,-b),lty="dashed")
  lines(x=c(-b,b),y=c(-b,-b),lty="dashed")
  lines(x=c(b,b),y=c(-b,b),lty="dashed")
  abline(v=0,lty="dashed")
  abline(h=0,lty="dashed")
  for (l in 1:nclass) {
    I <- m[[l]]$observations
    group.col = colors[l]
    require(scales)
    pcol = alpha(group.col, 0.9) # make points more transparent
    gcol = alpha(group.col, 0.9) # make ellipses more transparent
    meancol = alpha("white",0)
    points(data[I, dimens[1]], data[I, dimens[2]], 
        pch = symbols[l], col = pcol, cex=size, lwd=width)
    
    for (k in 1:(m[[l]]$G)) {
      mvn2plot(mu = m[[l]]$parameters$mean[, k], 
        sigma = m[[l]]$parameters$variance$sigma[, , k], k = 15,
        col = c(gcol, "grey60", meancol), lty=c(1,3), pch=1, lwd=c(1.2, 1.2))
    }
  }
}



### Generate scatterplot and show MclustDA results (for dim=2)
scatter.mclustda = function(data, out, scale){

  colors = c("magenta","darkcyan")
  symbols = c(3,2)
  dataNames = colnames(data)
  mod = out$models

  plotellipse(data, 1:2, 2, symbols, colors, models=mod, width=2, size=1.3, scale=scale)
  legend("bottomright",legend=c("ET","JA"),col=c("magenta", "darkcyan"),
          pch=c(3,2),cex=1.2)
}



### Plot scatterplots with specified feature subset "sset"
plot_minmax = function(sset){
  y.fil = readRDS("Data/leaf_path12_nonunique.rds")
  expLabel = colnames(y.fil)
  temp_2g = combn(23,2)
  combNames_2g = character(length=ncol(temp_2g))
  for(i in 1:length(combNames_2g)){
    combNames_2g[i] = paste("[",paste(expLabel[temp_2g[,i]],collapse=", "),"]",sep="")
  }
  index = match(sset, combNames_2g)
  n2 = combn(23,2)
  dat = y.fil[,n2[,index]]/log(2)
  out = MclustDA(dat,class=path.fil,G=1:2)
  fileName = paste("PeerJ/Graphs/log2 scale/2groups/",index,"_minmax.png",
    sep="")
  png(fileName,height=10,width=10,units="in",res=200)
  scatter.mclustda(dat, out, scale=log(2))
  dev.off()
}




### Generate scatterplot matrix with aspect ratio = 1 (dim>2)
plotellipse.asp1 = function(data, dimens, nclass, symbols, colors, models, width, size, ...) {
  m <- lapply(models, function(m) {
        m$parameters$mean <- array(m$parameters$mean[dimens, ], c(2, m$G))
        m$parameters$variance$sigma <- array(m$parameters$variance$sigma[dimens, 
            dimens, ], c(2, 2, m$G))
        m
  })
  xup = yup = ceiling(max(data)) # set height and width same for all scatterplots
  plot(data[, dimens], type = "n", asp=1,...)
  for (l in 1:nclass) {
    I <- m[[l]]$observations
    group.col = colors[l]
    require(scales)
    pcol = alpha(group.col, 0.9) # make points more transparent
    gcol = alpha(group.col, 0.9) # make ellipses more transparent
    meancol = alpha("white",0)
    points(data[I, dimens[1]], data[I, dimens[2]], 
        pch = symbols[l], col = pcol, cex=size, lwd=width)
    
    for (k in 1:(m[[l]]$G)) {
      mvn2plot(mu = m[[l]]$parameters$mean[, k], 
        sigma = m[[l]]$parameters$variance$sigma[, , k], k = 15,
        col = c(gcol, "grey60", meancol), lty=c(1,3), pch=1, lwd=c(1.2, 1.2))
    }
  }
}



### Generate scatterplot matrix with aspect ratio = 1 (dim=2)
scatter.asp1 = function(data, out){

  colors = c("magenta","darkcyan")
  symbols = c(3,2)
  dataNames = colnames(data)
  mod = out$models

  plotellipse.asp1(data, 1:2, 2, symbols, colors, models=mod, width=2, size=1.3)
  legend("bottomright",legend=c("ET","JA"),col=c("magenta", "darkcyan"),
          pch=c(3,2),cex=1.2)
}



### Generate scatterplot with LDA result and MclustDA result
plot.lda = function(sset, model="MclustDA", filename){
  setwd("/Users/wzhang/Project 2/PeerJ")
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
  if(model=="LDA"){
    out = MclustDA(dat,class=path.fil,modelType="EDDA")} else{
      out = MclustDA(dat,class=path.fil,G=1:2)}
  png(filename,height=7,width=11,units="in",res=200)
  scatter.asp1(dat, out)
  dev.off()
}


