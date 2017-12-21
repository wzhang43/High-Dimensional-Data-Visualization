


#####-------------------- Preamble --------------------------#####

### Load required packages and source graphing function
library(mclust)
source("Code/functions.R")

### Import data and group labels
setwd("/Users/wzhang/Project 2")
y.fil = readRDS("Data/leaf_path12_nonunique.rds")
path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")







####### Plot mclustDA fit on PCs #######

### Run PCA on combined data
out = prcomp(y.fil/log(2),center=T)
pcvar = (out$sdev)^2

## Percent of variation explained
sum(pcvar[1:2])/sum(pcvar) # 0.6649475
sum(pcvar[1:3])/sum(pcvar) # 0.7819958
sum(pcvar[1:4])/sum(pcvar) # 0.8560355

## PCs used for classification
pc.top2 = out$x[,1:2]
pc.top3 = out$x[,1:3]
pc.top4 = out$x[,1:4]

dat = pc.top2
out = MclustDA(dat,class=path.fil,G=1:2)
png("PeerJ/Graphs/pc2.png",height=10,width=10,units="in",res=200)
scatter.mclustda(dat, out, scale=log(2))
dev.off()




##################################################################
#####-------------------- Figure 1 --------------------------#####
##################################################################

### Plot scatterplot with max and min separation
plot_minmax("[4-7, 4-8]") # separation = 0.7082
plot_minmax("[2-1, 4-5]") # separation = 0.5355

##################################################################
#####------------------- End of Figure 1 --------------------#####
##################################################################









##################################################################
#####-------------------- Figure 2 --------------------------#####
##################################################################
  
### Import data, convert to log2 scale
temp = readRDS("Data/leaf_path12_nonunique.rds")
dat = temp/log(2)
  
### Manually input contrast IDs
trtID = c("1_1","2_1","2_2","2_3","3_1","3_2","4_1","4_2","4_3","4_4","4_5",
    "4_6","4_7","4_8","5_1","5_2","5_3","5_4","5_5","5_6","5_7","5_8","5_9",
    "6_1","6_2","7_1","7_2","7_3","8_1","8_2","9_1","10_1","10_2","10_3","10_4",
    "11_1","11_2","11_3","12_1","13_1","14_1","15_1","16_1","16_2","16_3","16_4",
    "17_1","17_2","17_3","17_4","17_5","17_6","17_7","17_8","18_1","18_2","19_1",
    "19_2","19_3","20_1","21_1","22_1","23_1","24_1","24_2","24_3","24_4","24_5",
    "24_6","24_7")

### Set lower and upper bounds for graph
n_leaf = ncol(dat)
yup = max(dat)+0.2
ylow = min(dat)-0.2
 
### Partition data by pathway group 
dat.eth = dat[which(path.fil==1),]
dat.jas = dat[which(path.fil==2),]

### Graph parameters
col1 = "grey"; col2 = "orange"
color = c(col1,rep(col2,2),col1,col2,col1,rep(col2,7),col1,rep(col2,8))
ltype = c(3,rep(1,2),3,1,3,rep(1,7),3,rep(1,8))
s = seq(n_leaf-1)

### Make lineplots
png("lineplot.png",width=17,height=7,units="in",res=300)
par(mfrow=c(1,2))
plot(1:n_leaf, dat.eth[1,], type="p",pch=2,cex=0.5,ylim=c(ylow,yup),ylab="log2 fold change",
  xlab="contrast ID",xaxt="n",lwd=0.7,cex.lab=1.2,cex.axis=1.2,cex.main=1.2,main="Ethylene Pathway Genes")
segments(s, dat.eth[1,s], s+1, dat.eth[1,s+1], col=color[s], lty=ltype[s])
for(k in 2:nrow(dat.eth)){
  lines(1:n_leaf, dat.jas[k,], type="p",pch=2,cex=0.5)
  segments(s, dat.eth[k,s], s+1, dat.eth[k,s+1], col=color[s], lty=ltype[s])
}
axis(1,at=1:n_leaf,labels=trtID[1:n_leaf])
  
plot(1:n_leaf, dat.jas[1,], type="p",pch=2,cex=0.5,ylim=c(ylow,yup),ylab="log2 fold change",
  xlab="contrast ID",xaxt="n",lwd=0.7,cex.lab=1.2,cex.axis=1.2,cex.main=1.2,main="Jasmonate Pathway Genes")
segments(s, dat.jas[1,s], s+1, dat.jas[1,s+1], col=color[s], lty=ltype[s])
for(k in 2:nrow(dat.jas)){
  lines(1:n_leaf, dat.jas[k,], type="p",pch=2,cex=0.5)
  segments(s, dat.jas[k,s], s+1, dat.jas[k,s+1], col=color[s], lty=ltype[s])
}
axis(1,at=1:n_leaf,labels=trtID[1:n_leaf])
par(mfrow=c(1,1))
dev.off()

##################################################################
#####------------------- End of Figure 2 --------------------#####
##################################################################







##################################################################
#####-------------------- Figure 3 --------------------------#####
##################################################################

plot.lda("[3-1, 4-2]", filename="Figure_3a.png")
plot.lda("[3-1, 4-2]", "LDA", "Figure_3b.png")


##################################################################
#####------------------- End of Figure 3 --------------------#####
##################################################################







##################################################################
#####----------------- Figure 4 and 5 -----------------------#####
##################################################################

### Figure 4 ###
n2 = combn(23,2)
data.199 = y.fil[,n2[,199]]
ja.199 = data.199[which(path.fil==2),]
et.199 = data.199[which(path.fil==1),]

## Which genes are outlying?
outlier.ja = names(sort(ja.199[,1],decreasing=T)[1:3])
outlier.et = names(c(sort(et.199[,2],decreasing=F)[1],sort(et.199[,2],decreasing=T)[1]))
outlier.name = union(outlier.ja, outlier.et)
outlier.199 = match(outlier.name,row.names(y.fil))

png("PeerJ/Graphs/log2 scale/example1_outliers.png",
  height=10,width=10,units="in",res=180)
mod.199 = MclustDA(data.199/log(2),class=path.fil,G=1:2)
plot.outlier(data.199,mod.199,outlier.199, scale=log(2))
dev.off()


### Figure 5 ###
n2 = combn(23,2)
data.93 = y.fil[,n2[,93]]
ja.93 = data.93[which(path.fil==2),]
et.93 = data.93[which(path.fil==1),]

## Which genes are outlying?
outlier.name = c(names(sort(ja.93[,1],decreasing=T)[1:2]), 
  names(sort(ja.93[,1])[1]), names(sort(et.93[,2])[1]))
outlier.93 = match(outlier.name,row.names(y.fil))

png("PeerJ/Graphs/log2 scale/example2_outliers.png",
  height=10,width=10,units="in",res=180)
mod.93 = MclustDA(data.93/log(2),class=path.fil,G=1:2)
plot.outlier(data.93,mod.93,outlier.93, scale=log(2))
dev.off()

##################################################################
#####------------- End of Figure 4 and 5 --------------------#####
##################################################################






##################################################################
#####----------------- Appendix Figures ---------------------#####
##################################################################

#### Graph scatterplot for top ranked features ####
### 4-subset scatterplots ###
n4 = combn(23,4)
top5_4g = readRDS("Results/4g_RCV_top5.rds")
for(i in 1:5){
  topInd = top5_4g[[i]]$index
  for(j in 1:5){
    dat = y.fil[,n4[,topInd[j]]]/log(2)
    out = MclustDA(dat,class=path.fil,G=1:2)
    fileName = paste("PeerJ/Graphs/4groups/Run ",i,"/4g_run",
      i,"_top",j,".png",sep="")
    png(fileName,height=10,width=10,units="in",res=200)
    splom.mclustda(dat, out, scale=log(2), width=1, size=0.9)
    dev.off()
  }
}


### 3-subset scatterplots ###
n3 = combn(23,3)
top5_3g = readRDS("Results/3g_RCV_top5.rds")
for(i in 1:5){
  topInd = top5_3g[[i]]$index
  for(j in 1:5){
    dat = y.fil[,n3[,topInd[j]]]/log(2)
    out = MclustDA(dat,class=path.fil,G=1:2)
    fileName = paste("PeerJ/Graphs/3groups/Run ",i,"/3g_run",
      i,"_top",j,".png",sep="")
    png(fileName,height=10,width=10,units="in",res=200)
    splom.mclustda(dat, out, scale=log(2), width=1, size=0.9)
    dev.off()
  }
}

### 3-subset scatterplots ###
n2 = combn(23,2)
top5_2g = readRDS("Results/2g_RCV_top5.rds")
for(i in 1:5){
  topInd = top5_2g[[i]]$index
  for(j in 1:5){
    dat = y.fil[,n2[,topInd[j]]]/log(2)
    out = MclustDA(dat,class=path.fil,G=1:2)
    fileName = paste("PeerJ/Graphs/2groups/Run ",i,"/2g_run",
      i,"_top",j,".png",sep="")
    png(fileName,height=10,width=10,units="in",res=200)
    col.ellipses.2g(dat, out, scale=log(2))
    dev.off()
  }
}

##################################################################
#####-------------- End of Appendix Figures -----------------#####
##################################################################





