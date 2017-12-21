

##---------- Step 0: Fit NB regression model to read counts ------------##

rm(list=ls());
## Load all data
file.leaf = "Data/leaf.rds";
leaf = readRDS(file.leaf);

file.seedling = "Data/seedling.rds";
seedling = readRDS(file.seedling);

file.tissue = "Data/tissue.rds";
tissue = readRDS(file.tissue);

## Merge data into a single data.matrix
counts = as.matrix(cbind(leaf$count, seedling$count, tissue$count));

## Renumber the labs and find corresponding GEO (dataset ID)
lab = c(leaf$lab, seedling$lab + 5, tissue$lab + 14);
trt = c(leaf$trt, seedling$trt, tissue$trt);
GEO = c(leaf$GEO[leaf$lab], seedling$GEO[seedling$lab], tissue$GEO[tissue$lab]);
data.frame(lab, GEO, trt);

## Read reference genes
ref.genes = read.csv("Data/OverlappedGeneFromThreeSetsTop1000.csv", header=TRUE,
  stringsAsFactors = FALSE)[,1];

## Load NBPSeq and fit nb glm
library(NBPSeq);

norm.factors.0 = estimate.norm.factors(counts);
norm.factors = estimate.norm.factors(counts[ref.genes,], colSums(counts));
nb.data = prepare.nb.data(counts, norm.factors=norm.factors);
lib.sizes = colSums(counts)
eff.lib.sizes = norm.factors*lib.sizes
saveRDS(eff.lib.sizes, file="Data/effectivelibsize.rds")

fit.it = function(lab.id) {
  nb.data.1 = nb.data[, lab==lab.id];
  x = model.matrix(~factor(trt[lab == lab.id]));
  fit.nb.glm(nb.data.1, x); 
}

for (i in 1:24) {
  nb.glm = fit.it(i);
  file.out = sprintf("../Output/lab.%d.nb.glm.rds", i);
  saveRDS(nb.glm, file=file.out);
}






##----------------- Step 1: Aggregate gene expression data --------------##
## Load counts data
file.leaf = "Data/leaf.rds";
leaf = readRDS(file.leaf);
str(leaf);
genes = rownames(leaf$count);

## Load NB glm models fitted to 24 data sets
nb.glm.list = list();
for (i in 1:24) {
  file.rds = sprintf("Output/lab.%d.nb.glm.rds", i);
  nb.glm.list[[i]] = readRDS(file.rds);
}

## Results from one lab
library(NBPSeq);
m = nb.glm.list[[1]];
beta.mat = m$beta[,-1]
beta.num = ncol(m$beta)-1

## Combine betas from 24 experiments
for(i in 2:24){
  temp = nb.glm.list[[i]]
  beta.mat = cbind(beta.mat,temp$beta[,-1])
  beta.num = c(beta.num,(ncol(temp$beta)-1))
}
rownames(beta.mat) = genes;

## Combine mus from 24 experiments
mu.mat = m$mu
for(i in 2:24){
  temp = nb.glm.list[[i]]
  mu.mat = cbind(mu.mat,temp$mu)
}
rownames(mu.mat) = genes;

## Filter out genes with NA mu/beta values
beta.mat.nna = beta.mat[complete.cases(beta.mat),]
mu.mat.nna = mu.mat[complete.cases(mu.mat),]

## Filter: min(mu) > 1
rmin = apply(mu.mat.nna,1,min)
rkeep = rmin>1
beta.fil = beta.mat.nna[rkeep,]
trtID = c("1-1","2-1","2-2","2-3","3-1","3-2","4-1","4-2","4-3","4-4","4-5",
    "4-6","4-7","4-8","5-1","5-2","5-3","5-4","5-5","5-6","5-7","5-8","5-9",
    "6-1","6-2","7-1","7-2","7-3","8-1","8-2","9-1","10-1","10-2","10-3","10-4",
    "11-1","11-2","11-3","12-1","13-1","14-1","15-1","16-1","16-2","16-3","16-4",
    "17-1","17-2","17-3","17-4","17-5","17-6","17-7","17-8","18-1","18-2","19-1",
    "19-2","19-3","20-1","21-1","22-1","23-1","24-1","24-2","24-3","24-4","24-5",
    "24-6","24-7")
colnames(beta.fil) = trtID
saveRDS(beta.fil,"Data/beta_mean_filtered.rds")
beta.fil = readRDS("Data/beta_mean_filtered.rds")





##---------------- Step 2: Subset data by pathway groups ---------------###

subset.path = function(file){
  ## Pathway gene list manually downloaded from AmiGO2
  tmp = read.csv(file,header=F)
  a = grep("^AT",tmp$V1,value=T) # keep only gene names starting with "AT"
  genelist = a[nchar(a)==9]

  ## Which genes are included in the 24 experiments?
  Gnames = intersect(genelist,row.names(beta.fil))

  ## Subset expression data with those genes
  Gbeta = beta.fil[row.names(beta.fil) %in% Gnames,]
  return(Gbeta)
}

beta.eth = subset.path("Data/eth.csv")
beta.jas = subset.path("Data/jas.csv")

saveRDS(beta.eth,"Output/beta_path1.rds")
saveRDS(beta.jas,"Output/beta_path2.rds")






##--------------- Step 3: Combine data from the two pathways ------------##

x1 = readRDS("Output/beta_path1.rds");
x2 = readRDS("Output/beta_path2.rds");

x = rbind(x1,x2)
path = c(rep(1,nrow(x1)),rep(2,nrow(x2)))

## Filter out genes with small log fold changes
N = nrow(x)
rkeep = numeric(N)
for(i in 1:N){
  rkeep[i] = sum(abs(x[i,])>=log(2))!=0
}
x.fil = x[which(rkeep==1),]
path.fil = path[which(rkeep==1)]

## Export pathways 1 and 2 for leaf tissue group
k = length(path.fil)
x12.fil = x.fil[,1:23]
saveRDS(x12.fil,"Data/leaf_path12_nonunique.rds")
saveRDS(path.fil,"Data/leaf_path12_label_nonunique.rds")





