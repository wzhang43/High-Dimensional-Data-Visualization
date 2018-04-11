
library(tspair)
library(mclust)
setwd("/Users/wzhang/Project 2/PeerJ")
source("Code/functions.R")

##### Use tspair to find top ranked feature pair #####

## Import data
y.fil = readRDS("Data/leaf_path12_nonunique.rds")
path.fil = readRDS("Data/leaf_path12_label_nonunique.rds")

## Transpose data for tsp method
d = t(y.fil)

## Use tspair to find top feature pair
tp = tspcalc(d, path.fil)
tp
# tsp object with: 1 TSPs
# Pair:		TSP Score	 Tie-Breaker		Indices
# TSP 1 : 	 0.3 		 NA 				 12 14 

## Plot the selected feature pairs
colnames(y.fil)[c(12,14)] # "4-6" "4-8"
plot_anypair("[4-6, 4-8]", "TS Pairs/46_48.png")
