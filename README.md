# Trevor_cig_smoke
Microarray data analysis for the smoking study

# SCRIPTS
## Header.R
import for all scripts using source

## Blood.R
Analysis of the blood data

## Nose.R
Analysis of the nose data


# FUNCTIONS

## f_lGetPCAClusterCount
### Desc:
Takes first two components of the PCA and counts possible clusters. The function does this by 
binning the vector of data into X bins, assigning a class label to each bin, counting how many
observations in each bin, and total number of bins with at least one observations. This is calculated
for both the components of the pca matrix, and the max number of bins with at least one observation, in
first or second dimension is reported.
### Args:
pr.out = principal component object returned by prcomp function
### Rets:
returns list with 2 elements: 
1 - cluster.count = possible number of clusters in the data  
2 - cluster.label = data.frame with cluster labels

