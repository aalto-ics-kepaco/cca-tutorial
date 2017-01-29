
rm(list=ls())

# Set the directory
setwd("/Users/vuurtio/Desktop/anno_iii/cca_survey/manuscript/examples")

#install.packages("PMA")
library(PMA)
source("http://bioconductor.org/biocLite.R")
biocLite("impute")

library(R.matlab)
data <- readMat("sparse_data.mat")

# standardize variables
x <- data$A
y <- data$B

out <- CCA(x,y,typex="standard",typez="standard",K=5, niter=80)
#print(out,verbose=FALSE)

writeMat("pmd_result.mat",wa = out$u, wb = out$v)





