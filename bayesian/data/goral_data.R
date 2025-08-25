library(readxl)
library(dplyr)
source('functions/sampling_functions/e2dist1.R')

buffer <- 0.564 # 계산 절차?

# DATA PREPROCESSING
goral <- read_excel('data/rawdata/G_ISRB.xlsx')
goral <- goral[,colSums(is.na(goral)) == 0]


goral <- goral %>% dplyr::select(-c('camera','userlabel')) %>% select(c('2023/01/01':'2023/03/29'))

#use more data : 174 days
#goral <- goral[,colSums(is.na(goral))==0]
n <- as.matrix(goral)
#n <- n[, 3:ncol(n)]
#storage.mode(n) <- "integer"



gtraps <- read.csv('data/rawdata/TDF_ISRB.csv')
coord.scale <- 1000 #km 단위로 변환

X <- as.matrix(gtraps[,2:3]) / coord.scale
X.scaled <- X[,1] - min(X[,1])
Xl.scaled <- min(X.scaled) - buffer
Xu.scaled <- max(X.scaled) + buffer

Y.scaled <- X[,2]-min(X[,2])
Yl.scaled <- min(Y.scaled) - buffer
Yu.scaled <- max(Y.scaled) + buffer

xlims.scaled <- c(Xl.scaled, Xu.scaled)
ylims.scaled <- c(Yl.scaled, Yu.scaled) 
areakm2.scaled <- (Xu.scaled - Xl.scaled) * (Yu.scaled - Yl.scaled)

X2 <- as.matrix(cbind(X.scaled, Y.scaled)) # scaled trap location matrix

goral <- list(n=n, X=X2, J=nrow(n), K=ncol(n), xlim=xlims.scaled, ylim=ylims.scaled, area=areakm2.scaled)

goral
