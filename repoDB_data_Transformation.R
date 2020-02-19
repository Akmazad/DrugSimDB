library(data.table)
library(tidyr)
library(purrr)
library(pbapply)
repoDB.dat <- fread("Data/repoDB_full.csv")
repoDB.dat <- select(repoDB.dat, c("drug_id", "ind_id", "status"))
repoDB.dat$status = ifelse(repoDB.dat$status == "Approved", 1, 0)
repoDB.dat = tidyr::unite(repoDB.dat, "Ind_Status", c("ind_id", "status"), sep="__")
newDat <- split(repoDB.dat$Ind_Status, repoDB.dat$drug_id)

findStatusofaPair <- function(x,y){
  d1 = data.frame(info=x) %>% 
        tidyr::separate(1,c("ind_id","status"), sep="__")
  # print(d1)
  d2 = data.frame(info=y) %>% 
    tidyr::separate(1,c("ind_id","status"), sep="__")
  # print(d2)
  com = dplyr::inner_join(d1,d2, by="ind_id")
  # print(com)
  if(nrow(com) > 0){
    tp_indx = which((com[,2] == 1) && (com[,3] == 1))
    if(length(tp_indx) != 0)
      return(1)
    else
      return(0)
  }
}

repoMat<- sapply(newDat[1:10], function(x) sapply(newDat[1:10], findStatusofaPair, x))
