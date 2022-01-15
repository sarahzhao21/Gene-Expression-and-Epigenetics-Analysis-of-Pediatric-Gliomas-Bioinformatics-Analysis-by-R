getwd()

rm(list=ls())
if(T){
  library(GEOquery)
  eSet <- getGEO('GSE36245', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE36245_eSet.Rdata')
}
load('GSE36245_eSet.Rdata')
length(eSet)
b0 = eSet[[1]]
raw_exprSet0=exprs(b0)
View(raw_exprSet0)


phe0=pData(b0)
View(phe0)
group_list0 <- phe0[,c(2,37,39)]
View(group_list0)
names(group_list0)[3] <- "Group"
names(group_list0)[2] <- "age"

library(dplyr)
group_list0 = filter(group_list0, age <= 21)
group_list0 <- filter(group_list0, Group=="WT"|Group=="K27M")
group_list0 <- group_list0[,-2]
View(group_list0)

raw_exprSet0 <- subset(raw_exprSet0, select=c(rownames(group_list0)))
View(raw_exprSet0)


#group_list=paste0('group',group_list,'h')
save(raw_exprSet0,group_list0,
     file='GSE36245_raw_exprSet.Rdata')
#raw_exprSet is a matrix, group_list is a dataframe



rm(list=ls())
if(T){
  library(GEOquery)
  eSet <- getGEO('GSE34824', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE34824_eSet.Rdata')
}
load('GSE34824_eSet.Rdata')
length(eSet)
b1 = eSet[[1]]
raw_exprSet1=exprs(b1)
View(raw_exprSet1)
#raw_exprSet[1:4,1:4]
#raw_exprSet=raw_exprSet[,9:37]

phe1=pData(b1)
View(phe1)

group_list1 <- phe1[,c(2,38)]
names(group_list1)[2] <- "Group"
library(dplyr)
group_list1 <- filter(group_list1, Group=="wt"|Group=="K27M")
View(group_list1)
group_list1[c(1,3:11,13:20),2] <- "WT"
raw_exprSet1 <- subset(raw_exprSet1, select=c(rownames(group_list1)))
View(raw_exprSet1)

#-------------combine
load('GSE36245_raw_exprSet.Rdata')
group_list_all <- rbind(group_list0, group_list1)
View(group_list_all)

raw_exprSet_all <- merge(raw_exprSet0,raw_exprSet1,by="row.names",all.x=FALSE, all.y=FALSE)
View(raw_exprSet_all)
rownames(raw_exprSet_all) <- raw_exprSet_all[,1]
raw_exprSet_all <- raw_exprSet_all[, -1]
identical(colnames(raw_exprSet_all), rownames(group_list_all))

save(raw_exprSet_all, group_list_all,file="GSE36245_34824_raw_exprSet.Rdata")


