getwd()

rm(list=ls())
if(T){
  library(GEOquery)
  eSet <- getGEO('GSE36278', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE36278_eSet.Rdata')
}
load('GSE36278_eSet.Rdata')
class(eSet)
length(eSet)
b = eSet[[1]]
raw_exprSet=exprs(b)
View(raw_exprSet)

phe=pData(b)
View(phe)
as.double(phe$`age at diagnosis (years):ch1`)
View(phe)
tem=phe[, c(2,37,42)]
names(tem)[2] <- "age"
names(tem)[3] <- "Group"
View(tem)
#tem <- as.data.frame(t(tem))
library(dplyr)
new_tem = filter(tem, age <= 23|is.na(age))
View(new_tem)
class(new_tem)
new_tem <- new_tem[-c(41:46),]
#new_tem <- as.data.frame(t(new_tem))
View(new_tem)
group_list = new_tem[,c(1,3)]
View(group_list)

group_list1 <- filter(group_list, Group=="WT"|Group=="K27M")
View(group_list1)

group_list2 <- filter(group_list, Group=="WT"|Group=="G34R")
View(group_list2)

raw_exprSet <- subset(raw_exprSet, select=c(rownames(group_list)))
View(raw_exprSet)

raw_exprSet1 <- subset(raw_exprSet, select=c(rownames(group_list1)))
View(raw_exprSet1)

raw_exprSet2 <- subset(raw_exprSet, select=c(rownames(group_list2)))
View(raw_exprSet2)



#group_list=paste0('group',group_list,'h')
save(raw_exprSet,group_list,raw_exprSet1,group_list1,raw_exprSet2,group_list2,
     file='GSE36278_raw_exprSet.Rdata')


