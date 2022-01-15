getwd()
rm(list=ls()) 
load(file='GSE36245_34824_raw_exprSet.Rdata')

exprSet0=raw_exprSet_all
View(exprSet0)
library(hgu133plus2.db)

ids0=toTable(hgu133plus2SYMBOL)
length(unique(ids0$symbol))
tail(sort(table(ids0$symbol)))
table(sort(table(ids0$symbol)))
plot(table(sort(table(ids0$symbol))))

table(rownames(exprSet0) %in% ids0$probe_id)
dim(exprSet0)
exprSet0=exprSet0[rownames(exprSet0) %in% ids0$probe_id,]
dim(exprSet0)

ids0=ids0[match(rownames(exprSet0),ids0$probe_id),]
head(ids0)
View(ids0)
exprSet0[1:5,1:5]

changeid <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}

new_exprSet0 <- changeid(exprSet0,ids0)

View(new_exprSet0)

group_list0=group_list_all

save(new_exprSet0,group_list0,
     file='GSE36245_34824_new_exprSet.Rdata')

#load(file='GSE26576_new_exprSet.Rdata')







#---------part 2-------------

load(file='GSE50021_raw_exprSet.Rdata')

exprSet1=raw_exprSet
View(exprSet1)
group_list1=group_list
group_list1
?read.table
datatable <- read.csv("GPL13938-11302.csv")
View(datatable)

ids1 = datatable[,c(1,12)]
View(ids1)
#ids=toTable(hgu133plus2SYMBOL)
#length(unique(ids$symbol))
#tail(sort(table(ids$symbol)))
#table(sort(table(ids$symbol)))
#plot(table(sort(table(ids$symbol))))

table(rownames(exprSet1) %in% ids1$ID)
dim(exprSet1)
exprSet1=exprSet1[rownames(exprSet1) %in% ids1$ID,]
dim(exprSet1)


ids1=ids1[match(rownames(exprSet1),ids1$ID),]
head(ids1)
View(ids1)
exprSet1[1:5,1:5]

changeid <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$Symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$ID),2]
  return(exprSet)
}

new_exprSet1 <- changeid(exprSet1,ids1)

View(new_exprSet1)
save(new_exprSet1,group_list1,
     file='GSE50021_new_exprSet.Rdata')

load(file='GSE50021_new_exprSet.Rdata')

#---------------combine two datasets
group_list_all <- rbind(group_list0,group_list1)
View(group_list_all)
?merge

new_exprSet_all <- merge(new_exprSet0,new_exprSet1,by="row.names",all.x=FALSE, all.y=FALSE)
View(new_exprSet_all)
class(new_exprSet_all)
rownames(new_exprSet_all) <- new_exprSet_all[,1]
new_exprSet_all <- new_exprSet_all[, -1]
identical(colnames(new_exprSet_all), rownames(group_list_all))

save(new_exprSet_all, group_list_all,file="GSE26576_50021_new_exprSet.Rdata")

