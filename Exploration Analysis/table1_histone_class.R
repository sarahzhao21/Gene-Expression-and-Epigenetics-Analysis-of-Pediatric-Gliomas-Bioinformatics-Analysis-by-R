getwd()

table1 <- read.csv("ts1.csv", header=T)
summary(table1)
str(table1)
View(table1)
library(ggplot2)
df1 <- table1[,c(4,7,10)]
View(df1)
library(tidyr)
df1 <- df1 %>% drop_na()
names(df1)[1] <- "Location"
names(df1)[2] <- "Age"
names(df1)[3] <- "Histone_Mutation"
View(df1)
a <- ggplot(data=df1, aes(x=Histone_Mutation, y=Age,
                             colour=Histone_Mutation))
a + geom_boxplot()
#a + geom_boxplot(size=1) + geom_point()
a + geom_jitter(size=0.5) + geom_boxplot(size=0.5, alpha=0.5)+
  facet_grid(Location~., scales="free") + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

df2 <- table1[,c(4,8,10)]
View(df2)
library(tidyr)
df2 <- df2 %>% drop_na()
names(df2)[1] <- "Location"
names(df2)[3] <- "Histone_Mutation"
View(df2)
b <- ggplot(data=df2, aes(x=Histone_Mutation, y=Survival,
                             colour=Histone_Mutation))
b + geom_boxplot()

b + geom_jitter(size=0.5) + geom_boxplot(size=0.5, alpha=0.5)+
  facet_grid(Location~., scales="free")+ 
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

#------------------------------------------------------------
df3 <- table1[,c(1,10)]
View(df3)
library(tidyr)
df3 <- df3 %>% drop_na()
rownames(df3) <- df3[,1]
names(df3)[2] <- "Group"

table2 <- read.csv("ts2_methyl.csv", header=T)
View(table2)
colnames(table2) <- table2[1,]
table2 <- table2[-1,]
rownames(table2) <- table2[,1]
table2 <- table2[,-1]

options(stringsAsFactors = F)
require(Biobase)
library("impute")
library(data.table)
#get rid of the N/A value
beta=as.matrix(table2)
beta=impute.knn(beta)
#prevent problem when do log, so add an 0.00001 to each value
betaData=beta$data
betaData=betaData+0.00001
table_methyl=betaData
View(table_methyl)
#merge_df <- merge(table_methyl, df3, by="row.names",all.x=FALSE, all.y=FALSE)
#View(merge_df)
df3=df3[match(rownames(table_methyl), rownames(df3)),]

table_methyl<- as.data.frame(t(table_methyl))
table_methyl <- data.matrix(table_methyl)
identical(colnames(table_methyl),rownames(df3))

library(dplyr)    
df3a <- filter(df3,Group=="Wild-type")
View(df3a)
df3b <- filter(df3,Group=="H3.3_K27M")
table(df3$Group)
df3c <- filter(df3,Group=="H3.1_K27M")
df3d <- filter(df3,Group=="H3.3_G34R")
df3e <- filter(df3,Group=="H3.3_G34V")

df3_new <- rbind.data.frame(df3a, df3b, df3c, df3d, df3e)
View(df3_new)
v <- c(rownames(df3_new))
v
table_methyl[,c(1:420)] <- table_methyl[,v]
View(table_methyl)
library(pheatmap)
pheatmap(table_methyl)
         


library(ChAMP)
myLoad=champ.filter(beta = table_methyl,pd = df3)
myLoad

save(myLoad, file = "step1_data1000_output.Rdata")



#---------------------------------------------------------------------------
library(survival)
library(survminer)
#library(ranger)
library(dplyr)
library(ggfortify)
install.packages("ggfortify")
install.packages("ranger")
surv_obj <- survfit(Surv(Survival, Alive) ~ 1, data=table1)
surv_obj
autoplot(surv_obj)

fit <- survfit(Surv(Survival, Alive) ~ Gender, data=table1)
ggsurvplot(fit, data=table1, pva1=TRUE)
autoplot(fit)
ggplot2::autoplot(fit)
