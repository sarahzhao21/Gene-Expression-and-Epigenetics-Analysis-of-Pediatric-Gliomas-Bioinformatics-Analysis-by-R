getwd()
rm(list=ls())
load(file='GBM_combine_DEG.Rdata')
source('functions.R')
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

View(nrDEG)
gene=rownames(nrDEG)
class(gene)
length(gene)
### get the universal genes and sDEG 

gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

###KEGG anaysis
kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)[,1:6]

#GSEA analysis
data(geneList, package="DOSE")
geneList
boxplot(geneList)
boxplot(nrDEG$logFC)
nrDEG$logFC
###we found that geneList is just nrDEG$logFC
geneList=nrDEG$logFC
names(geneList)=rownames(nrDEG)
geneList=sort(geneList, decreasing=T) #order the geneList in a proper order
#then we need to change the geneList ID into ENTREZID
genelist.df <- bitr(names(geneList), 
                    fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)
head(genelist.df)
tmp=data.frame(SYMBOL=names(geneList),
               logfc=as.numeric(geneList))
View(tmp)
tmp=merge(tmp,genelist.df, by='SYMBOL')
View(tmp)
#so we have the dataframe with SYMBOL, ENTREZID, and ENSEMBL
#now we creat a new geneList
geneList=tmp$logfc
View(geneList)
names(geneList)=genelist.df$ENTREZID
geneList=sort(geneList, decreasing=T) 
View(geneList)
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.1,
               verbose      = FALSE)
head(kk2)[,1:6]

gseaplot(kk2, geneSetID = "hsa04145")

save(kk, kk2 , file="GBM_Combine_KEGG.Rdata")

