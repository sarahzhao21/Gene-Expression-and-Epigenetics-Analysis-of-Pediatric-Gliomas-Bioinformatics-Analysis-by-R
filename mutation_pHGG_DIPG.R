phgg <- read.csv("pHGGs_Mutated_Genes.csv", header=T)
View(phgg)
names(phgg)[5] <- "Frequency"
library(stringr)
phgg$Frequency = str_split(as.character(phgg$Frequency),'%',simplify = T)[,1]
phgg <- phgg[1:30,]
phgg$Frequency <- as.numeric(as.character(phgg$Frequency))
phgg$Gene <- as.character(as.factor(phgg$Gene))
View(phgg)
library(ggplot2)
phgg_mutant <- ggplot(data=phgg, aes(x=reorder(Gene,-Frequency), y=Frequency, colours=Gene)) +
  geom_bar(stat="identity", fill="steelblue",) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  geom_text(aes(label=Frequency), vjust=1.6, color="white", size=2) +
  labs(title="Mutated Genes in pHGG", x="Gene")

phgg_mutant

#------------------DIPG

dipg <- read.csv("DIPG_Mutated_Genes.csv", header=T)
View(dipg)
names(dipg)[5] <- "Frequency"
library(stringr)
dipg$Frequency = str_split(as.character(dipg$Frequency),'%',simplify = T)[,1]
dipg <- dipg[1:30,]
dipg$Frequency <- as.numeric(as.character(dipg$Frequency))
dipg$Gene <- as.character(as.factor(dipg$Gene))
View(dipg)
library(ggplot2)
dipg_mutant <- ggplot(data=dipg, aes(x=reorder(Gene,-Frequency), y=Frequency)) +
  geom_bar(stat="identity", fill="steelblue",) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  geom_text(aes(label=Frequency), vjust=1.6, color="white", size=2) +
  labs(title="Mutated Genes in DIPG", x="Gene")

dipg_mutant
