
### Hierarchical clustering

#Option 1: All NSAF values and detected proteins
#Option 2: Filtered NSAF values
#Optional: remove day 0 and day 15

#####################################  OPTION 1: Averaged replicate values  #################################################################################

NSAF.avg <- read.csv("Documents/robertslab/labnotebook/data/NSAF/silo3and9_nozerovals_AVGs.csv")
rownames(NSAF.avg) <- paste("D",NSAF.avg$day, "T", NSAF.avg$temp, sep = "")
NSAF.avg <- t(NSAF.avg)

#Remove contaminant proteins
NSAF.avg <- subset(NSAF.avg, grepl(paste('CHOYP', collapse="|"), rownames(NSAF.avg)))
which(grepl('ALBU_BOVIN', rownames(NSAF.avg))) #ensure contaminants are gone

#dataframe
NSAF.trans <- data.frame(NSAF.avg)

##############################################################################################################################

#####################################  OPTION 2: Averaged and filtered values  #################################################################################

#NSAF.fil <- read.csv("Documents/robertslab/labnotebook/data/NSAF/silo3and9_nozerovals_noincnstprot.csv")
#rownames(NSAF.fil) <- paste("D",NSAF.fil$day, "T", NSAF.fil$temp, sep = "")
#NSAF.fil <- t(NSAF.fil)

#Remove contaminant proteins
#NSAF.fil <- subset(NSAF.fil, grepl(paste('CHOYP', collapse="|"), rownames(NSAF.fil)))
#which(grepl('ALBU_BOVIN', rownames(NSAF.fil))) #ensure contaminants are gone

#Because the rest of the code was set for NSAF.trans not NSAF.avg
#NSAF.trans <- data.frame(NSAF.fil)

#####################################  For genes as names  #############mn#######################################################

#set gene names as row names
#annotations <- read.csv("Documents/robertslab/labnotebook/data/allsilos-tag_and_annot.csv")
#annotations <- annotations[,c(1,67)]
#annotations$Protein.ID <- sub("\\|", ".", annotations$Protein.ID)

#merge <- merge(NSAF.trans, annotations, by.x = "row.names", by.y="Protein.ID", all.x=TRUE)
#merge$Gene.names <- as.character(unlist(merge$Gene.names))
#class(merge$Gene.names)
#merge$Names <- ifelse(merge$Gene.names == "", yes=merge$Row.names, no= merge$Gene.names)
#merge$Names <- ifelse(merge$Names == "None", yes=merge$Row.names, no= merge$Names)
#merge$Names <- ifelse(is.na(merge$Names) == TRUE, yes=merge$Row.names, no=merge$Names)
#any(is.na(merge$Names))

#seperate silos (for merge file)
#s3 <- merge[,c(15,2,4,6,8,10,12)]
#s9 <- merge[,c(15,3,5,7,9,11,13)]
#########################################  For Proteins as names  ################################################################

NSAF.trans$Names <- rownames(NSAF.trans)
rownames(NSAF.trans) <- NULL

#seperate silos
silo3 <- NSAF.trans[, c(grepl("23", colnames(NSAF.trans)))]
silo3$Names <- NSAF.trans$Names

silo9 <- NSAF.trans[, c(grepl("29", colnames(NSAF.trans)))]
silo9$Names <- NSAF.trans$Names

####################################################################################
#add silo/temp prefix into protein/gene names
silo3$Names <- sub("^", "3_", silo3$Names)
silo9$Names <- sub("^", "9_", silo9$Names)

#Make column names time only
colnames(silo3) <- sub("T.*", "", colnames(silo3))
colnames(silo9) <- sub("T.*", "", colnames(silo9))

#combine dataframes
silo3.9 <- rbind(silo3,silo9)

# do below lines if genes are names
#rownames(silo3.9) <-  make.names(silo3.9$ID, unique=TRUE)
#rownames(silo3.9) <- sub("X", "", rownames(silo3.9))

#do below if proteins are names
rownames(silo3.9) <- silo3.9$Names
silo3.9 <- silo3.9[,-c(7)]

#Optional: add in day 0 
#silo3.9$D0 <- NSAF.trans$D0T16
#library(dplyr)
#silo3.9 <- silo3.9 %>% select("D0", everything())

#biostats package
source("Documents/robertslab/labnotebook/analysis/scripts/biostats.R")

##############################################################################################################################

#Data is loaded in as silo3.9 with proteins as row names and biostats is loaded
#Now start cluster analysis!

#dissimilarity matrix and clustering
library(vegan)
nsaf.euc<-vegdist(silo3.9, method='euclidean')

library(cluster)
clust.avg<-hclust(nsaf.euc, method='average')

#agglomerate coefficent
coef.hclust(clust.avg) #bray=0.9349286; euclidean=0.9975303(w/out day 0)

#cophenetic correlation (want at least 0.75)
cor(nsaf.euc, cophenetic(clust.avg)) #bray=0.7690317; euclidean=0.9460521

#Scree plot
jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/euclidean/eu-screeplot.jpeg", width = 1000, height = 1000)
hclus.scree(clust.avg)
dev.off()

#Look for the elbow/inflection point on the scree plot and you can estimate number of clusters. 
#But  it seems that this information cannot be pulled from the scree plot.

#cut dendrogram at selected height
plot(clust.avg, labels=FALSE,xlab = "Protein")
rect.hclust(clust.avg, h=200)


jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/euclidean/eu-dendrogram.jpeg", width = 1000, height = 1000)
plot(clust.avg, labels=FALSE)
rect.hclust(clust.avg, h=250)
rect.hclust(clust.avg, h=200)
rect.hclust(clust.avg, h=150)
rect.hclust(clust.avg, h=100)
dev.off()

#this looks reasonable
clust.class<-cutree(clust.avg, h=200)
max(clust.class) 

#Cluster Freq table
silo3_9.freq <- data.frame(table(clust.class))

write.csv(silo3_9.freq, file = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/euclidean/200freq.csv", row.names = FALSE)

#Make df
silo3_9.clus <- data.frame(clust.class)
names <- rownames(silo3_9.clus)
silo3_9.clus <- cbind(names, silo3_9.clus)
rownames(silo3_9.clus) <- NULL
colnames(silo3_9.clus)[1] <- "ID"
colnames(silo3_9.clus)[2] <- "Cluster"

#merge for abundance values
silo3_9.all <- merge(silo3_9.clus, silo3.9, by.x = "ID", by.y = "row.names")

#this gives matrix of 2 columns, first with proteins second with cluster assignment
#Line plots for each cluster
library(ggthemes)
library(reshape)
library(ggplot2)

melted_all_s3_9<-melt(silo3_9.all, id.vars=c('ID', 'Cluster'))

Temperature <- paste(sep="", 2,(substr(melted_all_s3_9$ID, 0, 1)))

melted_all_s3_9$variable <- gsub("D","", melted_all_s3_9$variable)
melted_all_s3_9$variable <- as.factor(melted_all_s3_9$variable) 
class(melted_all_s3_9$variable)

jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/euclidean/250faceted-abund.jpeg", width = 1000, height = 1000)

ggplot(melted_all_s3_9, aes(x=variable, y=value, group=ID, color=Temperature)) +
  geom_line(alpha=0.8) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + 
  labs(x='Day', y='Normalized Spectral Abundance Factor')

dev.off()

#Seperate Silo from protein name
silo3_9.edit <- silo3_9.all
silo3_9.edit$Silo <- substr(silo3_9.edit$ID,1,1)

#Remove protein silo notation and organize
class(silo3_9.edit$ID)
silo3_9.edit$ID <- as.character(unlist(silo3_9.edit$ID))
silo3_9.edit$ID <- substr(silo3_9.edit$ID, 3, nchar(silo3_9.edit$ID))

library(dplyr)
silo3_9.edit <- silo3_9.edit %>% select(Silo, everything())
silo3_9.edit <- silo3_9.edit %>% select(Cluster, everything())

#save datasheet with abundance values, clusters, and annotations
#silo3_9.annot <- merge(silo3_9.edit, annotations, by.x = "ID", by.y = "Protein.ID")
#write.csv(silo3_9.annot, file = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/techreps-avgs/silo3and9_nozerovals_AVGs/all-prot-anno.csv", row.names = FALSE)

#Parse out unique proteins
library(data.table)
unique.prot <- silo3_9.edit[!(duplicated(silo3_9.edit[c("ID", "Cluster")]) | duplicated(silo3_9.edit[c("ID", "Cluster")], fromLast = TRUE)), ]

#determine if duplicates were removed
anyDuplicated(silo3_9.edit[,c("ID","Cluster")]) #returns first duplicated rows
anyDuplicated(unique.prot[,c("ID","Cluster")]) #returns 0 because there are no duplicates

#datasheet with unique proteins and annotations
#final.unique.prot <- merge(unique.prot, annotations, by.x = "ID", by.y = "Protein.ID")

#ensure protein numbers stayed equal between silos
sum(unique.prot$Silo == "3")
sum(unique.prot$Silo == "9")

write.csv(unique.prot, file = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/euclidean/250unique-prot.csv", row.names = FALSE)

#Now I need to have only one column for each day rather than one column per silo per day
#protein.names <-  data.frame(paste(s39.unq.abudance$Silo, "_", s39.unq.abudance$Protein, sep = ""))
#colnames(protein.names) <- "Protein"
#plot.unq.prot <-  silo3_9.all[which(silo3_9.all$Protein %in% protein.names$Protein),]

#Plot abudances of unique proteins
plot.unq <- unique.prot
plot.unq$S.ID <- paste(plot.unq$Silo, plot.unq$ID, sep="_")
plot.unq <- plot.unq[,-c(2:3)]

unq_melted_all_s3_9<-melt(plot.unq, id.vars=c('S.ID', 'Cluster'))
Temperature <- paste(sep="", 2,(substr(unq_melted_all_s3_9$S.ID, 0, 1)))

unq_melted_all_s3_9$variable <- gsub("D","", unq_melted_all_s3_9$variable)
unq_melted_all_s3_9$variable <- as.factor(as.numeric(unq_melted_all_s3_9$variable))
class(unq_melted_all_s3_9$variable)

jpeg(filename = "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/euclidean/250faceted-unq.jpeg", width = 1000, height = 1000)

ggplot(unq_melted_all_s3_9, aes(x=variable, y=value, group=S.ID, color = Temperature)) +
  geom_line(alpha=0.8) + theme_bw() +
  facet_wrap(~Cluster, scales='free_y') + 
  labs(x='Day', y='Normalized Spectral Abundance Factor')

dev.off()


write.csv(silo3_9.all, "Documents/robertslab/labnotebook/analysis/clustering/silo3_9-NSAF/euclidean/250-clust-protein.csv")
###########################

#if you want more annotations instead of just gene names
#annotations <- read.csv("Documents/robertslab/labnotebook/data/allsilos-tag_and_annot.csv")
#annotations <- annotations[,c(1,63)]
#annotations$Protein.ID <- sub("\\|", ".", annotations$Protein.ID)

#unq.anno <- merge(unique.prot, annotations, by.x="ID", by.y="Protein.ID")

#write.csv(unq.anno, "Documents/robertslab/labnotebook/analysis/clustering/NSAF-s3s9/unqprot-travg-anno.csv")

#library(dplyr)
#count(unq.anno, Entry)

