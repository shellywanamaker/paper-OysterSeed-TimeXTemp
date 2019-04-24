#title: "Clustering technical replicates and calculating average NSAF values"
#author: "Shelly Trigg"
#date: "1/11/2019"

#Load packages

library(vegan)
library(ggplot2)
library(dplyr)
library(gtools)


#Load Abacus data, parse out ADJNSAF values, and simplify column names to just sample number

#upload data file
ABACUSdata <- read.csv("~/Documents/GitHub/paper-OysterSeed-TimeXTemp/Data/ABACUS_output.tsv", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
#select only columns containing ADJNSAF and Protein ID
ABACUSdata <- ABACUSdata[,c(1,grep("ADJNSAF", colnames(ABACUSdata)))]

## change column names in ABACUSdata to just sampleID
colnames(ABACUSdata) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata))
colnames(ABACUSdata) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata))


#Load meta data file with temperature and day information

#upload meta data; this was a csv file I create from Rhonda's notebook entry: https://github.com/Ellior2/Ellior2.github.io/blob/master/_posts/2017-3-11-NMDS-analysis.md
meta_data <- read.csv("~/Documents/GitHub/paper-OysterSeed-TimeXTemp/Data/Sample_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
meta_data$silo <- substr(meta_data$Contents,5,5)
meta_data$day <- substr(meta_data$SampleName,5,6)
meta_data$SampleName <- gsub(pattern = "H","",meta_data$SampleName)
meta_data$SampleName <- gsub(pattern = "C","",meta_data$SampleName)
#create a temperature column
meta_data$temp <- "temp"
for(i in 1:nrow(meta_data)){
  if(meta_data$silo[i] == "2"){
    meta_data$temp[i] <- "23"
  }
  if(meta_data$silo[i] == "3"){
    meta_data$temp[i] <- "23"
  }
  if(meta_data$silo[i] == "9"){
    meta_data$temp[i] <- "29"
  }
  if(meta_data$silo[i] == "e"){
    meta_data$temp[i] <- "16"
  }
}


#Reformat Abacus data for PCA

#Transpose- switch rows and columns
tABACUSdata <- t.data.frame(ABACUSdata[,-1])
colnames(tABACUSdata) <- ABACUSdata[,1]
tABACUSdata <- cbind(data.frame(rownames(tABACUSdata)),tABACUSdata)
colnames(tABACUSdata)[1] <- "SampleID"

#add meta data to abacus data
tABACUSdata <- merge(meta_data[,c(1,2,7,8)],tABACUSdata, by = "SampleID")

#Remove Silo 2 and day 15
silo3and9 <- tABACUSdata[which(substr(tABACUSdata$SampleName,1,2) != "S2" & tABACUSdata$day != "15"),]
#make rownames from Sample ID column so that the NMDS knows what's what
rownames(silo3and9) <- silo3and9$SampleID
#order the data frame by day and temperature so coloring the points on the plot is easier
silo3and9 <- silo3and9[order(as.numeric(silo3and9$day),silo3and9$temp),]
#remove day 0 
silo3and9 <- silo3and9[which(silo3and9$day !=0),]


#Determine if any proteins have zero ADJNSAF vals for all samples; this would be because they were in Silo 2, but not in Silo 3 or 9

no_val_proteins <- silo3and9[,which(apply(silo3and9, 2, var) == 0)]
ncol(no_val_proteins)


#Remove proteins if they have a zero value in all samples

silo3and9_nozerovar <- silo3and9[,-c(1:4,which(colnames(silo3and9) %in% colnames(no_val_proteins)))]
#check to make sure it worked
ncol(silo3and9)-ncol(silo3and9_nozerovar)


#PCA on log transformed values

#For proteins with a zero value in any sample, replace with very small value
silo3and9_nozerovar[silo3and9_nozerovar == 0.0000] <- 0.1000
#log transform NSAF values
silo3and9_log <- log(silo3and9_nozerovar,2)
#create PCA object
pca_log <- prcomp(silo3and9_log, center = F, scale = F)
#add meta data to PCA data
pca_log_meta <- cbind(silo3and9$day, silo3and9$temp, data.frame(paste(silo3and9$day,silo3and9$temp, sep = "_")),pca_log$x)
#rename columns
colnames(pca_log_meta)[1:3] <- c("day","temp","SampleName")
#plot PCA
jpeg("~/Documents/GitHub/paper-OysterSeed-TimeXTemp/Figures/Supplementary/SupplementaryFigure1.jpg")
ggplot(pca_log_meta, aes(PC1, PC2)) + geom_point(aes(col = day, shape = temp)) + theme_bw() + ggtitle("PCA of log ADJNSAF values with zeros replaced with 0.1")
dev.off()

#Calculate the technical replicate average NSAF value for each protein

#create an empty data frame to get filled in by loop
df_all_avg <- data.frame()
#loop through the data and calculate the mean between replicates for each protein
for (i in seq(1,nrow(silo3and9_nozerovar),2)){
  #this calculates the mean for each odd number row and the row following it
  df_all_avg_row <- apply(silo3and9_nozerovar[c(i,i+1),],2,mean)
  #this sequencially appends each row together 
  df_all_avg <- rbind(df_all_avg, df_all_avg_row)
}
#add column names to mean NSAF data
colnames(df_all_avg) <- colnames(silo3and9_nozerovar)
#add Sample ID column by only including technical replicates with numeric names (excluding all technical replicate names that include "A")
df_all_avg$SampleID <- rownames(silo3and9_nozerovar[-grep("A",rownames(silo3and9_nozerovar)),])


#export technical replicate NSAF mean values data for all proteins as columns with zero values converted to 0.1
silo3and9_meta <- silo3and9[,c(1,3,4)]
new_data_all <- merge(silo3and9_meta,df_all_avg, by = "SampleID")
new_data_all <- new_data_all[order(new_data_all[,"day"], new_data_all[,"temp"]),-1]
write.csv(new_data_all, "~/Documents/GitHub/paper-OysterSeed-TimeXTemp/Data/silo3and9_nozerovals_NSAF_AVGs.csv", row.names = FALSE, quote = FALSE)

#export technical replicate NSAF mean values data for all proteins as rows with zero values
new_data_all_t <- data.frame(t(new_data_all[,-c(1:2)]))
colnames(new_data_all_t) <- paste0("Day",new_data_all[,1],"_",new_data_all[,2],"C")
new_data_all_t[new_data_all_t == 0.1] <- 0
new_data_all_t$proteinID <- rownames(new_data_all_t)
write.csv(new_data_all_t,"~/Documents/GitHub/paper-OysterSeed-TimeXTemp/Data/silo3and9_NSAF_AVGs.csv", row.names = FALSE, quote =FALSE )
