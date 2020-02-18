#title: "Clustering technical replicates and calculating average NSAF values"
#author: "Shelly Trigg"
#date: "1/11/2019"

#Load packages

library(ggplot2)
library(dplyr)
library(gtools)


#Load Abacus data, parse out ADJNSAF values, and simplify column names to just sample number

#upload data file
ABACUSdata <- read.csv("../../Data/ABACUS_output.tsv", sep = "\t", header=TRUE, stringsAsFactors = FALSE)
#select only columns containing ADJNSAF and Protein ID
ABACUSdata <- ABACUSdata[,c(1,grep("ADJNSAF", colnames(ABACUSdata)))]

## change column names in ABACUSdata to just sampleID
colnames(ABACUSdata) <- gsub(pattern = "X20161205_SAMPLE_", "", colnames(ABACUSdata))
colnames(ABACUSdata) <- gsub(pattern = "_ADJNSAF", "", colnames(ABACUSdata))


#Load meta data file with temperature and day information

#upload meta data; this was a csv file I create from Rhonda's notebook entry: https://github.com/Ellior2/Ellior2.github.io/blob/master/_posts/2017-3-11-NMDS-analysis.md
meta_data <- read.csv("../../Data/Sample_metadata.csv", header = TRUE, stringsAsFactors = FALSE)
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

#create dpf column
meta_data$dpf <- "dpf"
for(i in 1:nrow(meta_data)){
  if(meta_data$day[i] == "0"){
    meta_data$dpf[i] <- "19"
  }
  if(meta_data$day[i] == "3"){
    meta_data$dpf[i] <- "21"
  }
  if(meta_data$day[i] == "5"){
    meta_data$dpf[i] <- "23"
  }
  if(meta_data$day[i] == "7"){
    meta_data$dpf[i] <- "25"
  }
  if(meta_data$day[i] == "9"){
    meta_data$dpf[i] <- "27"
  }
  if(meta_data$day[i] == "11"){
    meta_data$dpf[i] <- "29"
  }
  if(meta_data$day[i] == "13"){
    meta_data$dpf[i] <- "31"
  }
  if(meta_data$day[i] == "15"){
    meta_data$dpf[i] <- "33"
  }
}  


#Reformat Abacus data for PCA

#Transpose- switch rows and columns
tABACUSdata <- t.data.frame(ABACUSdata[,-1])
colnames(tABACUSdata) <- ABACUSdata[,1]
tABACUSdata <- cbind(data.frame(rownames(tABACUSdata)),tABACUSdata)
colnames(tABACUSdata)[1] <- "SampleID"

#add meta data to abacus data
tABACUSdata <- merge(meta_data[,c(1,2,7,8,9)],tABACUSdata, by = "SampleID")

#Remove Silo 2 and day 15
silo3and9 <- tABACUSdata[which(substr(tABACUSdata$SampleName,1,2) != "S2" & tABACUSdata$day != "15"),]
#make rownames from Sample ID column so that the NMDS knows what's what
rownames(silo3and9) <- silo3and9$SampleID
#order the data frame by day and temperature so coloring the points on the plot is easier
silo3and9 <- silo3and9[order(as.numeric(silo3and9$day),silo3and9$temp),]


#Determine if any proteins have zero ADJNSAF vals for all samples; this would be because they were in Silo 2, but not in Silo 3 or 9

no_val_proteins <- silo3and9[,which(apply(silo3and9, 2, var) == 0)]
ncol(no_val_proteins)


#Remove proteins if they have a zero value in all samples

silo3and9_nozerovar <- silo3and9[,-c(1:5,which(colnames(silo3and9) %in% colnames(no_val_proteins)))]
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
pca_log_meta <- cbind(silo3and9$dpf, silo3and9$temp, data.frame(paste(silo3and9$dpf,silo3and9$temp, sep = "_")),pca_log$x)
#rename columns
colnames(pca_log_meta)[1:3] <- c("dpf","temp","SampleName")

#plot PCA with days ordered from earliest to latest
jpeg("../../Figures/AdditionalFile1/FigureS1.jpg", height = 5, width = 7, units = "in", res = 300 )
ggplot(pca_log_meta, aes(PC1, PC2)) + geom_point(aes(col = dpf, shape = temp)) + theme_bw()
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
silo3and9_meta <- silo3and9[,c(1,5,4)]
new_data_all <- merge(silo3and9_meta,df_all_avg, by = "SampleID")
new_data_all <- new_data_all[order(new_data_all[,"dpf"], new_data_all[,"temp"]),-1]
write.csv(new_data_all, "silo3and9_nozerovals_NSAF_AVGs.csv", row.names = FALSE, quote = FALSE)

#RUN PCA ON AVERAGE NSAF VALUES
#log transform NSAF values
new_data_all_log <- log(new_data_all[,-c(1:2)],2)
#create PCA object
pca_log <- prcomp(new_data_all_log, center = T, scale = T)
summary_pca_log <- summary(pca_log)

#add meta data to PCA data
pca_log_meta <- cbind(new_data_all$dpf, new_data_all$temp, data.frame(paste(new_data_all$dpf,new_data_all$temp, sep = "_")),pca_log$x)
#rename columns
colnames(pca_log_meta)[1:3] <- c("dpf","temp","SampleName")

#plot PCA of averaged NSAF values with days ordered from earliest to latest
jpeg("../../Figures/MainText/Figure2b.jpg", height = 5, width = 7, units = "in", res = 300 )
ggplot(pca_log_meta, aes(PC1, PC2,label = dpf, color = temp, size = 2)) + geom_text(fontface = "bold") + scale_colour_manual(values = c("brown4","cyan3", "magenta3")) + theme_bw() + ylab(paste("PC2"," (",formatC(summary_pca_log$importance[5] * 100,digits=2,format="f"),"%)", sep = "")) + xlab(paste("PC1"," (",formatC(summary_pca_log$importance[2] * 100,digits=2,format="f"),"%)", sep = ""))
dev.off()

#extract loadings for PC1 and PC2
PC1_great <- data.frame(pca_log$rotation)
PC1_great <- PC1_great[which(PC1_great$PC1 > 0),1:2]

#export plot
jpeg("../../Figures/AdditionalFile1/FigureS2a.jpg", height = 5, width = 5, units = "in", res = 300 )
plot(head(PC1_great[order(PC1_great, decreasing = TRUE),1],100), xlab = "protein", ylab = "PC1 Loadings")
abline(h=0.0236, col = "red")
dev.off()


PC1_great_prots <- rownames(PC1_great[which(PC1_great$PC1 > 0.0236),])

PC1_less <- data.frame(pca_log$rotation)
PC1_less <- PC1_less[which(PC1_less$PC1 < 0),1:2]

jpeg("../../Figures/AdditionalFile1/FigureS2b.jpg", height = 5, width = 5, units = "in", res = 300 )
plot(tail(PC1_less[order(PC1_less$PC1, decreasing = TRUE),1],100), xlab = "protein", ylab = "PC1 Loadings")
abline(h=-0.02355, col = "red")
dev.off()

PC1_less_prots <- rownames(PC1_less[which(PC1_less$PC1 < -0.02355),])

#now PC2
PC2_great <- data.frame(pca_log$rotation)
PC2_great <- PC2_great[which(PC2_great$PC2 > 0),1:2]

jpeg("../../Figures/AdditionalFile1/FigureS2c.jpg", height = 5, width = 5, units = "in", res = 300 )
plot(head(PC2_great[order(PC2_great$PC2, decreasing = TRUE),2],100), xlab = "protein", ylab = "PC2 Loadings")
abline(h=0.0245, col = "red")
dev.off()

PC2_great_prots <- rownames(PC2_great[which(PC2_great$PC2 > 0.0245),])

PC2_less <- data.frame(pca_log$rotation)
PC2_less <- PC2_less[which(PC2_less$PC2 < 0),1:2]

jpeg("../../Figures/AdditionalFile1/FigureS2d.jpg", height = 5, width = 5, units = "in", res = 300 )
plot(tail(PC2_less[order(PC2_less$PC2, decreasing = TRUE),2],100), xlab = "protein", ylab = "PC2 Loadings")
abline(h=-0.0262, col = "red")
dev.off()

PC2_less_prots <- rownames(PC2_less[which(PC2_less$PC2 < -0.0262),])


global_PC12 <- c(PC1_great_prots,PC1_less_prots, PC2_great_prots, PC2_less_prots)
write.table(global_PC12, "../../Data/global_PC12_proteins.txt", quote = FALSE, row.names = FALSE)





