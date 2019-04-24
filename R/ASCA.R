#load libraries
library(dplyr)
library(tidyr)
library(MetStaT)
library(ggplot2)
library(heatmap3)


#Read in NSAF data (this is Supplementary Data 3)
data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/nmds_R/silo3and9_nozerovals_AVGs.csv", stringsAsFactors = FALSE)
#make a dataframe with just CHOYP proteins
choyp_data <- data[which(data$day !=0),c(1:2,grep("CHOYP", colnames(data)))]


#create matrix to pass to ASCA command, excluding the silo and time info
ASCA_X <- as.matrix(choyp_data[,-c(1:2)])
#create matrix to pass to ASCA command with only the silo and time info
ASCA_F <- as.matrix(choyp_data[,c(1:2)])
#perform ASCA
ASCA <- ASCA.Calculate(ASCA_X, ASCA_F, equation.elements = "1,2,12", scaling = FALSE)

ASCA.GetSummary(ASCA)

### Plot PCAs from ASCA; This first plot is the time (days) effect PCA**
#plot PCA for factor 1, which is time in this case
ASCA.PlotScoresPerLevel(ASCA, ee = "1", pcs = "1,2")


##This next plot is the temperature effect PCA**
#plot PCA for factor 2, which is temperature in this case
ASCA.PlotScoresPerLevel(ASCA, ee = "2", pcs = "1,2")


##This next plot is the time x temp interaction effect PCA**
#plot PCA for factor interaction, which is time x temp in this case
timextemp_PC12 <- data.frame(ASCA$`12`$svd$t[,c(1,2)])
timextemp_PC12 <- cbind(data.frame(ASCA$`12`$level.combinations$row.patterns), timextemp_PC12)
colnames(timextemp_PC12)<- c("day","temp","PC1","PC2")
timextemp_PC12$day <- as.character(timextemp_PC12$day)
timextemp_PC12$temp <- as.character(timextemp_PC12$temp)
ggplot(timextemp_PC12, aes(PC1, PC2)) + geom_point(aes(col = day, shape = temp, size = 3)) + theme_bw() + ggtitle("PC1 vs PC2 for time x temperature interaction effect") + theme(plot.title = element_text(face = "bold")) + xlab(paste("PC1"," (",formatC(ASCA$`12`$svd$var.explained[1] * 100,digits=2,format="f"),"%)", sep = "")) + ylab(paste("PC2"," (",formatC(ASCA$`12`$svd$var.explained[2] * 100,digits=2,format="f"),"%)", sep = ""))

### Analysis of proteins affected by temperature
#Because the temperature effect PCA show the most separation between 23C and 29C in PC2, we will look at those loadings.
#PC2 loadings for temperature effect
#extract protein names from ASCA data; these will be combined with loadings values; the order is maintained
protnames <- data.frame(colnames(ASCA$data))
#combine protein names with ASCA loadings for PC2 for temperature, since this component showed the greatest separation between temperatures
d <- cbind(protnames, ASCA$`2`$svd$v[,1])
#rename the columns
colnames(d) <- c("protein", "PC1loadings")
#make a dataframe of proteins with PC loadings greater than zero for the loadings plot
d_great <- d[which(d$PC1loadings > 0),]
#make a dataframe of proteins with PC loadings less than zero for the loadings plot
d_less <- d[which(d$PC1loadings < 0),]
#plot PC1 loadings
plot(d_great[order(d_great$PC1loadings, decreasing = TRUE),2], xlab = "protein", ylab = "PC1 Loadings")
abline(h=0.03, col = "red")
plot(d_less[order(d_less$PC1loadings, decreasing = TRUE),2],xlab = "protein", ylab = "PC1 Loadings")
abline(h=-0.025, col = "red")

#To pull out proteins affected by temperature based on their influence in seperating treatment groups on PC2 of the temperature PCA, I picked an absolute value loadings threshold of 0.025. This means any protein that had a loadings value > 0.025 or < -0.025 was selected.
#make list of proteins with temperature PC1 loadings values >= 0.025
cutd <- d[which(d$PC1loadings >= 0.03 | d$PC1loadings <= -0.025),]
#make a list of cutoff proteins with normalized abundance
cut_data <-choyp_data[,which(colnames(choyp_data) %in% cutd$protein)]
rownames(cut_data) <- paste(choyp_data$day, choyp_data$temp, sep="_")

write.csv(cut_data, "~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_all_proteins_avgADJNSAF/ASCA_avgNSAFvals_allProteins_noDay0/ASCA_TempAffectedProteins.csv")
write.csv(cutd,"~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_all_proteins_avgADJNSAF/ASCA_avgNSAFvals_allProteins_noDay0/ASCA_TempAffectedProteins_loadings.csv", row.names = FALSE, quote =FALSE)

#Number of proteins affected by temperature at loadings value > 0.03 or < -0.025

nrow(cutd)



#Heatmap of proteins affected by temperature based on temperature effect PC2 loadings value cutoff
#reorder the cut data so that samples are grouped by temperature
cut_data_ord <- data[,c(1:3,which(colnames(data) %in% cutd$protein))]
cut_data_ord <- cut_data_ord[order(cut_data_ord$temp,cut_data_ord$day),]
rownames(cut_data_ord) <- paste(cut_data_ord$day, cut_data_ord$temp, sep="_")
cut_data_ord <- cut_data_ord[,-c(1:3)]
#replace CHOYP names with "entry names"
#transposed to eventually get row names as columns to match entry names to. 
#t.data.frame returns a matrix so I had convert back to data frame for all downstream rearranging because data frames are much easier to work with
cut_data_ord_t <- data.frame(t.data.frame(cut_data_ord))
#make row names which contain the protein IDs a column
cut_data_ord_t <- cbind(rownames(cut_data_ord_t), cut_data_ord_t)
#add the column name
colnames(cut_data_ord_t)[1] <- "Protein.ID"
#convert the Protein ID column to character
cut_data_ord_t[,1] <- as.character(cut_data_ord_t[,1])
#add the Entry names
cut_data_ord_t <- merge(data_w_uniprot[,c("Protein.ID","Entry.name.new")],cut_data_ord_t, by = "Protein.ID")
#make the entry names the row names
rownames(cut_data_ord_t) <- cut_data_ord_t[,"Entry.name.new"]
#remove the protein ID and entry name columns since the row names are added
cut_data_ord_t <- cut_data_ord_t[,-c(1:2)]
#this step removes the x from the column names while preserving the order of colnames
colnames(cut_data_ord_t) <- rownames(cut_data_ord)
#plot the heat map
heatmap3(as.matrix(cut_data_ord_t), Colv = NA, cexRow = 0.5)
