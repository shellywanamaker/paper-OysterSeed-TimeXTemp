#load libraries
library(plyr)
library(tidyr)
#install.packages("ontologyIndex")
#install.packages("ontologySimilarity")
library(ontologyIndex)
library(ontologySimilarity)
library(GSEABase)
library(reshape2)


#Cytoscape needs two files to build a network: 1. Node Attribute file and 2. Edge attribute file. 
#Node attribute file can be a list of proteins with information about each protein like alternative names, fold change and pvalue information, etc. 
#All of this data can be used to change the appearace of the nodes in the network. 
#For example, you can color nodes by their foldchange, node sizes can be based on their p-value, etc.

##########################################################
###Getting node attributes containing same day comparisons
###########################################################

#read in uniprot mapping
uniprot <- read.csv("/Volumes/web/metacarcinus/Cgigas/all_giga-uniprot-blastP-out.nopipe.annotations.tab", sep ="\t", header = FALSE, stringsAsFactors = FALSE)
#rename uniprot columns
colnames(uniprot) <- c("protein_ID","Entry", "Entry_name", "perc_ident_match", "align_len", "num_mismatch", "num_gaps","querStart", "querEnd", "subjStart", "subjEnd", "evalue", "bitscore","Entry.1","Entry_name.1", "Protein_names", "Gene_names", "Organism", "Protein_length","Pathway", "GO_bp", "GO","GO_IDs", "Protein_fams")
#convert columns that should be numeric to numeric
str(uniprot)
uniprot[,4:13] <- lapply(uniprot[,4:13], as.numeric)
#NAs introduced by coercion
#this is because the entries that are unmappaed have NA values
str(uniprot)

#read in same day log FC and pval data
sameday_logFC_pval <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/sumNUMSPECSTOT_plus1_ratioFC_logFC_pval_DAYSCOMPARED.csv", stringsAsFactors = FALSE)
colnames(sameday_logFC_pval)[1] <- "protein_ID"
#combine uniprot and foldchange data
sameday_logFC_pval_uniprot<- merge(sameday_logFC_pval, uniprot, by = "protein_ID", all.x = TRUE)
#exclude proteins that didn't map to uniprot DB
sameday_logFC_pval_uniprot_mapped <- sameday_logFC_pval_uniprot[-grep("unmapped", sameday_logFC_pval_uniprot$Entry),]

#make a list of proteins that mapped to uniprot DB @ eval of < 10^10
proteins_evalpass <- sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped$evalue <= 10^-10),"protein_ID"]
write.csv(proteins_evalpass, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/proteins_evalpass.csv", quote =FALSE, row.names = FALSE)
#save data frame without proteins that have evalue > 10^-10
sameday_logFC_pval_uniprot_mapped <- sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped$evalue <= 10^-10),]


#####select only proteins with adj Chi sq. pvalue <= 0.1####

#create a list of all column names with adj.Chisq.pval 
#adjChiSqpvalColumns <- colnames(sameday_logFC_pval_uniprot_mapped)[grep("adj.ChiSq.pval",colnames(sameday_logFC_pval_uniprot_mapped))]

#build a list of protiens with adj Chi sq. pvalue <= 0.1
#create empty data frame the loop will add too
#all_sig_pro <- data.frame()
#for (i in 1:length(adjChiSqpvalColumns)){ # for each name in adj.Chisq.pval column name list
#  column <- adjChiSqpvalColumns[i] # create a variable for indexed column name
  #make a data frame containing protein IDs for all proteins in indexed column that have adj.Chisq.pval <=0.1
#  sig_pro <- data.frame(sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped[,column] <= 0.1),1],stringsAsFactors = FALSE)
  #iteratively add protein lists to initial data frame
#  all_sig_pro <- rbind(all_sig_pro, sig_pro)
#}

#count how many unique proteins are in the list of proteins with adj Chi sq. pvalue <= 0.1
#nrow(unique(all_sig_pro))
#[1] 153

#make a data frame of just unique proteins so we can select these from the foldchange/pvalue data
#all_sig0.1_pro <- unique(all_sig_pro)
#colnames(all_sig0.1_pro)[1]<- "protein_ID"
#all_sig0.1_pro$Chi <- "PropTest"

#read in ASCA data
#ASCA_tempdata <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/ASCA/ASCA_all_proteins_avgADJNSAF/ASCA_avgNSAFvals_allProteins_noDay0/ASCA_TempAffectedProteins_loadings.csv", stringsAsFactors = FALSE)
#colnames(ASCA_tempdata)[1]<- "protein_ID"
#ASCA_tempdata$ASCA <- "ASCA"
#ASCA_tempdata <- ASCA_tempdata[,-2]

#all_sig0.1_ASCA_pro <- merge(all_sig0.1_pro, ASCA_tempdata, by = "protein_ID", all = TRUE)
#all_sig0.1_ASCA_pro$method <- paste(all_sig0.1_ASCA_pro$Chi, all_sig0.1_ASCA_pro$ASCA, sep = "")
#all_sig0.1_ASCA_pro$method <- gsub("NA","",all_sig0.1_ASCA_pro$method)
#select sig proteins @ p.adj 0.1 from logFC list
#all_sig0.1_pro_logFC_pval <- merge(sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped$evalue <= 10^-10),],all_sig0.1_ASCA_pro[,-c(2,3)], by = "protein_ID")


#select sig proteins @ p.adj 0.1 from logFC list
#all_sig0.1_pro_logFC_pval <- merge(sameday_logFC_pval_uniprot_mapped[which(sameday_logFC_pval_uniprot_mapped$evalue <= 10^-10),],all_sig0.1_pro)

#read in data from stats verification 
sig0.1_ASCA_clustering_data <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/VerifyStatsProteinSelection/all_sig0.1_ASCA_clust_pro.csv", stringsAsFactors = FALSE)
nrow(sig0.1_ASCA_clustering_data)
#merge with uniprot and same day fc data
all_sig0.1_pro_logFC_pval <- merge(sameday_logFC_pval_uniprot_mapped,sig0.1_ASCA_clustering_data, by = "protein_ID")

#read in avg NSAF same day logFC data
avgNSAF_logFC <- read.csv("~/Documents/GitHub/OysterSeedProject/analysis/TotNumSpecRatio_FC_Pval/avgADJNSAF_logFC_DAYSCOMPARED.csv", stringsAsFactors = FALSE)

#merge avgNSAF_logFC with all_sig0.1_pro_logFC_pval

all_sig0.1_pro_logFC_pval <- merge(all_sig0.1_pro_logFC_pval, avgNSAF_logFC, by = "protein_ID")

##Export protein list as node attributes table to upload in cytoscape
#remove extra columns (e.g. bit score, map length, GO terms)

all_sig0.1_pro_logFC_pval$gene_root <- gsub("_.*","", all_sig0.1_pro_logFC_pval$Entry_name)

gene_freq <- data.frame(table(gsub("_.*","", all_sig0.1_pro_logFC_pval$Entry_name)))
gene_freq <- gene_freq[which(gene_freq$Freq > 1),]

#remove choyp_gene part of protein_ID
vers <- gsub("^.*?\\.","", all_sig0.1_pro_logFC_pval$protein_ID)

#remove everything after ".m."
#vers <- gsub("\\.m\\..*$","", vers)

#add version number to gene name
all_sig0.1_pro_logFC_pval$gene_vers <- paste(all_sig0.1_pro_logFC_pval$gene, vers, sep = ".")

#only append version numbers to genes that are in the list more than once
for(i in 1:nrow(all_sig0.1_pro_logFC_pval)){
  if(all_sig0.1_pro_logFC_pval$gene_root[i] %in% gene_freq$Var1){
    all_sig0.1_pro_logFC_pval$gene[i] <- paste(all_sig0.1_pro_logFC_pval$gene_root[i], vers[i], sep = ".")
  }
  else{
    all_sig0.1_pro_logFC_pval$gene[i] <- all_sig0.1_pro_logFC_pval$gene_root[i]
  }
}


#all_sig0.1_pro_logFC_pval <- merge(all_sig0.1_pro_logFC_pval, all_sig0.1_pro_logFC_pval_abbrv[,c("protein_ID", "gene")], by = "protein_ID")

#remove extra columns and non-adj.pvals
all_sig0.1_pro_logFC_pval_abbrv <- all_sig0.1_pro_logFC_pval[,-c(1,20:46,48,50,51,grep("_ChiSq",colnames(all_sig0.1_pro_logFC_pval)))]
#add column to describe node type for cytoscape
all_sig0.1_pro_logFC_pval_abbrv$type <- "protein"


#save node attribute file
write.table(all_sig0.1_pro_logFC_pval_abbrv,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSq_ASCA_clust/all_sig0.1ASCAclust_pro_logFC_pval_abbrv_Evalcutoff_NodeAttb.txt", sep = "\t",quote = FALSE, row.names = FALSE)



############################################################

##########################################################
###Getting edge attributes; these are the GO term relationships
###########################################################
####I want to use GO slim terms instead of GO terms so the network is less busy
###So I need to map my GO terms to GO slims while keeping the protein info.
###That way, I am able to incorporate the protein data into the networks

#Get GO IDs from all_sig0.1_pro_logFC_pval
sig0.1_pro_GO <- all_sig0.1_pro_logFC_pval[,c("gene","GO_IDs")]
sig0.1_pro_GOid_term <- data.frame()
for (i in 1:nrow(sig0.1_pro_GO)){
  sig0.1_pro_GOid_term_row <- data.frame(t(data.frame(strsplit(as.character(sig0.1_pro_GO$GO_IDs[i]),'; ', fixed = TRUE))))
  sig0.1_pro_GOid_term <- rbind.fill(sig0.1_pro_GOid_term,sig0.1_pro_GOid_term_row)
}

#add protein IDs back to GO IDs
sig0.1_pro_GOid_term <- cbind(all_sig0.1_pro_logFC_pval[,"gene"], sig0.1_pro_GOid_term)
#this results in a table with protein ID listed in one column and each next column contains a GO ID that was listed with the protein in Uniprot DB.
#convert all factors to  characters
sig0.1_pro_GOid_term <- data.frame(lapply(sig0.1_pro_GOid_term, as.character), stringsAsFactors = FALSE)
#there are two proteins that don't have GO IDs

#reshape data so that all GO ID columns are gathered in one column called "GO" 
STACKED_sig0.1_pro_GOid_term <- tidyr::gather(sig0.1_pro_GOid_term,"gene","GO", 2:ncol(sig0.1_pro_GOid_term))
#exlude middle column which just contains the string "gene" in each row
STACKED_sig0.1_pro_GOid_term <- STACKED_sig0.1_pro_GOid_term[,c(1,3)]
#remove duplicate rows
STACKED_sig0.1_pro_GOid_term <- unique(STACKED_sig0.1_pro_GOid_term)
colnames(STACKED_sig0.1_pro_GOid_term)[1] <- "gene"
#remove any rows where GO column has NA value. 
STACKED_sig0.1_pro_GOid_term <- STACKED_sig0.1_pro_GOid_term[which(!is.na(STACKED_sig0.1_pro_GOid_term$GO)),]
#this resulting data frame has two columns "gene" and "GO"
write.csv(STACKED_sig0.1_pro_GOid_term, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSq_ASCA_clust/intermediate_files/STACKED_sig0.1_pro_GOid_term.csv", quote = FALSE, row.names = FALSE)
###Next map all GO IDs to GO slim terms


#make list of unique GO terms without a protein ID column
sig0.1_sig_GOids <- unique(STACKED_sig0.1_pro_GOid_term$GO)
#Use GSEA to generate list of all GO Slim BP, MP, and CC
#DF will have "term", "GOid", "GOcategory"

#BP first
#goslims with GSEA
myCollection <- GOCollection(sig0.1_sig_GOids)
#I downloaded goslim_generic.obo from http://geneontology.org/docs/go-subset-guide/
#then i moved it to the R library for GSEABase in the extdata folder
fl <- system.file("extdata", "goslim_generic.obo", package="GSEABase")
slim <- getOBOCollection(fl)
slims <- data.frame(goSlim(myCollection, slim, "BP"))
slims$GOid <- rownames(slims)
slims$Term <- as.character(slims$Term)
rownames(slims) <- NULL

GSEA_bp <- slims[,c("Term", "GOid")]
GSEA_bp$GOcategory <- "BP"

#MP next
slims <- data.frame(goSlim(myCollection, slim, "MF"))
slims$GOid <- rownames(slims)
slims$Term <- as.character(slims$Term)
rownames(slims) <- NULL

GSEA_MF <- slims[,c("Term", "GOid")]
GSEA_MF$GOcategory <- "MF"

GSEA_BP_MF <- rbind(GSEA_bp, GSEA_MF)
colnames(GSEA_BP_MF) <- c("Term", "GOslim", "GOcategory")

#next make a list of all GO ids and any other GO id that is related to it
# first make a list of all ancestor terms that correspond to GO IDs in my list

#load GO data from OntologyX package
data(go)

###this creates a list of ancestor GO IDs for each GO ID in my list
##if the GO ID is not in the "go" object from the package, the entry is "NULL"
#https://jonlefcheck.net/2013/05/20/continuing-a-for-loop-in-r-after-an-error/

term_prop <- list() # create an empty list that will get filled in by loop
for(i in 1:length(sig0.1_sig_GOids)){ # for each line in my GO IDs list
  temp_term_prop <- try(go$ancestors[[sig0.1_sig_GOids[i]]], TRUE) #make a list of all ancester GO IDs for each GO ID in my list
  if(isTRUE(class(temp_term_prop)=="try-error")) {next} else {term_prop[[i]] = temp_term_prop} # if the "go" data doesn't contain my GO ID, go on to the next GO ID. 
}

#create an empty data frame the length the list of ancestor list made above and with two columns
ancestors <- data.frame(matrix(0,length(term_prop),2))
#names the two columns
colnames(ancestors) <- c("orig.GO","GOan")

#make an empty data frame to get filled in by the loop
ances_test <- data.frame()
for(i in 1:length(term_prop)){ #for each GO ID in the ancestors list (which is the same length and order as the sig0.1_sig_GOids list)
  ancestors$orig.GO[i] <- sig0.1_sig_GOids[i] # fill in the orig. GO ID column 
  ancestors$GOan[i] <- paste(term_prop[[i]], collapse = "_") # fill in the GO ancestor column with all GO ancestor IDs corresponding to original GO term separated by an underscore
  ancestors_row <- data.frame(t(data.frame(strsplit(as.character(ancestors$GOan[i]),'_', fixed = TRUE)))) #spread ancestor IDs out across multiple columns
  ances_test <- rbind.fill(ances_test,ancestors_row) #add each row to the data frame
}
#convert factors to characters
ances_test <- data.frame(lapply(ances_test, as.character), stringsAsFactors = FALSE)
#add original GO IDs back to ancestor GO IDs
ances_test <- cbind(data.frame(sig0.1_sig_GOids, stringsAsFactors = FALSE), ances_test)

###list GO IDs not in 'go' object (don't have ancestors)
length(ances_test[which(is.na(ances_test$X1)),1])
#12
ances_test[which(is.na(ances_test$X1)),1]
# [1] "GO:0062023" "GO:0103025" "GO:0102102" "GO:0103046" "GO:0106077" "GO:0061842" "GO:0102131" "GO:0090736" "GO:1905905"
#[10] "GO:0061844" "GO:1905907" "GO:0106036"

#reshape data so that all ancestor GO ID columns are gathered in one column called "Ancterm" 
STACKED_sig0.1_anc_term <- tidyr::gather(ances_test,"Ancestor","Ancterm", 2:ncol(ances_test))
STACKED_sig0.1_anc_term <- STACKED_sig0.1_anc_term[,c(1,3)]
STACKED_sig0.1_anc_term <- unique(STACKED_sig0.1_anc_term)
STACKED_sig0.1_anc_term <- STACKED_sig0.1_anc_term[which(!is.na(STACKED_sig0.1_anc_term$Ancterm)),]

# repeat the above code but this time to make a list of all child terms that correspond to GO IDs in my list
term_propC <- list()
for(i in 1:length(sig0.1_sig_GOids)){
  temp_term_propC <- try(go$children[[sig0.1_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_term_propC)=="try-error")) {next} else {term_propC[[i]] = temp_term_propC}
}

children <- data.frame(matrix(0,length(term_propC),2))
colnames(children) <- c("orig.GO","GOch")

child_test <- data.frame()
for(i in 1:length(term_propC)){
  children$orig.GO[i] <- sig0.1_sig_GOids[i]
  children$GOch[i] <- paste(term_propC[[i]], collapse = "_")
  children_row <- data.frame(t(data.frame(strsplit(as.character(children$GOch[i]),'_', fixed = TRUE))))
  child_test <- rbind.fill(child_test,children_row)
}
#convert factors to characters
child_test <- data.frame(lapply(child_test, as.character), stringsAsFactors = FALSE)
#add original GO IDs back to child GO IDs
child_test <- cbind(data.frame(sig0.1_sig_GOids, stringsAsFactors = FALSE), child_test)

###count GO IDs not in 'go' object (these GO IDs don't have child terms)
length(child_test[which(is.na(child_test$X1)),1])
#621

#reshape data so that all child GO ID columns are gathered in one column called "Chterm" 
STACKED_sig0.1_ch_term <- tidyr::gather(child_test,"Child","Chterm", 2:ncol(child_test))
STACKED_sig0.1_ch_term <- STACKED_sig0.1_ch_term[,c(1,3)]
STACKED_sig0.1_ch_term <- unique(STACKED_sig0.1_ch_term)
STACKED_sig0.1_ch_term <- STACKED_sig0.1_ch_term[which(!is.na(STACKED_sig0.1_ch_term$Chterm)),]


# repeat the above code but this time to make a list of all parent terms that correspond to GO IDs in my list
term_propP <- list()
for(i in 1:length(sig0.1_sig_GOids)){
  temp_term_propP <- try(go$parents[[sig0.1_sig_GOids[i]]], TRUE)
  if(isTRUE(class(temp_term_propP)=="try-error")) {next} else {term_propP[[i]] = temp_term_propP}
}

parents <- data.frame(matrix(0,length(term_propP),2))
colnames(parents) <- c("orig.GO","GOpar")

par_test <- data.frame()
for(i in 1:length(term_propP)){
  parents$orig.GO[i] <- sig0.1_sig_GOids[i]
  parents$GOpar[i] <- paste(term_propP[[i]], collapse = "_")
  parents_row <- data.frame(t(data.frame(strsplit(as.character(parents$GOpar[i]),'_', fixed = TRUE))))
  par_test <- rbind.fill(par_test,parents_row)
}

#convert factors to characters
par_test <- data.frame(lapply(par_test, as.character), stringsAsFactors = FALSE)
#add original GO IDs back to ancestor GO IDs
par_test <- cbind(data.frame(sig0.1_sig_GOids, stringsAsFactors = FALSE), par_test)

###count GO IDs not in 'go' object (these GO IDs don't have parent terms)
length(par_test[which(is.na(par_test$X1)),1])
#638

#reshape data so that all child GO ID columns are gathered in one column called "Parterm" 
STACKED_sig0.1_par_term <- tidyr::gather(par_test,"Parent","Parterm", 2:ncol(par_test))
STACKED_sig0.1_par_term <- STACKED_sig0.1_par_term[,c(1,3)]
STACKED_sig0.1_par_term <- unique(STACKED_sig0.1_par_term)
STACKED_sig0.1_par_term <- STACKED_sig0.1_par_term[which(!is.na(STACKED_sig0.1_par_term$Parterm)),]

#change colnames to match before we bind these data frames together
colnames(STACKED_sig0.1_par_term)[2] <- "term"
colnames(STACKED_sig0.1_ch_term)[2] <- "term"
colnames(STACKED_sig0.1_anc_term)[2] <- "term"

#create an additional column to describe the term
STACKED_sig0.1_par_term$term_type <- "parent"
STACKED_sig0.1_ch_term$term_type <- "child"
STACKED_sig0.1_anc_term$term_type <- "ancestor"

#create another data frame with the original GO IDs in case some of these are in fact already GO slim IDs
sig0.1_orig_term <- cbind(as.data.frame(sig0.1_sig_GOids), as.data.frame(sig0.1_sig_GOids))
colnames(sig0.1_orig_term)[2] <- "term"
sig0.1_orig_term$term_type <- "orig"

#combine all ancestor, child, parent, and original GO IDs into one data frame
par_ch_anc <- rbind(STACKED_sig0.1_anc_term, STACKED_sig0.1_ch_term, STACKED_sig0.1_par_term, sig0.1_orig_term)
#this creates a data frame with three columns:
# sig0.1_sig_GOIDs, which are the original GO IDs from my list
# term, which are the Ancestor, parent, child, or original GO IDs
# type, which specifies "Ancestor", "parent", "child" or, "orig" for each GO ID

nrow(par_ch_anc)
#[1] 21365

# ###################################
# ###Attempt to improve GO ID mapping to GSEA GO slim IDs by including the GO slim alt.id. This would help in the event that an ancestor, parent, child, or the GO id itself matches a GO slim alt.id rather than a GO slim ID. 
#this code is commented out because it didn't help improve the mapping; no child, parent, ancestor or original ID mapped to an alt. id.
# ##################################
# 
# #Get data
# #I downloaded goslim_generic.obo from http://www.geneontology.org/GO_slims/goslim_generic.obo
# 
# #build an index where one column is the GO slim ID and the second column is the GO slim alt. id; The GO slim ID can be listed in the first column multiple times if it has multiple GO slim alt. ids.
# #made an alt id file in command line by head -2637 ~/Downloads/goslim_generic.obo | grep "id\:" > 
# ~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1_ASCA/intermediate_files/GOslimsIDs_alts.csv; then opening in excel and adding the slim term to alt terms in another column
# ##I probably could have made a loop in R to do this, but in the interest of time I did not
# 
# 
# #read in data
# GOslims_alts <- read.csv("~/Desktop/GOslimsIDs_alts.csv", stringsAsFactors = FALSE)
# colnames(GOslims_alts)[1] <- "term"
# par_ch_anc_GOslims_alts <- merge(unique(par_ch_anc[,1:2]), GOslims_alts, by = "term")
# par_ch_anc_GOslims_uniq <- unique(par_ch_anc_GOslims_alts[,2:3])
# 
# #create a frequency table of GO slim IDs in par_ch_anc_slimBP data frame
# par_ch_anc_GOslims_uniq_freq <- data.frame(table(par_ch_anc_GOslims_uniq$GOslim)) 
# colnames(par_ch_anc_GOslims_uniq_freq)[1] <- "GOid" #rename the first column
# par_ch_anc_GOslims_uniq_freq$GOid <- as.character(par_ch_anc_GOslims_uniq_freq$GOid) #change the class of the first column from factor to character
# 
# #view the par_ch_anc_slim_freq table with the slims table to compare how many original GO ids map to GO Slim IDs
# #View(merge(slims[,c("GOid", "Count")], par_ch_anc_GOslims_uniq_freq, by = "GOid", all.x = TRUE))
# 
# #par_ch_anc_GOslims_uniq_cat <- merge(par_ch_anc_GOslims_uniq, GSEA_BP_MF, by = "GOslim")
# 
# ###I don't get more mapping from including the alternate IDs; maybe this is because they are obsolete IDs
# ##################################################


#extract only unqiue "all terms" (e.g. if an original term and an ancestor term are the same, don't list twice)
par_ch_anc_uniq <- unique(par_ch_anc[,c("sig0.1_sig_GOids","term")])
nrow(par_ch_anc_uniq)
#18417


#rename term column so merge will work (REMEMBER this term column is a list of all terms ever for each sig0.1_sig_GOids)
colnames(par_ch_anc_uniq)[2] <- "GOslim"

par_ch_anc_uniq_BP_MF <- merge(par_ch_anc_uniq, GSEA_BP_MF, by = "GOslim")
#count unique GO ids remaining
length(unique(par_ch_anc_uniq_BP_MF$sig0.1_sig_GOids))
#912

colnames(par_ch_anc_uniq_BP_MF) <- c("GOslim", "GO", "GOslimTerm", "GOcategory")

#merge par_ch_anc_slimBP and STACKED_sig0.1_pro_GOid_term and remove proteins with 'Biological Process' GO slim term since these are uninformative
goslim_protein <- merge(STACKED_sig0.1_pro_GOid_term,par_ch_anc_uniq_BP_MF, by = "GO")
#count the number of unique proteins that have GO Slim terms
length(unique(goslim_protein$gene))
#[1] 261

#take a look a proteins that did not map to GO slim terms
all_sig0.1_pro_logFC_pval[which(!(all_sig0.1_pro_logFC_pval$gene %in% (unique(goslim_protein$gene)))),]
#there are 23 of them and they either didn't have a GO id or only has a cell component term


#how many proteins have "biological process" term (which is a really vague, non-descriptive term)
nrow(goslim_protein[grep("GO:0008150|GO:0003674|GO:0005575", goslim_protein$GOslim),])
#[1] 1879

#remove vague terms
goslim_protein <- goslim_protein[-grep("GO:0008150|GO:0003674|GO:0005575", goslim_protein$GOslim),]

########
##Using GO semantic similarity to define relationships among GO slim terms. This is similar methodology to REVIGO
#####

#create a list of just unique GO Slim IDs excluding the biological process ID GO:0008150
sig0.1_GOslims <- unique(goslim_protein$GOslim)
#output this list so that we can upload it to REVIGO to see how these terms can be further slimmed/categorized
#write.csv(sig0.1_GOslims, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSqPval0.1proteins_GOsemsim_edges/sig0.1_GOslimIDs.csv", quote = FALSE, row.names = FALSE)

#make a list of GO slim IDs in a format can be used in the OntologySimilarity function
OX_sig0.1_GOslims <- list()
for(i in 1:length(sig0.1_GOslims)){
  temp_OX_sig0.1_GOslims <- try(go$id[[sig0.1_GOslims[i]]], TRUE)
  if(isTRUE(class(temp_OX_sig0.1_GOslims)=="try-error")) {next} else {OX_sig0.1_GOslims[[i]] = temp_OX_sig0.1_GOslims}
}

#creating a GO_IC formatted file for the OntologySimilarity function
GO_IC_data <- get_term_info_content(go, term_sets = OX_sig0.1_GOslims)

# create a GO semantic similarity matrix from our GO slim IDs using get_sim_grid function in the OntologySimilarity package 
sim_matrix <- get_sim_grid(
  ontology=go, 
  information_content=GO_IC_data,
  term_sets=OX_sig0.1_GOslims)
#add column and row names to the semantic similarity matrix
rownames(sim_matrix) <- sig0.1_GOslims
colnames(sim_matrix) <- sig0.1_GOslims
#convert lower triangle of matrix to NA val including the diagonal
sim_matrix[lower.tri(sim_matrix, diag = TRUE)] <- NA
#reshape the data so that each GO ID combination and semantic similarity value is listed on a different row 
term_term <- melt(sim_matrix)
#make a new data frame with GO ID combinations with semantic similarity values greater than 0.5. This also excludes NAs.
term_term_0.5 <- term_term[which(term_term$value > 0.5),]

#####convert GO IDs to terms####

#first create an empty data frame the length of the term_term_0.5 data frame and with 3 columns
longterm_term <- data.frame(matrix(0,nrow(term_term_0.5),3))
#name the columns
colnames(longterm_term) <- c("Var1","Var2", "value")

#Loop through each row of the term_term_0.5 data frame
for(i in 1:nrow(term_term_0.5)){
  longterm_term$Var1[i] <- GSEA_BP_MF[which(GSEA_BP_MF$GOslim == term_term_0.5$Var1[i]),"Term"] #fill in new data frame with column 1 GO ID's GO term
  longterm_term$Var2[i] <- GSEA_BP_MF[which(GSEA_BP_MF$GOslim == term_term_0.5$Var2[i]),"Term"]#fill in new data frame with column 2 GO ID's GO term
  longterm_term$value[i] <- term_term_0.5$value[i]#fill in new data frame with the semantic similarity value for the GO combination
}
#create column with interaction type information for edge attribute file
longterm_term$type <- "term-term"
#rename column
colnames(longterm_term)[3] <- "semsimValue"

#create data frame containing protein-term data to merge with term-term data for edge attribute file
prot_term <- unique(goslim_protein[,c("gene","GOslim","GOslimTerm")])
prot_term$type <- "protein-term"
prot_term <- prot_term[,-2]
colnames(prot_term) <- c("Var1", "Var2", "type")
#merge term-term and protein-term data frames
edge_attb <- rbind.fill(longterm_term,prot_term)
#write out edge attribute file
write.table(edge_attb, "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSq_ASCA_clust/edge_attb_semsim0.5_sig0.1ASCAClust_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)

write.table(edge_attb[which(edge_attb$Var1 %in% GSEA_BP_MF[grep("BP",GSEA_BP_MF$GOcategory),"Term"] & edge_attb$Var2 %in% GSEA_BP_MF[grep("BP",GSEA_BP_MF$GOcategory),"Term"]),], "~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSq_ASCA_clust/edge_attb_GOBPsemsim0.5_sig0.1ASCAClust_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#####Calculate magnitude FC#############
#make list of unique GO Slims associated with these proteins
sig0.1_GOslims_terms_uniq <- unique(edge_attb[grep("protein", edge_attb$type),"Var2"])
days <- c(3,5,7,9,11,13)

go_term_magFC <- data.frame()
for (i in 1:length(sig0.1_GOslims_terms_uniq)){
  mag_df <- data.frame()
  pattern <- paste("^",sig0.1_GOslims_terms_uniq[i],"$", sep = "")
  protein_list <- prot_term[grep(pattern,prot_term$Var2),"Var1"]
  for(j in days){
    mag_TNS_sum <- sum(abs(all_sig0.1_pro_logFC_pval[which(all_sig0.1_pro_logFC_pval$gene %in% protein_list),paste("D_",j,"_logFC",sep = "")]))
    mag_NSAF_sum <- sum(abs(all_sig0.1_pro_logFC_pval[which(all_sig0.1_pro_logFC_pval$gene %in% protein_list),paste("D_",j,"_logFC_NSAF",sep = "")]))
    mag_df[1,c(paste("D",j,"_TNSmagFC", sep = ""),paste("D",j,"_NSAFmagFC", sep = ""))] <- c(mag_TNS_sum, mag_NSAF_sum)
  }
  mag_df$numprots <- length(protein_list)
  go_term_magFC <- rbind(go_term_magFC,mag_df)
}

go_term_magFC$term <- sig0.1_GOslims_terms_uniq
#find min and max magnitude fold change
min(unlist(go_term_magFC[,grep("magFC", colnames(go_term_magFC))]))
#[1] 0.005004154
max(unlist(go_term_magFC[,grep("magFC", colnames(go_term_magFC))]))
#[1] 57.97934; 118.5166

#there is quite a difference between max and min magnitude foldchanges

#find the spread of the number of proteins per term
summary(go_term_magFC$numprots)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    3.00    7.00   12.31   16.00   63.00 

#find the spread of the magnitude foldchanges
TNS_mags <- data.frame(unlist(go_term_magFC[,grep("TNSmagFC", colnames(go_term_magFC))]))
summary(TNS_mags[,1])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.005   1.074   2.617   4.940   6.207  57.979 
NSAF_mags <- data.frame(unlist(go_term_magFC[,grep("NSAFmagFC", colnames(go_term_magFC))]))
summary(NSAF_mags[,1])

#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.00525   1.36347   4.88123  12.58723  18.88961 118.51659 

#based on the above summary statistics, I will normalize the magnitude foldchanges by the number of proteins in each term

norm2prot <- data.frame()
for(i in 1:nrow(go_term_magFC)){
  row <- go_term_magFC[i,-grep("numprots|term", colnames(go_term_magFC))]/go_term_magFC$numprots[i]
  norm2prot <- rbind(norm2prot,row)
}
norm2prot$term <- go_term_magFC$term

go_term_magFC <- merge(go_term_magFC[,c("term", "numprots")], norm2prot, by = "term")
colnames(go_term_magFC) <- gsub("mag","norm_mag",colnames(go_term_magFC))
all_norm_mags <- data.frame(unlist(go_term_magFC[,grep("norm_mag", colnames(go_term_magFC))]))
summary(all_norm_mags[,1])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.005004 0.239779 0.343192 0.416699 0.488441 2.462816 
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.005004 0.270388 0.436059 0.682195 0.823171 9.619448 


###calulate sum and mean of TNS_normMagFC and NSAF_normMagFC

for(i in days){
  go_term_magFC[,paste0("D",i,"_TNS_NSAF_sum")] <- apply(go_term_magFC[,grep(paste("D",i,"_",sep = ""), colnames(go_term_magFC))], 1, sum)
  go_term_magFC[,paste0("D",i,"_TNS_NSAF_mean")] <- apply(go_term_magFC[,grep(paste("D",i,"_",sep = ""), colnames(go_term_magFC))], 1, mean)
}

all_norm_mag_means <- data.frame(unlist(go_term_magFC[,grep("mean", colnames(go_term_magFC))]))
summary(all_norm_mag_means[,1])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01846 0.44533 0.71696 0.90959 1.11938 6.43590 

write.table(go_term_magFC,"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSq_ASCA_clust/GOnode_attb_semsim0.5_sig0.1ASCAClust_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#only BP 
write.table(go_term_magFC[which(go_term_magFC$term %in% GSEA_BP_MF[grep("BP",GSEA_BP_MF$GOcategory),"Term"]),],"~/Documents/GitHub/OysterSeedProject/analysis/UniprotAnnotations_NetworkAnalysis/SameDayFCtoTemp/ChiSq_ASCA_clust/BPGOnode_attb_semsim0.5_sig0.1ASCAClust_EValcutoff.txt", sep = "\t", row.names = FALSE, quote = FALSE)



###looking at the distribution of norm mags

GO_TNSnormMag_dist <- data.frame(unlist(norm2prot[,grep("TNS", colnames(norm2prot))]))
GO_TNSnormMag_dist<- GO_TNSnormMag_dist[order(GO_TNSnormMag_dist$unlist.norm2prot...grep..TNS...colnames.norm2prot....),]
GO_NSAFnormMag_dist <- data.frame(unlist(norm2prot[,grep("NSAF", colnames(norm2prot))]))

plot(density(GO_TNSnormMag_dist$unlist.norm2prot...grep..TNS...colnames.norm2prot....))
