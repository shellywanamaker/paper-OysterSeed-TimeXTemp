Literature Comparison
================
Shelly Trigg
3/23/2020

load libraries

``` r
library(readxl)
```

read in data

``` r
## Literature data


#Dineshram GBC 2016 table 6
temp_prots <- read_xlsx("downloaded_published_data/gcb13249-sup-0003-tables4.xlsx", skip = 2)


#pediveliger specific proteins
inSil_prots <- read.csv("downloaded_published_data/2019_Foulon_etal_inSilTRX_Cgigas.csv", header = TRUE, stringsAsFactors = FALSE)


## Read in all proteins detected in proteomics and BLAST'd to Uniprot DB **for C.gigas**. Need the Uniprot IDs from this table for matching the temperature- and time-affected proteins to.

# I did this blast on roadrunner here: /home/srlab/Documents/Shelly/Cgigas/Cgiga-uniprot-blastP-out.tab 
### ran the following commands:
##### makeblastdb -in UP000005408_29159.fasta -out UNIPROT_Cgigas_db -dbtype prot
##### blastp -query Cg_Giga_cont_AA.fa -db UNIPROT_Cgigas_db -out Cgiga-uniprot-blastP-out.tab -num_threads 8 -max_hsps 1 -max_target_seqs 1 -outfmt 6
### all input and output can be found here: https://gannet.fish.washington.edu/metacarcinus/Cgigas/BLAST2gigas_uniprot/

tempxtime_prots <- read.csv("Cgiga-uniprot-blastP-out.tab", sep = "\t", header = FALSE, stringsAsFactors = FALSE)


## Read in ASCA-identified temperature- and time-affected proteins
ASCA_Temp_clade <- read.csv("../../SupplementalFiles/Supplemental_Table_5.csv", stringsAsFactors = FALSE)

ASCA_Time_clade <- read.csv("../../SupplementalFiles/Supplemental_Table_4.csv", stringsAsFactors = FALSE)

PCA_d21d27_clade <- read.csv("../../SupplementalFiles/Supplemental_Table_3.csv", stringsAsFactors = FALSE)
```

format data

``` r
#change column name to match all dfs
colnames(tempxtime_prots)[1] <- "protein_ID"

colnames(ASCA_Temp_clade)[1] <- "protein_ID"
colnames(temp_prots)[3] <- "UniprotID"

#replace '|' with '.' and '-' with '.' in protein ID column of tempxtime_prots df
tempxtime_prots$protein_ID <- gsub("\\|","\\.", tempxtime_prots$protein_ID)
tempxtime_prots$protein_ID <- gsub("-","\\.", tempxtime_prots$protein_ID)
#replace 'tr|' with nothing in UniprotID column
tempxtime_prots$V2 <- gsub("tr\\|","", tempxtime_prots$V2)

#format column 2 to uniprot ID and entry name columns
tempxtime_prots$UniprotID <- gsub("\\|.*","", tempxtime_prots$V2)
tempxtime_prots$Entry_name <- gsub(".*\\|","", tempxtime_prots$V2)

#filter for proteins with low evalue
tempxtime_prots <- tempxtime_prots[which(tempxtime_prots$V11 <= 1 * 10^(-10)),]
```

merge Uniprot data with proteome data

``` r
#merge ASCA data with Uniprot data
ASCA_Temp_clade <- merge(ASCA_Temp_clade[which(ASCA_Temp_clade$ASCA_PC_threshold_pass=="Yes"),], tempxtime_prots[,c("protein_ID", "UniprotID", "Entry_name")], by = "protein_ID")

ASCA_Time_clade <- merge(ASCA_Time_clade[which(ASCA_Time_clade$ASCA_PC_threshold_pass=="Yes"),], tempxtime_prots[,c("protein_ID", "UniprotID", "Entry_name")], by = "protein_ID")

PCA_d21d27_clade <- merge(PCA_d21d27_clade[which(PCA_d21d27_clade$PC_threshold_pass=="Yes"),], tempxtime_prots[,c("protein_ID", "UniprotID", "Entry_name")], by = "protein_ID")
```

check overlap in ASCA temperature-affected
proteins

``` r
# compare Dineshram temp affected proteins with ASCA temp affected proteins
ASCA_temp_Dineshram_temp_common <- merge(ASCA_Temp_clade[which(ASCA_Temp_clade$ASCA_PC_threshold_pass=="Yes"),], temp_prots[,c(1,3:4,grep("Average:T", colnames(temp_prots)))], by = "UniprotID")

nrow(ASCA_temp_Dineshram_temp_common)
```

    ## [1] 42

``` r
#threshold for proteins that were more highly expressed at high temperature
ASCA_temp_Dineshram_temp_common_1.2x <- ASCA_temp_Dineshram_temp_common[ASCA_temp_Dineshram_temp_common[,16]>=1.2|ASCA_temp_Dineshram_temp_common[,17]>=1.2|ASCA_temp_Dineshram_temp_common[,18]>=1.2|ASCA_temp_Dineshram_temp_common[,19]>=1.2,]


#print the number of overlapping proteins that show increased abundance at 29C in our experiment
nrow(ASCA_temp_Dineshram_temp_common_1.2x[which(ASCA_temp_Dineshram_temp_common_1.2x$ClusterColor == "darkorange3"),])
```

    ## [1] 11

``` r
ASCA_temp_Dineshram_temp_common_1.2x[which(ASCA_temp_Dineshram_temp_common_1.2x$ClusterColor == "darkorange3"),c("protein_ID", "UniprotID", "Accession no.", "Protein_names")]
```

    ##                         protein_ID UniprotID Accession no.
    ## 4   CHOYP_LOC100367954.2.2.m.66596    K1P9U4  CGI_10005881
    ## 5            CHOYP_ADD.3.5.m.17639    K1PEX5  CGI_10006848
    ## 7   CHOYP_LOC101173335.4.4.m.49816    K1PFT9  CGI_10006016
    ## 10           CHOYP_RPS24.1.8.m.571    K1PUV4  CGI_10001493
    ## 11          CHOYP_NF70.1.4.m.31159    K1PWQ2  CGI_10018067
    ## 13     CHOYP_contig_043280.m.49983    K1PZS2  CGI_10005951
    ## 14     CHOYP_contig_044078.m.50900    K1Q086  CGI_10019530
    ## 21  CHOYP_LOC100705966.1.1.m.45957    K1QJR4  CGI_10019738
    ## 22 CHOYP_LOC100375029.6.10.m.36981    K1QM61  CGI_10009700
    ## 23 CHOYP_LOC100375029.8.10.m.60484    K1QM61  CGI_10009700
    ## 26  CHOYP_LOC100696604.1.1.m.40638    K1QP17  CGI_10010975
    ##                                                                                                                                           Protein_names
    ## 4                                                                                                                        C2 domain-containing protein 2
    ## 5                                                                                                         Protein hu-li tai shao (Adducin-like protein)
    ## 7                                                                                                                                             Myophilin
    ## 10                                                                                                                            40S ribosomal protein S24
    ## 11                                                                                                                  60 kDa neurofilament protein (NF60)
    ## 13                                                                                                                                                 <NA>
    ## 14                                                                                Ankyrin-2 (ANK-2) (Ankyrin-B) (Brain ankyrin) (Non-erythroid ankyrin)
    ## 21                                                                                                                                                 <NA>
    ## 22                                                                                                                   Thymosin beta (Tetrathymosin beta)
    ## 23                                                                                                                   Thymosin beta (Tetrathymosin beta)
    ## 26 Caprin-1 (Cytoplasmic activation- and proliferation-associated protein 1) (GPI-anchored protein p137) (GPI-p137) (p137GPI) (RNA granule protein 105)

check overlap in ASCA time-affected
proteins

``` r
ASCA_time_inSil <- merge(ASCA_Time_clade[which(ASCA_Time_clade$ASCA_PC_threshold_pass=="Yes"),], inSil_prots, by = "UniprotID")
nrow(ASCA_time_inSil)
```

    ## [1] 7

``` r
ASCA_time_inSil
```

    ##   UniprotID                          protein_ID ASCA_PC1_loadings
    ## 1    K1PX96         CHOYP_contig_026477.m.30166       -0.04863141
    ## 2    K1QBI0 CHOYP_BRAFLDRAFT_201924.1.1.m.30855       -0.04778439
    ## 3    K1QYU0 CHOYP_BRAFLDRAFT_129759.2.2.m.30082       -0.04436572
    ## 4    K1QZ56         CHOYP_contig_024582.m.27943       -0.04025016
    ## 5    K1R0T7             CHOYP_CO5A2.1.1.m.28030       -0.04077688
    ## 6    K1R8S4      CHOYP_LOC100485485.1.1.m.27961       -0.04385534
    ## 7    K1RKM7      CHOYP_LOC100372505.2.2.m.26614       -0.04052811
    ##   ASCA_PC2_loadings ASCA_PC_threshold_pass   ClusterColor SimpleClusterColor
    ## 1       -0.02363237                    Yes darkslategray1         light blue
    ## 2       -0.02204828                    Yes darkslategray1         light blue
    ## 3       -0.02043362                    Yes darkslategray1         light blue
    ## 4       -0.01862222                    Yes darkslategray1         light blue
    ## 5       -0.01786796                    Yes darkslategray1         light blue
    ## 6       -0.02790943                    Yes darkslategray1         light blue
    ## 7       -0.01810830                    Yes darkslategray1         light blue
    ##             Protein_names Protein_fams Gene_names                 GO_IDs
    ## 1                    <NA>         <NA>       <NA>                   <NA>
    ## 2   Perlucin-like protein                         GO:0005576; GO:0030246
    ## 3   Perlucin-like protein                         GO:0005576; GO:0030246
    ## 4                    <NA>         <NA>       <NA>                   <NA>
    ## 5 Collagen-like protein 7               MIMI_L669             GO:0019012
    ## 6   Perlucin-like protein                         GO:0005576; GO:0030246
    ## 7   Perlucin-like protein                         GO:0005576; GO:0030246
    ##   sig_bp_GO sig_bp_GOslim   Entry_name Ensembl.Gene.ID Protein.ID
    ## 1      <NA>          <NA> K1PX96_CRAGI    CGI_10004853   EKC21005
    ## 2      <NA>          <NA> K1QBI0_CRAGI    CGI_10010615   EKC18891
    ## 3      <NA>          <NA> K1QYU0_CRAGI    CGI_10006921   EKC42167
    ## 4      <NA>          <NA> K1QZ56_CRAGI    CGI_10026725   EKC38958
    ## 5      <NA>          <NA> K1R0T7_CRAGI    CGI_10010375   EKC27351
    ## 6      <NA>          <NA> K1R8S4_CRAGI    CGI_10006922   EKC42168
    ## 7      <NA>          <NA> K1RKM7_CRAGI    CGI_10006919   EKC42165
    ##                                Name Cell..Loc.
    ## 1 Hypothetical protein CGI_10004853   Ext 0.42
    ## 2             Aggrecan core protein   Ext 0.79
    ## 3             Perlucin-like protein   Ext 0.99
    ## 4 Hypothetical protein CGI_10026725   Ext 0.62
    ## 5           Collagen-like protein 7   Mem 0.72
    ## 6             Perlucin-like protein   Ext 0.99
    ## 7             Perlucin-like protein   Ext 0.95
    ##                                                       Group
    ## 1                                      Hypothetical protein
    ## 2 Calcification-related protein and calcium-binding protein
    ## 3 Calcification-related protein and calcium-binding protein
    ## 4                                      Hypothetical protein
    ## 5                                        Structural protein
    ## 6 Calcification-related protein and calcium-binding protein
    ## 7 Calcification-related protein and calcium-binding protein

``` r
ASCA_temp_inSil <- merge(ASCA_Temp_clade[which(ASCA_Temp_clade$ASCA_PC_threshold_pass=="Yes"),], inSil_prots, by = "UniprotID")
nrow(ASCA_temp_inSil)
```

    ## [1] 0

``` r
PCA_inSil <- merge(PCA_d21d27_clade[which(PCA_d21d27_clade$ASCA_PC_threshold_pass=="Yes"),], inSil_prots, by = "UniprotID")

nrow(PCA_inSil)
```

    ## [1] 0
