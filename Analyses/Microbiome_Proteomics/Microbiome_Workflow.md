# Microbiome Proteomics Workflow

## Overview

This workflow extends the existing oyster metamorphosis proteomics analysis to include examination of proteomic expression in non-target species, specifically the microbiome. This enables analysis of host-microbe interactions during _Crassostrea gigas_ development.

## Prerequisites

- Existing mass spectrometry data (.raw files)
- Original oyster proteome database (`Cg-Giga_cont_AA.fa`)
- Sample metadata (`Sample_metadata.csv`)
- Computational resources for extended database searches

## Workflow Steps

### Step 1: Database Preparation

#### 1.1 Obtain Microbial Proteomes

Download representative proteomes for common marine microorganisms:

```bash
# Example databases to download (adjust based on research focus)
# Marine bacteria
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Bacteria/

# Marine archaea  
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Archaea/

# Marine fungi
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/
```

#### 1.2 Create Combined Database

```bash
# Navigate to Data directory
cd ../../Data/

# Create organism-prefixed databases
sed 's/>/>HOST_/' Cg-Giga_cont_AA.fa > temp_host.fa

# Process microbial proteomes (example for bacteria)
for file in marine_bacteria/*.faa; do
    organism=$(basename "$file" .faa)
    sed "s/>/>BACTERIA_${organism}_/" "$file" >> temp_bacteria.fa
done

# Combine all databases
cat temp_host.fa temp_bacteria.fa temp_archaea.fa temp_fungi.fa > Combined_Host_Microbiome_DB.fa

# Clean up temporary files
rm temp_*.fa

# Index database for COMET
# (Database indexing will be handled automatically by COMET)
```

### Step 2: Mass Spectrometry Search

#### 2.1 Convert Raw Files (if not already done)

```bash
# Use existing workflow from Proteomics_Data_Processing/DDA-data-Analyses.md
cd ../Analyses/Microbiome_Proteomics/

# Convert .raw to .mzXML files using existing protocol
for file in ../../Data/*.raw; do
    no_path=${file##*/}
    no_ext=${no_path%.raw}
    WINEPREFIX=~/.wine32 ReAdW.2016010.msfilereader.exe "$file" "$no_ext".mzXML
done
```

#### 2.2 Database Search with COMET

```bash
# Run COMET search with microbiome parameters
comet.exe -P Comet_params_microbiome.txt *.mzXML

# This will generate .pep.xml files for each sample
```

### Step 3: Protein Inference and Quantification

#### 3.1 Process with Trans-Proteomic Pipeline

```bash
# Convert to pepXML format and calculate probabilities
for file in *.pep.xml; do
    # Peptide validation
    xinteract -N"$file".interact.pep.xml -p0.05 -l7 "$file"
    
    # Protein inference
    ProteinProphet "${file}.interact.pep.xml" "${file}.interact.prot.xml"
done
```

#### 3.2 Quantification with ABACUS

```bash
# Run ABACUS with microbiome-specific parameters
java -jar abacus.jar -p Abacus_params_microbiome.txt *.prot.xml

# This generates combined quantification results
```

### Step 4: Data Analysis

#### 4.1 Load Analysis Framework

```r
# In R, load the microbiome analysis functions
source("Microbiome_Analysis_Example.R")

# Load required libraries (automatically loaded by the script)
```

#### 4.2 Process Results

```r
# Load ABACUS output
microbiome_data <- read.csv("./abacus_output/microbiome_Abacus_output.tsv", 
                           sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# Load metadata
metadata <- read.csv("../../Data/Sample_metadata.csv", stringsAsFactors = FALSE)

# Run microbiome analysis
results <- run_microbiome_analysis(
  combined_nsaf_file = "./abacus_output/microbiome_Abacus_output.tsv",
  metadata_file = "../../Data/Sample_metadata.csv",
  output_dir = "./microbiome_results"
)
```

#### 4.3 Organism Classification and Analysis

```r
# Classify proteins by organism type
protein_classification <- classify_proteins(rownames(microbiome_data))

# Separate data by organism
organism_data <- separate_organisms(microbiome_data)

# Calculate diversity metrics
diversity_results <- list()
for (org in names(organism_data)) {
  if (nrow(organism_data[[org]]) > 5) {
    diversity_results[[org]] <- calculate_organism_diversity(organism_data[[org]])
  }
}

# Calculate host-microbiome ratios
host_microbe_ratios <- calculate_host_microbe_ratios(microbiome_data)
metadata$host_microbe_ratio <- host_microbe_ratios[metadata$SampleID]
```

### Step 5: Statistical Analyses

#### 5.1 PCA Analysis by Organism

```r
# Perform PCA for each organism type
pca_results <- list()
for (org in names(organism_data)) {
  if (nrow(organism_data[[org]]) > 10) {
    log_data <- log2(organism_data[[org]] + 0.1)
    pca_results[[org]] <- prcomp(t(log_data), center = TRUE, scale = TRUE)
  }
}
```

#### 5.2 Extended ASCA Analysis

```r
# Extend existing ASCA to include organism type
# Load existing ASCA code and modify to include organism factor

# Create design matrix with organism information
design_extended <- metadata
design_extended$organism_dominance <- ifelse(
  design_extended$host_microbe_ratio > median(design_extended$host_microbe_ratio, na.rm = TRUE),
  "host_dominant", "microbe_dominant"
)

# Run ASCA with additional factors
# (This extends the existing ASCA approach from PCA_ASCA_Functional_Analysis.Rmd)
```

#### 5.3 Host-Microbiome Correlations

```r
# Analyze protein correlations between host and microbiome
if ("host" %in% names(organism_data) && "bacteria" %in% names(organism_data)) {
  correlation_results <- analyze_host_microbe_correlations(
    organism_data$host, 
    organism_data$bacteria
  )
  
  # Save significant correlations
  write.csv(correlation_results$significant_correlations,
            "./microbiome_results/host_microbiome_correlations.csv")
}
```

### Step 6: Visualization and Interpretation

#### 6.1 Create Microbiome-Specific Plots

```r
# Organism composition over time
organism_summary <- data.frame(
  sample = rep(colnames(microbiome_data), length(names(organism_data))),
  organism = rep(names(organism_data), each = ncol(microbiome_data)),
  abundance = as.vector(sapply(organism_data, colSums))
)

ggplot(organism_summary, aes(x = sample, y = abundance, fill = organism)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~time_point) +
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(title = "Organism Composition Across Time Points",
       x = "Sample", y = "Relative Abundance")

# Host-microbiome ratio over development
ggplot(metadata, aes(x = dpf, y = host_microbe_ratio, color = as.factor(temp))) +
  geom_point(size = 3) +
  geom_smooth(method = "loess") +
  scale_color_manual(values = c("cyan3", "magenta3")) +
  theme_minimal() +
  labs(title = "Host-Microbiome Balance During Development",
       x = "Days Post Fertilization", y = "Host:Microbiome Protein Ratio")
```

#### 6.2 Functional Analysis

```r
# Perform organism-specific GO/KEGG enrichment
# Use appropriate databases for each organism type

# Identify metabolic pathway interactions between host and microbiome
# Cross-reference with existing functional analysis from original study
```

### Step 7: Integration with Existing Results

#### 7.1 Validate Host Protein Consistency

```r
# Load original host-only results
original_results <- read.csv("../PCA_ASCA_Functional_Analysis/data_file.csv")

# Compare host protein quantification between analyses
# Ensure microbiome extension doesn't compromise host analysis quality

# Calculate correlation between host-only and combined analysis for host proteins
host_only_proteins <- rownames(original_results)
host_combined_proteins <- organism_data$host

# Validate consistency (should show high correlation)
consistency_check <- cor(host_only_data, host_combined_data, use = "complete.obs")
```

#### 7.2 Extended Interpretation

```r
# Add microbiome context to existing developmental findings
# Examine how microbiome changes relate to:
# - Metamorphic transitions
# - Temperature responses  
# - Physiological processes identified in original study

# Create integrated visualizations combining host and microbiome data
```

## Expected Outputs

### 1. Quantitative Results
- `organism_classification_summary.csv`: Protein counts by organism
- `*_diversity_metrics.csv`: Diversity metrics for each organism type
- `host_microbiome_correlations.csv`: Significant host-microbe protein correlations
- `combined_NSAF_results.tsv`: Full quantification matrix with organism annotations

### 2. Statistical Analysis Results
- Extended PCA plots showing organism-specific clustering
- ASCA results including organism/microbiome factors
- Correlation matrices and networks

### 3. Visualizations
- Organism composition stacked bar plots over time
- Host-microbiome ratio trends across development
- Correlation heatmaps between host and microbial proteins
- Integration plots with existing oyster development data

### 4. Functional Analysis
- Organism-specific pathway enrichment results
- Host-microbiome pathway interaction networks
- Temporal dynamics of microbial functional categories

## Quality Control and Validation

### 1. Technical Validation
- Compare host protein results with original analysis (should show >0.9 correlation)
- Check organism classification accuracy using known markers
- Validate microbial identifications with literature/databases

### 2. Biological Validation
- Confirm ecological relevance of identified microorganisms
- Cross-reference with 16S rRNA data if available
- Validate functional predictions with known host-microbe interactions

### 3. Statistical Validation
- Use appropriate multiple testing corrections
- Apply conservative thresholds for microbial protein identification
- Validate temporal patterns across replicates

## Integration with Original Study

This microbiome extension:
- **Preserves** all existing host-focused analyses and results
- **Adds** microbiome dimension to developmental proteomics
- **Enables** holobiont-level understanding of oyster metamorphosis
- **Provides** new insights into environmental adaptation mechanisms
- **Maintains** compatibility with existing statistical and visualization approaches

## Computational Requirements

- **Runtime**: 3-5x increase due to larger database
- **Memory**: 2-3x increase for database handling
- **Storage**: Additional 50-100GB for microbial databases
- **Processing**: Extended analysis time for multi-organism processing

## Literature Integration

Results can be compared with:
- Marine microbiome studies in molluscs
- Host-microbiome interactions during development
- Temperature effects on microbial communities
- Proteomics studies of marine microbial communities

## Conclusion

This workflow provides a comprehensive approach to extend existing oyster proteomics analysis to include microbiome examination, enabling investigation of host-microbe interactions during metamorphosis while preserving all original analytical capabilities.