# Microbiome Proteomics Analysis Extension

## Overview

This document outlines the methodology to extend the existing oyster proteomics pipeline to analyze proteomic expression in non-target species, specifically the microbiome associated with _Crassostrea gigas_ during metamorphosis.

## Background

The original study focuses on oyster protein expression during metamorphosis using mass spectrometry data searched against the _C. gigas_ proteome. However, the same mass spectrometry data contains spectral information from microbial proteins present in the oyster samples. By expanding the search database and analysis pipeline, we can examine:

1. **Host-microbe interactions** during oyster development
2. **Microbial community changes** across temperature treatments
3. **Temporal dynamics** of microbiome protein expression
4. **Environmental influences** on microbiome composition and function

## Methodology Extension

### 1. Multi-Organism Database Construction

#### Database Components
- **Host proteins**: Existing _C. gigas_ proteome (`Cg-Giga_cont_AA.fa`)
- **Marine bacteria**: Common marine bacterial proteomes
- **Archaea**: Marine archaeal proteomes  
- **Fungi**: Marine fungal proteomes
- **Viruses**: Marine viral proteomes
- **Contaminants**: Existing contaminant database

#### Recommended Microbial Proteomes
Based on common marine microbiome studies, include proteomes from:

**Bacteria:**
- _Vibrio_ spp. (common marine pathogens)
- _Pseudoalteromonas_ spp. (marine bacteria)
- _Alteromonas_ spp. (marine bacteria)
- _Colwellia_ spp. (psychrophilic marine bacteria)
- _Shewanella_ spp. (marine bacteria)
- _Pseudomonas_ spp. (ubiquitous bacteria)
- _Bacillus_ spp. (spore-forming bacteria)

**Archaea:**
- _Nitrosopumilus_ spp. (marine ammonia-oxidizing archaea)
- _Thaumarchaeota_ representatives

**Fungi:**
- Marine yeasts and filamentous fungi

#### Database Construction Protocol

```bash
# Create combined database
cat Cg-Giga_cont_AA.fa > Combined_Host_Microbiome_DB.fa
cat marine_bacteria_proteomes/*.fa >> Combined_Host_Microbiome_DB.fa
cat marine_archaea_proteomes/*.fa >> Combined_Host_Microbiome_DB.fa
cat marine_fungi_proteomes/*.fa >> Combined_Host_Microbiome_DB.fa
cat contaminant_proteins.fa >> Combined_Host_Microbiome_DB.fa

# Add organism prefixes for identification
sed 's/>/>HOST_/' Cg-Giga_cont_AA.fa > temp_host.fa
sed 's/>/>BACTERIA_/' marine_bacteria.fa > temp_bacteria.fa
sed 's/>/>ARCHAEA_/' marine_archaea.fa > temp_archaea.fa
sed 's/>/>FUNGI_/' marine_fungi.fa > temp_fungi.fa

# Combine with prefixes
cat temp_host.fa temp_bacteria.fa temp_archaea.fa temp_fungi.fa contaminants.fa > Combined_Host_Microbiome_DB.fa
```

### 2. Modified Search Parameters

#### COMET Parameters
Use existing `Pcomet.params.high-low.txt` as base with modifications:

```
# Database path
database_name = Combined_Host_Microbiome_DB.fa

# Increased precursor mass tolerance for diverse organisms
peptide_mass_tolerance = 20.00

# Consider additional modifications common in microbes
variable_mod01 = 15.9949 M 0 3 -1 0 0
variable_mod02 = 42.0106 * 0 1 -1 0 0  # acetylation
variable_mod03 = 0.984016 N 0 1 -1 0 0  # deamidation

# Enzyme specificity - consider non-specific cleavage for microbes
search_enzyme_number = 0  # non-specific
allowed_missed_cleavage = 3
```

### 3. Identification and Quantification Strategy

#### Protein Identification Thresholds
- **Host proteins**: Existing thresholds (>95% probability)
- **Microbial proteins**: Conservative thresholds (>99% probability) to reduce false positives
- **Unique peptides**: Minimum 2 unique peptides per microbial protein

#### NSAF Calculation Modifications
Calculate separate NSAF values for:
- Host proteins only
- Microbial proteins only  
- Combined proteome

#### Quality Control Metrics
- Monitor host:microbe protein ratios
- Assess contamination levels
- Validate microbial identifications

### 4. Data Analysis Pipeline

#### Taxonomic Classification
```r
# Classify proteins by organism type
classify_proteins <- function(protein_ids) {
  host_proteins <- protein_ids[grepl("^HOST_", protein_ids)]
  bacterial_proteins <- protein_ids[grepl("^BACTERIA_", protein_ids)]
  archaeal_proteins <- protein_ids[grepl("^ARCHAEA_", protein_ids)]
  fungal_proteins <- protein_ids[grepl("^FUNGI_", protein_ids)]
  
  return(list(
    host = host_proteins,
    bacteria = bacterial_proteins,
    archaea = archaeal_proteins,
    fungi = fungal_proteins
  ))
}
```

#### Statistical Analysis Approaches

**1. Host-Microbe Ratio Analysis**
```r
# Calculate host:microbe protein abundance ratios
calculate_host_microbe_ratio <- function(nsaf_data) {
  host_abundance <- sum(nsaf_data[grepl("^HOST_", rownames(nsaf_data)), ])
  microbe_abundance <- sum(nsaf_data[!grepl("^HOST_", rownames(nsaf_data)), ])
  return(host_abundance / microbe_abundance)
}
```

**2. Separate PCA for Host and Microbiome**
```r
# Perform separate PCA analyses
host_pca <- prcomp(host_protein_data, center = TRUE, scale = TRUE)
microbe_pca <- prcomp(microbe_protein_data, center = TRUE, scale = TRUE)
```

**3. Combined Analysis with Organism Factor**
```r
# Include organism type as factor in ASCA
asca_combined <- ASCA(
  X = combined_protein_data,
  design = design_matrix,
  factors = c("time", "temperature", "organism_type")
)
```

### 5. Functional Analysis Extensions

#### Microbiome-Specific Annotations
- Use microbial-specific databases (KEGG, COG, Pfam)
- Focus on metabolic pathways relevant to marine environments
- Analyze antimicrobial resistance genes
- Examine biofilm formation pathways

#### Host-Microbe Interaction Analysis
- Identify co-expressed host-microbe protein modules
- Analyze pathway crosstalk between host and microbiome
- Examine temporal coordination of host-microbe responses

### 6. Visualization Strategies

#### Multi-Organism Plots
```r
# Stacked bar plots showing organism contributions
ggplot(organism_abundances, aes(x = sample, y = abundance, fill = organism)) +
  geom_bar(stat = "identity") +
  facet_wrap(~time_point)

# Correlation networks between host and microbial proteins
corrplot(host_microbe_correlations, method = "circle")
```

### 7. Quality Control and Validation

#### Controls and Thresholds
- Compare results with 16S rRNA sequencing if available
- Validate major microbial identifications with targeted approaches
- Monitor batch effects and contamination
- Use negative controls to assess background

#### Interpretation Guidelines
- Focus on consistent microbial signals across replicates
- Consider ecological relevance of identified microbes
- Validate novel findings with independent methods
- Account for potential host protein homology in microbes

## Expected Outcomes

1. **Temporal microbiome dynamics** during oyster metamorphosis
2. **Temperature effects** on microbiome composition and function
3. **Host-microbe interaction patterns** at the protein level
4. **Microbial biomarkers** associated with development stages
5. **Environmental adaptation** strategies of the holobiont

## Integration with Existing Analyses

This extension complements existing analyses by:
- Maintaining all current host-focused results
- Adding microbiome dimension to developmental studies
- Providing holobiont perspective on temperature responses
- Enabling systems-level understanding of oyster biology

## Computational Requirements

- **Increased search time**: 3-5x longer due to larger database
- **Memory requirements**: 2-3x increase for database indexing
- **Storage**: Additional space for microbial proteome databases
- **Analysis time**: Extended due to multi-organism processing

## Literature Integration

This approach enables comparison with:
- Marine microbiome studies
- Host-microbe interaction research
- Developmental microbiome dynamics
- Temperature effects on microbial communities

## Conclusion

Extending the proteomics pipeline to include microbiome analysis provides a comprehensive view of the oyster holobiont during metamorphosis, revealing previously hidden host-microbe interactions and environmental responses.