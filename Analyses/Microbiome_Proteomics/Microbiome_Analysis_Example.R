# Microbiome Proteomics Analysis Script
# Extension of the existing oyster proteomics pipeline to include microbiome analysis
# Author: Computational Analysis Extension
# Date: 2024

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(corrplot)
library(pheatmap)
library(vegan)
library(MetStaT)

# ==============================================================================
# FUNCTIONS FOR MICROBIOME PROTEOMICS ANALYSIS
# ==============================================================================

#' Classify proteins by organism type based on prefixes
#' @param protein_ids Vector of protein IDs with organism prefixes
#' @return List with proteins classified by organism type
classify_proteins <- function(protein_ids) {
  classification <- list(
    host = protein_ids[grepl("^HOST_", protein_ids)],
    bacteria = protein_ids[grepl("^BACTERIA_", protein_ids)],
    archaea = protein_ids[grepl("^ARCHAEA_", protein_ids)],
    fungi = protein_ids[grepl("^FUNGI_", protein_ids)],
    contaminants = protein_ids[grepl("^CONTAM_", protein_ids)]
  )
  
  # Add summary statistics
  classification$summary <- data.frame(
    organism = names(classification)[1:5],
    count = sapply(classification[1:5], length),
    stringsAsFactors = FALSE
  )
  
  return(classification)
}

#' Calculate host-to-microbiome abundance ratios
#' @param nsaf_data Matrix of NSAF values with proteins as rows, samples as columns
#' @return Vector of host:microbe ratios per sample
calculate_host_microbe_ratios <- function(nsaf_data) {
  host_proteins <- rownames(nsaf_data)[grepl("^HOST_", rownames(nsaf_data))]
  microbe_proteins <- rownames(nsaf_data)[!grepl("^HOST_|^CONTAM_", rownames(nsaf_data))]
  
  host_abundance <- colSums(nsaf_data[host_proteins, , drop = FALSE])
  microbe_abundance <- colSums(nsaf_data[microbe_proteins, , drop = FALSE])
  
  # Avoid division by zero
  microbe_abundance[microbe_abundance == 0] <- min(microbe_abundance[microbe_abundance > 0]) * 0.1
  
  ratios <- host_abundance / microbe_abundance
  return(ratios)
}

#' Separate NSAF data by organism type
#' @param nsaf_data Matrix of NSAF values
#' @return List of NSAF matrices separated by organism
separate_organisms <- function(nsaf_data) {
  protein_classification <- classify_proteins(rownames(nsaf_data))
  
  organism_data <- list()
  for (org in names(protein_classification)[1:4]) {  # Exclude contaminants and summary
    if (length(protein_classification[[org]]) > 0) {
      organism_data[[org]] <- nsaf_data[protein_classification[[org]], , drop = FALSE]
    }
  }
  
  return(organism_data)
}

#' Calculate organism-specific diversity metrics
#' @param nsaf_data Matrix of NSAF values for one organism type
#' @return Data frame with diversity metrics per sample
calculate_organism_diversity <- function(nsaf_data) {
  # Convert to relative abundance (proportion)
  rel_abundance <- apply(nsaf_data, 2, function(x) x / sum(x))
  
  # Calculate Shannon diversity
  shannon <- apply(rel_abundance, 2, function(x) {
    x <- x[x > 0]  # Remove zeros
    -sum(x * log(x))
  })
  
  # Calculate Simpson diversity (1 - Simpson index)
  simpson <- apply(rel_abundance, 2, function(x) {
    x <- x[x > 0]
    1 - sum(x^2)
  })
  
  # Calculate richness (number of detected proteins)
  richness <- apply(nsaf_data > 0, 2, sum)
  
  diversity_metrics <- data.frame(
    sample = colnames(nsaf_data),
    shannon = shannon,
    simpson = simpson,
    richness = richness,
    stringsAsFactors = FALSE
  )
  
  return(diversity_metrics)
}

#' Perform correlation analysis between host and microbial proteins
#' @param host_data Matrix of host protein NSAF values
#' @param microbe_data Matrix of microbial protein NSAF values
#' @return Correlation matrix and significant correlations
analyze_host_microbe_correlations <- function(host_data, microbe_data) {
  # Ensure samples match
  common_samples <- intersect(colnames(host_data), colnames(microbe_data))
  host_subset <- host_data[, common_samples, drop = FALSE]
  microbe_subset <- microbe_data[, common_samples, drop = FALSE]
  
  # Calculate correlations
  correlations <- cor(t(host_subset), t(microbe_subset), use = "complete.obs")
  
  # Find significant correlations (example threshold)
  correlation_threshold <- 0.7
  significant_pairs <- which(abs(correlations) > correlation_threshold, arr.ind = TRUE)
  
  if (nrow(significant_pairs) > 0) {
    significant_correlations <- data.frame(
      host_protein = rownames(correlations)[significant_pairs[, 1]],
      microbe_protein = colnames(correlations)[significant_pairs[, 2]],
      correlation = correlations[significant_pairs],
      stringsAsFactors = FALSE
    )
  } else {
    significant_correlations <- data.frame(
      host_protein = character(0),
      microbe_protein = character(0),
      correlation = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    correlation_matrix = correlations,
    significant_correlations = significant_correlations
  ))
}

# ==============================================================================
# EXAMPLE ANALYSIS WORKFLOW
# ==============================================================================

#' Main analysis function demonstrating microbiome proteomics workflow
#' This function shows how to integrate microbiome analysis with existing pipeline
run_microbiome_analysis <- function(combined_nsaf_file, metadata_file, output_dir = ".") {
  
  cat("Starting microbiome proteomics analysis...\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ===========================================================================
  # 1. DATA LOADING AND PREPROCESSING
  # ===========================================================================
  
  cat("1. Loading data...\n")
  
  # Load combined NSAF data (host + microbiome)
  # This would be the output from ABACUS analysis with combined database
  # For demonstration, we'll show the expected structure
  
  # Example structure (replace with actual data loading):
  # combined_data <- read.csv(combined_nsaf_file, row.names = 1, stringsAsFactors = FALSE)
  # metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  # For demonstration purposes, create example data structure
  cat("   Note: This is a demonstration of expected data structure\n")
  cat("   Replace with actual combined ABACUS output loading\n")
  
  # ===========================================================================
  # 2. ORGANISM CLASSIFICATION
  # ===========================================================================
  
  cat("2. Classifying proteins by organism type...\n")
  
  # Example of how to classify proteins (use with real data)
  # protein_classification <- classify_proteins(rownames(combined_data))
  # 
  # cat("   Protein counts by organism:\n")
  # print(protein_classification$summary)
  
  # ===========================================================================
  # 3. ORGANISM-SPECIFIC ANALYSES
  # ===========================================================================
  
  cat("3. Performing organism-specific analyses...\n")
  
  # Separate data by organism
  # organism_data <- separate_organisms(combined_data)
  
  # Calculate diversity metrics for each organism type
  # diversity_results <- list()
  # for (org in names(organism_data)) {
  #   if (nrow(organism_data[[org]]) > 5) {  # Only if sufficient proteins
  #     diversity_results[[org]] <- calculate_organism_diversity(organism_data[[org]])
  #   }
  # }
  
  # ===========================================================================
  # 4. HOST-MICROBIOME RATIO ANALYSIS
  # ===========================================================================
  
  cat("4. Analyzing host-microbiome abundance ratios...\n")
  
  # Calculate ratios
  # host_microbe_ratios <- calculate_host_microbe_ratios(combined_data)
  # 
  # # Add ratios to metadata
  # metadata$host_microbe_ratio <- host_microbe_ratios[metadata$SampleID]
  
  # ===========================================================================
  # 5. STATISTICAL ANALYSES
  # ===========================================================================
  
  cat("5. Performing statistical analyses...\n")
  
  # PCA analysis for each organism type
  # pca_results <- list()
  # for (org in names(organism_data)) {
  #   if (nrow(organism_data[[org]]) > 10) {
  #     log_data <- log2(organism_data[[org]] + 0.1)
  #     pca_results[[org]] <- prcomp(t(log_data), center = TRUE, scale = TRUE)
  #   }
  # }
  
  # ASCA analysis including organism type as factor
  # This would extend the existing ASCA approach
  
  # ===========================================================================
  # 6. CORRELATION ANALYSIS
  # ===========================================================================
  
  cat("6. Analyzing host-microbiome correlations...\n")
  
  # Host-microbiome correlations
  # if ("host" %in% names(organism_data) && "bacteria" %in% names(organism_data)) {
  #   correlation_results <- analyze_host_microbe_correlations(
  #     organism_data$host, 
  #     organism_data$bacteria
  #   )
  # }
  
  # ===========================================================================
  # 7. VISUALIZATION
  # ===========================================================================
  
  cat("7. Creating visualizations...\n")
  
  # Example visualization code (use with real data)
  create_microbiome_visualizations <- function() {
    
    # Organism composition plot
    # ggplot(organism_summary, aes(x = sample, y = abundance, fill = organism)) +
    #   geom_bar(stat = "identity", position = "stack") +
    #   facet_wrap(~time_point) +
    #   scale_fill_viridis_d() +
    #   theme_minimal() +
    #   labs(title = "Organism Composition Across Time Points",
    #        x = "Sample", y = "Relative Abundance")
    
    # Host-microbiome ratio plot
    # ggplot(metadata, aes(x = dpf, y = host_microbe_ratio, color = as.factor(temp))) +
    #   geom_point(size = 3) +
    #   geom_smooth(method = "loess") +
    #   scale_color_manual(values = c("cyan3", "magenta3")) +
    #   theme_minimal() +
    #   labs(title = "Host-Microbiome Ratio Over Time",
    #        x = "Days Post Fertilization", y = "Host:Microbiome Protein Ratio")
    
    # Correlation heatmap
    # if (exists("correlation_results")) {
    #   pheatmap(correlation_results$correlation_matrix[1:20, 1:20],  # Subset for visibility
    #            cluster_rows = TRUE, cluster_cols = TRUE,
    #            color = colorRampPalette(c("blue", "white", "red"))(100),
    #            main = "Host-Microbiome Protein Correlations")
    # }
  }
  
  # ===========================================================================
  # 8. FUNCTIONAL ANALYSIS
  # ===========================================================================
  
  cat("8. Performing functional analysis...\n")
  
  # Organism-specific GO/KEGG enrichment
  # This would use organism-specific annotation databases
  
  # Pathway interaction analysis between host and microbiome
  # Cross-reference metabolic pathways
  
  # ===========================================================================
  # 9. INTEGRATION WITH EXISTING ANALYSES
  # ===========================================================================
  
  cat("9. Integrating with existing host-focused analyses...\n")
  
  # Compare host protein results with original analysis
  # Validate consistency of host protein quantification
  
  # Add microbiome dimension to existing PCA/ASCA results
  # Extend functional enrichment with microbiome pathways
  
  # ===========================================================================
  # 10. OUTPUT RESULTS
  # ===========================================================================
  
  cat("10. Saving results...\n")
  
  # Save organism classification results
  # write.csv(protein_classification$summary, 
  #           file.path(output_dir, "organism_classification_summary.csv"))
  
  # Save diversity metrics
  # for (org in names(diversity_results)) {
  #   write.csv(diversity_results[[org]], 
  #             file.path(output_dir, paste0(org, "_diversity_metrics.csv")))
  # }
  
  # Save correlation results
  # if (exists("correlation_results")) {
  #   write.csv(correlation_results$significant_correlations,
  #             file.path(output_dir, "host_microbiome_correlations.csv"))
  # }
  
  cat("Analysis complete! Check output directory for results.\n")
  
  return(list(
    message = "Microbiome analysis framework ready",
    note = "Replace example code with actual data processing",
    output_directory = output_dir
  ))
}

# ==============================================================================
# INTEGRATION FUNCTIONS
# ==============================================================================

#' Integrate microbiome results with existing oyster analysis
#' @param original_results Path to original PCA/ASCA results
#' @param microbiome_results Microbiome analysis results
#' @return Integrated analysis results
integrate_with_existing_analysis <- function(original_results, microbiome_results) {
  
  # Load original host-only results
  # Compare host protein consistency between analyses
  # Add microbiome dimensions to existing plots
  # Create combined statistical models
  
  cat("Integration functions ready for implementation\n")
  cat("This would combine host-only and host+microbiome analyses\n")
  
  return(list(
    status = "Integration framework prepared",
    notes = "Implement based on specific research questions"
  ))
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

if (FALSE) {  # Set to TRUE to run example (with real data)
  
  # Example usage of the microbiome analysis pipeline
  results <- run_microbiome_analysis(
    combined_nsaf_file = "../../Data/Combined_Microbiome_NSAF.csv",
    metadata_file = "../../Data/Sample_metadata.csv",
    output_dir = "./microbiome_results"
  )
  
  print(results)
}

# Print framework ready message
cat("\n================================================================================\n")
cat("MICROBIOME PROTEOMICS ANALYSIS FRAMEWORK LOADED\n")
cat("================================================================================\n\n")
cat("This script provides the analytical framework for extending\n")
cat("the existing oyster proteomics pipeline to include microbiome analysis.\n\n")
cat("Key functions available:\n")
cat("- classify_proteins(): Separate proteins by organism type\n")
cat("- calculate_host_microbe_ratios(): Quantify host:microbiome balance\n") 
cat("- separate_organisms(): Split NSAF data by organism\n")
cat("- calculate_organism_diversity(): Compute diversity metrics\n")
cat("- analyze_host_microbe_correlations(): Find protein correlations\n")
cat("- run_microbiome_analysis(): Main analysis pipeline\n\n")
cat("To use with real data:\n")
cat("1. Generate combined host-microbiome database\n")
cat("2. Run COMET + ABACUS with microbiome parameters\n")
cat("3. Load combined NSAF results into this framework\n")
cat("4. Execute analysis functions with real data\n\n")
cat("================================================================================\n")