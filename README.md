# An integrated analysis reveals temperature-influenced proteomic variation throughout early development in the Pacific oyster



## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [Issues](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/issues)
- [Citation](#citation)

## Overview
The Pacific oyster has major ecological and economic importance serving as a biofilter and habitat in coastal ecosystems, and contributing over $190M to annual marine aquaculture revenue. However, little is known about the landscape of protein expression during early development, a time when mass mortality is common which can negatively impact industry and ecosystems. To better characterize physiological pathways and associated networks active during oyster development we performed a developmental time series proteomics analysis of larval cultures reared at 23°C and 29°C. These temperatures were selected based on the reports from the aquaculture industry that differential performance is observed in oysters at these temperatures. In addition to observing larger sized animals at 29°C, we found differentially abundant proteins that point to upregulated systemic structural remodeling during early development at 23°C and upregulated growth during late development at 29°C. The differential proteomes suggest that at 23°C more cellular energy is being diverted to maintenances processes where at 29°C more cellular energy is being used for growth. This proteomics analysis combined with development observations offers greater clarity on environmental conditions that can improve aquaculture production.

## Repo Contents

- ### [Analyses]():
	1. [**Proteomics\_Data\_Processing:**](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/Proteomics_Data_Processing)
		- [Abacus_parameters.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/Proteomics_Data_Processing/Abacus_parameters.txt):  ABACUS parameter file used in ABACUS analysis of mass spectrometry data
		- [DDA-data-Analyses.md](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/Proteomics_Data_Processing/DDA-data-Analyses.md):  DDA proteomics analysis workflow that converts files from .raw to .mzXML, searches .mzXML files against C. gigas database, calculates statistics associated with peptide and protein IDs using the Trans Proteomic Pipeline, and correlates protein inferences using ABACUS.  
		- [Pcomet.params.high-low.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/Proteomics_Data_Processing/Pcomet.params.high-low.txt):  COMET 2016.01 'comet.params.high-low' search parameters used in COMET analysis of mass spectrometry data download from http://comet-ms.sourceforge.net/parameters/parameters_201601/comet.params.high-low
		
	2. [**Technical\_Replicates\_PCA:**](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/Technical_Replicates_PCA)
		- [Technical\_Replicates\_PCA.Rproj](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/Technical_Replicates_PCA/Technical_Replicates_PCA.Rproj): R project used to calculate average NSAF values for proteins and run PCA to cluster technical replicates
		- [ClusteringTechnicalReplicates\_and_CalculateAvgNSAFs.R](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/Technical_Replicates_PCA/ClusteringTechnicalReplicates_and_CalculateAvgNSAFs.Rmd): R code used to calculate average NSAF values for proteins and run PCA to cluster technical replicates
	
	3. [**PCA\_ASCA\_Functional\_Analysis**](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/PCA_ASCA_Functional_Analysis):
		- [PCA\_ASCA\_Functional\_Analysis.Rproj](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/PCA_ASCA_Functional_Analysis/PCA_ASCA_Functional_Analysis.Rproj): R project used to run PCA, ASCA, and functional analysis on average protein NSAF values. 
		- [PCA\_ASCA\_Functional\_Analysis.Rmd](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/PCA_ASCA_Functional_Analysis/PCA_ASCA_Functional_Analysis.Rmd): R code used to perform PCA, ASCA, and functional analysis on average protein NSAF values. This produces Figures 2-6, Figures S2-3, and TableS3-6.
		- [Cg\_Giga\_cont\_AA.fa\_BLASTP\_uniprot\_swprot2019.ipynb](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/PCA_ASCA_Functional_Analysis/Cg_Giga_cont_AA.fa_BLASTP_uniprot_swprot2019.ipynb):  Jupyter notebook used to align protein sequences to the Uniprot-KB Swiss-Prot database using BLASTp 
	
	4. [Lit_Compare](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/Lit_Compare):
		- 
	
- ### [Data:](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Data) 
	1. [ABACUS\_output.tsv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/Abacus_output.tsv):  Output file from ABACUS analysis  
	2. [all_giga-uniprot-blastP-out.nopipe.annotations.tab](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/all_giga-uniprot-blastP-out.nopipe.annotations.tab): Proteins mapped to Uniprot database output from [Cg\_Giga\_cont\_AA.fa\_BLASTP\_uniprot\_swprot2019.ipynb](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/PCA_ASCA_Functional_Analysis/Cg_Giga_cont_AA.fa_BLASTP_uniprot_swprot2019.ipynb)
	3. [Average\_adjNSAF\_nozerovals.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/Average_adjNSAF_nozerovals.csv):  Average technical replicate NSAF values for each protein with zero values converted to 0.1 (1/8 of the lowest NSAF value). 
	4. [Cg\-Giga\_cont\_AA.fa](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/Cg-Giga_cont_AA.fa):  Protein sequence file ('Cg-Giga\_cont\_AA.fa') containing FASTA sequence file of the C. gigas proteome (downloaded from http://gigaton.sigenae.org as 'contigs.fasta.transdecoder.pep.gz') and common contaminants downloaded from the crapOME (Mellacheruvu et al. Nat Methods 2013). This file is used by proteomics spectral analysis and by UniProt mapping.
	5. [Sample_metadata.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/Sample_metadata.csv):  Sample ID and treatment information



- ### [Figures:](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Figures)
	1. [Figure 1](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Figure1.jpg): Size distribution based on sorting screen size of oysters at 24 days post fertilization when settlement rate was assessed
	2. [Figure 2](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Figure2.jpg): Global protein abundance in juvenile oyster over time under heat stress
	3. [Figure 3](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Figure3.jpg): Time influence on proteomes
	4. [Figure 4](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Figure4.pdf): Summary of biological processes represented by enriched GO terms within each clade for time-influenced proteins
	5. [Figure 5](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Figure5.jpg): Temperature influence on proteomes
	6. [Figure 6](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Figure6.pdf): Summary of biological processes represented by enriched GO terms within each clade for temperature-influenced proteins

- ### [Tables](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Tables)
	1. [Table 1: ASCA\_permTest\_summary.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Tables/ASCA_permTest_summary.txt):  Summarized results from ASCA permutation test
	2. [Table 2](): Proteins that commonly show increased abundance in response to high temperature.
	
- ### [Additional Files:](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/AdditionalFiles)
	1. [Additional File 1: BiovolumeAfterSettlement.xlsx](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/BiovolumeAfterSettlement.xlsx):  Biovolume of oyster seed (mL) settled after 6 days of exposure to high temperatures.
	2. [Additional File 2: Figure S1](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/FigureS1.jpg): PCA of all technical replicate samples
	3. [Additional File 3: PreliminaryProteomeCharacterization.xlsx](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/PreliminaryProteomeCharacterization.xlsx): Data and plot showing overlap between different temperature proteomes at each timepoint underlying Figure 2a
	4. [Additional File 4: PCA\_model\_summary.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/PCA_model_summary.txt):  PCA results summary
	5. [Additional File 5: Figure S2](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/FigureS2.jpg): Analysis of proteins influenced by temperature at 21 and 27 dpf qualitatively identified through PCA
	6. [Additional File 6: Table S1.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/TableS3.csv): 
	7. [Additional File 7: ASCA\_model\_summary.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/ASCA_model_summary.txt):  ASCA results summary
	8. [Additional File 8: Figure S3](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/FigureS3.jpg): ANOVA-simultaneous component analysis plots of PC loadings for all proteins
	9. [Additional File 9: Table S2.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/TableS2.csv):
	10. [Additional File 10: Table S3.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/AdditionalFiles/TableS3.csv):



# Citation
Trigg, S. A., Mitchell, K. M., Elliot, R., Euladiene, B., Vadopalas, B., Timmins-Schiffman, E. B., and Roberts, S. B. (2019) An integrated analysis reveals temperature-influenced proteomic variation throughout early development in the Pacific oyster. 
