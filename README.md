# An integrated analysis reveals temperature-influenced proteomic variation throughout early development in the Pacific oyster



## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [Issues](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/issues)
- [Citation](#citation)

# Overview
The Pacific oyster has major ecological and economic importance serving as a biofilter and habitat in coastal ecosystems, and contributing over $190M to annual marine aquaculture revenue. However, little is known about the landscape of protein expression during early development, a time when mass mortality is common which can negatively impact industry and ecosystems. To better characterize physiological pathways and associated networks active during oyster development we performed a developmental time series proteomics analysis of larval cultures reared at 23°C and 29°C. These temperatures were selected based on the reports from the aquaculture industry that differential performance is observed in oysters at these temperatures. In addition to observing larger sized animals at 29°C, we found differentially abundant proteins that point to upregulated systemic structural remodeling during early development at 23°C and upregulated growth during late development at 29°C. The differential proteomes suggest that at 23°C more cellular energy is being diverted to maintenances processes where at 29°C more cellular energy is being used for growth. This proteomics analysis combined with development observations offers greater clarity on environmental conditions that can improve aquaculture production.

# Repo Contents
- **Manuscript files:**
	- Additional files:
		1. Additional file 1:
		2. Additional file 2:
		3. Additional file 3:
		4. Additional file 4:
		5. Additional file 5: 
		6. Additional file 6:
		7. Additional file 7:
		8. Additional file 8:

- **[Analyses]()**:
	1. [Proteomics\_Data\_Processing:](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/Proteomics_Data_Processing)
		- [Abacus_parameters.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Scripts/Abacus_parameters.txt):  ABACUS parameter file used in ABACUS analysis of mass spectrometry data
		- [DDA-data-Analyses.md](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Scripts/DDA-data-Analyses.md):  DDA proteomics analysis workflow that converts files from .raw to .mzXML, searches .mzXML files against C. gigas database, calculates statistics associated with peptide and protein IDs using the Trans Proteomic Pipeline, and correlates protein inferences using ABACUS.  
		- [Pcomet.params.high-low.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Scripts/Pcomet.params.high-low.txt):  COMET 2016.01 'comet.params.high-low' search parameters used in COMET analysis of mass spectrometry data download from http://comet-ms.sourceforge.net/parameters/parameters_201601/comet.params.high-low
		
	2. [Technical\_Replicates\_PCA:](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/Technical_Replicates_PCA)
		- [Technical\_Replicates\_PCA.Rproj](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/Technical_Replicates_PCA/Technical_Replicates_PCA.Rproj): R project used to run PCA to cluster technical replicates and calculate average NSAF values
		- [ClusteringTechnicalReplicates\_and_CalculateAvgNSAFs.R](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/Technical_Replicates_PCA/ClusteringTechnicalReplicates_and_CalculateAvgNSAFs.R): R code used to cluster technical replicates by PCA and calulate average NSAF values
		- [silo3and9\_nozerovals\_NSAF\_AVGs.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/Technical_Replicates_PCA/silo3and9_nozerovals_NSAF_AVGs.csv): Average technical replicate NSAF values for each protein 
	
	3. [ASCA](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Analyses/ASCA):
		- [ASCA.Rproj](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA.Rproj):
		- [ASCA\_silo3\_9\_log2avgNSAF.Rmd](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_silo3_9_log2avgNSAF.Rmd):
		- [ASCA\_model\_summary.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_model_summary.txt):
		- [ASCA\_permTest\_summary.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_permTest_summary.txt):
		- [ASCA\_time\_PC1andPC2\_loadings.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_time_PC1andPC2_loadings.csv):
		- [ASCA\_TimeProts\_avgNSAF.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_TimeProts_avgNSAF.csv):
		- [ASCA\_TimeProts\_autoscaled\_avgNSAF.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_TimeProts_autoscaled_avgNSAF.csv):
		- [ASCA\_temperature\_PC1\_loadings.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_temperature_PC1_loadings.csv):
		- [ASCA\_TempProts\_avgNSAF.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_TempProts_avgNSAF.csv):
		- [ASCA\_TempProts\_autoscaled\_avgNSAF.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Analyses/ASCA/ASCA_TempProts_autoscaled_avgNSAF.csv):

	4. [Functional_Analysis]():
		- []():
	
	
	
	
	
	 
- **[R](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/R)**: R code
	1. 	[ClusteringTechnicalReplicates\_and\_CalculateAvgNSAFs.R](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/R/ClusteringTechnicalReplicates_and_CalculateAvgNSAFs.R):   

	2. [General\_proteome\_characterization.Rmd](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/R/General_proteome_characterization.Rmd):  R code used to  convert the GO identifiers associated with each protein to GO slim terms using the MGI GO Slim database (http://www.informatics.jax.org/gotools/MGI_GO_Slim.html), and plot pie chart of biological processes.  
	3. [clustering-dendrogram\_s3\_9.R](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/R/clustering-dendrogram_s3_9.R):  R code used for hierarchical clustering do identify differentially abundant proteins
	4. [ASCA.R](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/R/ASCA.R):  R code used to run ASCA and identify differentially abundant proteins
	5. [LogFC\_ChiSquaredTest\_SameDaySamples.Rmd](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/R/LogFC_ChiSquaredTest_SameDaySamples.Rmd): R code used for calculating log foldchange and Chi squared proportions test on all protein 
	6. [VerifyStatsProteinSelection.Rmd](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/R/VerifyStatsProteinSelection.Rmd): R code that selects proteins with < 0.1 adjusted Chi Squared p-value, makes venn diagram (Fig. 2) plot to show the overlap between all statistical methods, plots protein abundance for each statistical method as facetted line plots and as heatmaps (Supplementary Figs. 2-4)
	7. [CreateNodeAndGoSemSimEdgeAttr.R](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/R/CreateNodeAndGoSemSimEdgeAttr.R): R code that retrieves GO slim identifiers for temperature-affected proteins using GSEA, calculates and thresholds GO semantic similarity scores for GO slims, calculates the total magnitude foldchange NSAF and total number of spectra for each term based on all proteins associated with each term, generates node and edge attribute files for Cytoscape.


- **[Data](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Data)**: 
	1. [BiovolumeAfterSettlement.xlsx](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/BiovolumeAfterSettlement.xlsx):  Biovolume of oyster seed (mL) settled after 6 days of exposure to high temperatures.
	2. [Cg\-Giga\_cont\_AA.fa](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/Cg-Giga_cont_AA.fa):  Protein sequence database file ('Cg-Giga\_cont\_AA.fa') containing FASTA sequence file of the C. gigas proteome (downloaded from http://gigaton.sigenae.org as 'contigs.fasta.transdecoder.pep.gz') and common contaminants downloaded from the crapOME (Mellacheruvu et al. Nat Methods 2013). This file is used by proteomics spectral analysis and by UniProt mapping.
	3. [ABACUS\_output.tsv.gz](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/Abacus_output.tsv):  Output file from ABACUS analysis  
	4. [Sample_metadata.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/Sample_metadata.csv):  Sample ID and treatment information
	5. [silo3and9\_NSAF\_AVGs.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/silo3and9_NSAF_AVGs.csv):  Average NSAF values for all proteins (listed as rows) detected in all samples (listed as columns)
	6. [silo3and9\_nozerovals\_NSAF\_AVGs.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/silo3and9_nozerovals_NSAF_AVGs.csv):  Average NSAF values for all proteins (listed as columns) detected in all samples (listed as rows) with zero value converted to 0.1. 
	7. [PreliminaryProteomeCharacterization.xlsx](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/PreliminaryProteomeCharacterization.xlsx): Data and plot showing overlap between different temperature proteomes at each timepoint underlying Figure 2a.  
	8. [HC_UniquelyClusteredProteins.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/HC_UniquelyClusteredProteins.csv):  Proteins identified by hierarchical clustering as showing differential abundance between temperature treatments. 
	9. [ASCA\_model\_summary.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/ASCA_model_summary.txt):  ASCA analysis summary
	10. [ASCA\_TempAffectedProteins\_loadings.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/ASCA_TempAffectedProteins_loadings.csv):  ASCA proteins and their PC1 loadings values that were selected by PC1 loadings cutoff.
	11. [ASCA_TempAffectedProteins.csv](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/ASCA_TempAffectedProteins.csv):  ASCA selected proteins proteins (listed as columns) and their average NSAF values.
	12. [sumNUMSPECSTOT\_plus1\_ratioFC\_logFC\_pval\_DAYSCOMPARED](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/sumNUMSPECSTOT_plus1_ratioFC_logFC_pval_DAYSCOMPARED.csv):  Log base 2 fold change and Chi squared raw and corrected p-values of all proteins detected across silos
	13. [all_sig0.1ASCAclust\_pro\_logFC\_pval\_abbrv\_Evalcutoff\_NodeAttb.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/all_sig0.1ASCAclust_pro_logFC_pval_abbrv_Evalcutoff_NodeAttb.txt):  Node attribute file for proteins containing unique names for proteins (Uniprot entry gene name with protein name suffix appended to distinguish proteoforms), log fold change values, and FDR-corrected Chi squared proportions test p-values
	14. [BPGOnode\_attb\_semsim0.5\_sig0.1ASCAClust\_EValcutoff.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/BPGOnode_attb_semsim0.5_sig0.1ASCAClust_EValcutoff.txt): Node attribute file for GO slim terms containing total number spectra magnitude foldchange ('TNSmagFC') for each timepoint, total NSAF magnitude foldchange ('NSAFmagFC') for each timepoint, and number of proteins per term.
	15. [edge\_attb\_GOBPsemsim0.5\_sig0.1ASCAClust\_EValcutoff.txt](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Data/edge_attb_GOBPsemsim0.5_sig0.1ASCAClust_EValcutoff.txt): Edge attribute file for GO terms containing term-term relationships and their semantic similarity scores.
- **[Scripts](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Scripts)**: 
	[Cg\_Giga\_cont\_AA.fa\_BLASTP\_uniprot\_swprot2019.ipynb](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Scripts/Cg_Giga_cont_AA.fa_BLASTP_uniprot_swprot2019.ipynb):  Jupyter notebook used to align protein sequences to the Uniprot-KB Swiss-Prot database using BLASTp 
- **[Figures](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Figures)**:
	1. [MainText](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Figures/MainText): 
		1. [Figure 1](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/MainText/Figure1.png): Oyster seed performance at the time of settlement
		2. [Figure 2](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/MainText/Figure2.png): Proteins detected across treatment groups and time points
		3. [Figure 3](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/MainText/Figure3.png): Temperature-affected proteins identified across three statistical methods
		4. [Figure 4](): Biological processes associated with temperature affected proteins
	2. [Supplementary](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/tree/master/Figures/Supplementary):
		1. [Supplementary Figure 1](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure1.jpg): PCA of all technical replicate samples
		2. [Supplementary Figure 2](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure2.jpg): Pie charts showing proportions of GO slim biological processes represented by proteins detected among all a) 23°C proteomes and b) 29°C proteomes
		3.[Supplementary Figure 3](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure2.png): Hierarchical clustering analysis of protein abundances over time for both temperature treatments
		4. [Supplementary Figure 4](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure3.png): ANOVA-simultaneous component analysis of protein abundances over time
		5. [Supplementary Figure 5](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure4.jpg): Proteins selected by cluster analysis
		6. [Supplementary Figure 6](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure5.jpg): Abundance of proteins selected by ANOVA-simultaneous component analysis
		7. [Supplementary Figure 7](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure6.jpg): Abundances of proteins selected by Chi-squared proportions test
		8. [Supplementary Figure 8](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/blob/master/Figures/Supplementary/SupplementaryFigure8.jpg): Biological processes likely affected by temperature-altered protein abundances at each ay post-fertilization (dpf) timepoint 



# Citation
Trigg, S. A., Mitchell, K. M., Elliot, R., Euladiene, B., Vadopalas, B., Timmins-Schiffman, E. B., and Roberts, S. B. (2019) An integrated analysis reveals temperature-influenced proteomic variation throughout early development in the Pacific oyster. 
