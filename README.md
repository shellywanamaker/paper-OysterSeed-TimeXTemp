# Integrated analysis reveals temperature-influenced proteomic variation throughout early development in the Pacific oyster



## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [Issues](https://github.com/shellytrigg/paper-OysterSeed-TimeXTemp/issues)
- [Citation](#citation)

# Overview
The Pacific oyster has major ecological and economic importance serving as a biofilter and habitat in coastal ecosystems, and contributing over $190M to annual marine aquaculture revenue. However, little is known about the landscape of protein expression during early development, a time when mass mortality is common which can negatively impact industry and ecosystems. To better characterize physiological pathways and associated networks active during oyster development we performed a developmental time series proteomics analysis of larval cultures reared at 23°C and 29°C. These temperatures were selected based on the reports from the aquaculture industry that differential performance is observed in oysters at these temperatures. In addition to observing larger sized animals at 29°C, we found differentially abundant proteins that point to upregulated systemic structural remodeling during early development at 23°C and upregulated growth during late development at 29°C. The differential proteomes suggest that at 23°C more cellular energy is being diverted to maintenances processes where at 29°C more cellular energy is being used for growth. This proteomics analysis combined with development observations offers greater clarity on environmental conditions that can improve aquaculture production.

# Repo Contents
- **[R]()**: R code
	1. 	[ClusteringTechnicalReplicates\_and\_CalculateAvgNSAFs.R]():   R code used to cluster technical replicates and calulate average NSAF values

	2. [General\_proteome\_characterization.Rmd]():  R code used to  convert the GO identifiers associated with each protein to GO slim terms using the MGI GO Slim database (http://www.informatics.jax.org/gotools/MGI_GO_Slim.html), and plot pie chart of biological processes.  
	3. [clustering-dendrogram\_s3\_9.R]():  R code used for hierarchical clustering do identify differentially abundant proteins
	4. [ASCA.R]():  R code used to run ASCA and identify differentially abundant proteins
	5. [LogFC\_ChiSquaredTest\_SameDaySamples.Rmd](): R code used for calculating log foldchange and Chi squared proportions test on all protein 
	6. [VerifyStatsProteinSelection.Rmd](): R code that selects proteins with < 0.1 adjusted Chi Squared p-value, makes venn diagram (Fig. 2) plot to show the overlap between all statistical methods, plots protein abundance for each statistical method as facetted line plots and as heatmaps (Supplementary Figs. 2-4)
	7. [CreateNodeAndGoSemSimEdgeAttr.R](): R code that retrieves GO slim identifiers for temperature-affected proteins using GSEA, calculates and thresholds GO semantic similarity scores for GO slims, calculates the total magnitude foldchange NSAF and total number of spectra for each term based on all proteins associated with each term, generates node and edge attribute files for Cytoscape.


- **[Data]()**: 
	1. [BiovolumeAfterSettlement.xlsx]():  Biovolume of oyster seed (mL) settled after 6 days of exposure to high temperatures.
	2. [Cg\-Giga\_cont\_AA.fa]():  Protein sequence database file ('Cg-Giga\_cont\_AA.fa') containing FASTA sequence file of the C. gigas proteome (downloaded from http://gigaton.sigenae.org as 'contigs.fasta.transdecoder.pep.gz') and common contaminants downloaded from the crapOME (Mellacheruvu et al. Nat Methods 2013). This file is used by proteomics spectral analysis and by UniProt mapping.
	2. [ABACUS\_output.tsv.gz]():  Output file from ABACUS analysis  
	3. [Sample_metadata.csv]():  
	4. [silo3and9\_NSAF\_AVGs.csv]():  Average NSAF values for all proteins (listed as rows) detected in all samples (listed as columns)
	5. [silo3and9\_nozerovals\_NSAF\_AVGs.csv]():  Average NSAF values for all proteins (listed as columns) detected in all samples (listed as rows) with zero value converted to 0.1. 
	6. [PreliminaryProteomeCharacterization.xlsx](): Data and plot showing overlap between different temperature proteomes at each timepoint underlying Figure 2a.  
	6. [HC_UniquelyClusteredProteins.csv]():  Proteins identified by hierarchical clustering as showing differential abundance between temperature treatments. 
	7. [ASCA\_model\_summary.txt]():  ASCA analysis summary
	8. [ASCA\_TempAffectedProteins\_loadings.csv]():  ASCA proteins and their PC1 loadings values that were selected by PC1 loadings cutoff.
	9. [ASCA_TempAffectedProteins.csv]():  ASCA selected proteins proteins (listed as columns) and their average NSAF values.
	10. [sumNUMSPECSTOT\_plus1\_ratioFC\_logFC\_pval\_DAYSCOMPARED]():  Log base 2 fold change and Chi squared raw and corrected p-values of all proteins detected across silos
	11. [all_sig0.1ASCAclust\_pro\_logFC\_pval\_abbrv\_Evalcutoff\_NodeAttb.txt]():  Node attribute file for proteins containing unique names for proteins (Uniprot entry gene name with protein name suffix appended to distinguish proteoforms), log fold change values, and FDR-corrected Chi squared proportions test p-values
	12. [BPGOnode\_attb\_semsim0.5\_sig0.1ASCAClust\_EValcutoff.txt](): Node attribute file for GO slim terms containing total number spectra magnitude foldchange ('TNSmagFC') for each timepoint, total NSAF magnitude foldchange ('NSAFmagFC') for each timepoint, and number of proteins per term.
	13. [edge\_attb\_GOBPsemsim0.5\_sig0.1ASCAClust\_EValcutoff.txt](): Edge attribute file for GO terms containing term-term relationships and their semantic similarity scores.
- **[Scripts]()**: 
	1. [DDA-data-Analyses.md]():  DDA proteomics analysis workflow that converts files from .raw to .mzXML, searches .mzXML files against C. gigas database, calculates statistics associated with peptide and protein IDs using the Trans Proteomic Pipeline, and correlates protein inferences using ABACUS.  
	2. [Pcomet.params.high-low.txt]():  COMET 2016.01 'comet.params.high-low' search parameters used in COMET analysis of mass spectrometry data download from http://comet-ms.sourceforge.net/parameters/parameters_201601/comet.params.high-low
	3. [Abacus_parameters.txt]():  ABACUS parameter file used in ABACUS analysis of mass spectrometry data
	4. [Cg\_Giga\_cont\_AA.fa\_BLASTP\_uniprot\_swprot2019.ipynb]():  Jupyter notebook used to align protein sequences to the Uniprot-KB Swiss-Prot database using BLASTp 
- **[Figures]()**:
	- MainText 
		1. [Figure 1](): Oyster seed performance at the time of settlement
		2. [Figure 2](): Proteins detected across treatment groups and time points
		3. [Figure 3](): Temperature-affected proteins identified across three statistical methods
		4. [Figure 4](): Biological processes associated with temperature affected proteins
	- Supplementary_Figures
		1. [Supplementary Figure 1](): PCA of all technical replicate samples
		2. [Supplementary Figure 2](): Hierarchical clustering analysis of protein abundances over time for both temperature treatments
		3. [Supplementary Figure 3](): ANOVA-simultaneous component analysis of protein abundances over time
		4. [Supplementary Figure 4](): Proteins selected by cluster analysis
		5. [Supplementary Figure 5](): Abundance of proteins selected by ANOVA-simultaneous component analysis
		6. [Supplementary Figure 6](): Abundances of proteins selected by Chi-squared proportions test



# Citation
Trigg, S. A., Mitchell, K. M., Elliot, R., Euladiene, B., Vadopalas, B., Timmins-Schiffman, E. B., and Roberts, S. B. (2019) Integrated analysis reveals temperature-influenced proteomic variation throughout early development in the Pacific oyster. 
