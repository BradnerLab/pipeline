#RNA-SEQ PIPELINE README

‘’’
Description of processing steps and output steps in Bradner Lab RNA-SEQ pipeline
‘’’
I. Summary
II. Input Data
III. Normalization and Gene Expression Quantification
IV. Data Files
V. Quality Control
VI. Pairwise Analysis

I.Summary:

This pipeline takes transcriptome aligned RNA-Seq data and quantifies common gene name level expression values with or without cell count normalized ERCC spike-in addition. Additional replicate quality control metrics as well as pair-wise sample comparisons are performed.

II. Input Data:

Currently pipeline supports any plain text input table with header where the first column is a gene name and subsequent columns are individual samples (e.g. genes.fpkm.table out table from Cufflinks). Units are in coverage and transcript length normalized units (e.g. FPKM or RPKM).

A name for the project is used in the analysis and provides the prefix for all output files (replacing * in the subsequent text)

Sample names must reflect a group structure (e.g., A_1, A_2, B_1, B_2) where all samples in a group share a prefix. These will be treated as replicates (A_1 and A_2) and (B_1 and B_2).

III. Normalization and Gene Expression Quantification:

1. To aid in transformations and normalization, a lower bound of 0.01 is set on all expression levels
2. If ERCC spike-ins are added, a loess normalization is applied using all ERCC probes as a subset. The file ERCC_Controls_Analysis.txt contains information about the ERCC mix including individual spike-in concentrations.
3. Expression values are then filtered based on expression level, keeping only genes with an expression >1 FPKM in at least one sample.
4. Mean expression values are summarized across all genes across group.

IV. Data Files:

1. *_all_fpkm_exprs_raw.txt - This is the raw input
2. *_all_fpkm_exprs_norm.txt - If ERCC normalizations are used, this is cell count normalized data
3. *_all_fpkm_means.txt - Mean summarized expression data across groups (if ERCC spike-ins are used, then derived from *_all_fpkm_exprs_norm.txt)
4. *_exprs_fpkm_means.txt - Mean summarized expression data across groups filtered for expressed genes (if ERCC spike-ins are used, then derived from *_all_fpkm_exprs_norm.txt)

V. Quality Control:

1. *_all_fpkm_exprs_raw_scatter.png - all sample pairwise scatter plots of raw expression data
2. *_all_fpkm_exprs_norm_scatter.png - all sample pairwise scatter plots of ERCC normalized expression data
3. *_spike_raw.pdf - scatter plot showing expression levels of ERCC spike-ins (y-axis) versus their input concentration (x-axis). A loess regression for each sample is plotted as a line. This should illustrate differences in sample level spike-in abundance.
4. *_spike_norm.pdf - scatter plot showing expression levels of ERCC spike-ins (y-axis) versus their input concentration (x-axis) post Loess Normalization. A loess regression for each sample is plotted as a line. If normalization is performed correctly, all lines should overlap.
5. *_exprs_boxplot.pdf - Boxplots showing the distribution of all expression values in each sample for raw data (left) or ERCC cell count normalized data (right). Helpful in quickly determining if spike-ins were added in equivalent amoutns to replicates and in identifying global changes in gene expression.
6. *_replicate correlations.pdf - Replicate scatter plots of normalized expression values filtered by expressed genes. All pairwise comparisons of replicates are plotted by group. Helpful in quickly identifying poor quality samples.

VI. Pairwise Analysis:

For each pairwise comparison of groups (e.g. group A vs. group B), the following outputs are produced.

1. *_A_vs_B_exprs_matrix.txt - Table of all expressed genes with the mean expression in group A and B, the log2 fold change, and the p-Value of their difference by a two-tailed t-test.
2. *_A_vs_B.cls and *_A_vs_B.gct - GSEA formatted input files (.cls and .gct) to perform leading edge analysis on any pairwise comparison
3. *_A_vs_B.pdf - Multipage analysis summary to identify differential gene expression and global changes in gene expression.
	i. Volcano scatter plot of all expressed genes with the log2 fold change in expression (x-axis) plotted vs. the significance of the difference in gene expression (y-axis). Genes that exceed a log2 1 fold change and a p-value of 0.05 are colored blue and red. For each direction, top 10 most differential genes in terms of magnitude of fold change and significance are highlighted.
	ii. Scatter plot of expressed genes in A (y-axis) vs. B (x-axis). Differential genes from same criteria as i are colored
	iii. Waterfall plot of all genes ranked by fold change in A vs. B. No significance criteria is applied here.
	iv. Rank ordered plots of gene expression to detect global change in gene expression or normalization artifacts. Left: Genes are ranked by expression in A and plotted (A: grey, B: red). A loess regression is shown for expression values in B (red line). Right: the opposite pairwise comparison. Genes are ranked by expression in B and plotted (B: grey, A: red). A loess regression is shown for expression values in A (red line).




