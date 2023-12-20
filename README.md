# Breast-Cancer-Analysis
### Deseq2: Gene Expression Data Preprocessing and Differential Gene Analysis
- **Purpose**: To preprocess breast cancer gene expression data and perform differential gene analysis.
- **Key Steps**:
  - Loading relevant R packages, such as `DESeq2`, `ggplot2`, etc.
  - Reading and preprocessing gene expression data.
  - Standardizing and analyzing differential genes using DESeq2.
  - Selecting genes with significant differential expression.

### kegg_go: Enrichment Analysis (Upregulated Genes)
- **Purpose**: To perform enrichment analysis on upregulated differential genes in breast cancer.
- **Key Steps**:
  - Reading differential gene data and categorizing it.
  - Using the `org.Hs.eg.db` package for gene ID conversion.
  - Performing KEGG and GO enrichment analysis on upregulated genes.
  - Visualizing the results of the enrichment analysis using `ggplot2` and other tools.

### cox: Enrichment Analysis (Downregulated Genes)
- **Purpose**: To perform enrichment analysis on downregulated differential genes in breast cancer.
- **Key Steps**:
  - Repeating the enrichment analysis steps for upregulated genes, but focusing on downregulated genes.
  - Conducting KEGG and GO enrichment analysis and visualizing the results.

### Volcano Plot_Heatmap: Volcano Plot and Heatmap Creation
- **Purpose**: To visualize the results of the differential gene analysis.
- **Key Steps**:
  - Using a volcano plot to intuitively display the significance and expression changes of differential genes.
  - Using a heatmap to show the expression patterns of specific differential genes across different samples.
  - Saving the visualization results in PDF format.

### Overall Workflow:
1. **Data Preprocessing and Differential Gene Analysis**: Starting from raw data, the gene expression data is cleaned, standardized, and differential expression genes are identified.
2. **Enrichment Analysis**: Conducting KEGG and GO enrichment analyses separately for upregulated and downregulated genes to explore the potential functions and pathways associated with these genes in biological processes.
3. **Data Visualization**: Presenting the results of differential gene analysis in a visual format through volcano plots and heatmaps for a better understanding and interpretation of the data.

Together, these codes provide researchers with a complete toolkit for analyzing and understanding the differences in gene expression in breast cancer, which is crucial for discovering potential biomarkers and understanding disease mechanisms.

The dataset can be downloaded from this link: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018