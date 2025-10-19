# DiffGeneExpression-DESeq2
# Differential Gene Expression Analysis with DESeq2

A comprehensive R pipeline for analyzing RNA-seq count data using DESeq2 to identify differentially expressed genes between control and experimental conditions. This analysis includes quality control, statistical testing, visualization, and results export.

## Overview

This project performs differential expression analysis on RNA-seq data, comparing control samples against experimental samples. The pipeline includes:
- Automated data filtering
- Quality control visualizations (PCA, dispersion plots)
- Statistical testing for differential expression
- Multiple visualization options (MA plots, volcano plots)
- Flexible filtering criteria
- Automated results export

## Requirements

### R Version
- R version 4.0 or higher recommended
- RStudio (optional but recommended)

### R Packages

```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install DESeq2
BiocManager::install("DESeq2")

# Required packages that come with R
# - stats
# - graphics
# - grDevices
```

### System Requirements
- Minimum 4GB RAM (8GB+ recommended for large datasets)
- ~100MB disk space for results and plots

## Input Data

### File Format: `counts.csv`

The input file should be a CSV with:
- **First column**: Gene identifiers (will be used as row names)
- **Subsequent columns**: Read counts for each sample
- **Header row**: Sample names (required)

#### Example Format:
```
Gene_ID,Control_1,Control_2,Control_3,Sample_1,Sample_2,Sample_3
GENE001,150,200,175,300,400,350
GENE002,50,60,55,45,40,48
GENE003,1200,1150,1180,2300,2400,2350
GENE004,5,8,6,10,12,9
```

### Data Requirements:
- **Count data**: Raw, unnormalized integer counts (NOT TPM, FPKM, or normalized values)
- **No missing values**: All cells should contain numeric values
- **Consistent naming**: Sample names in header must be unique

### Experimental Design:
The code assumes:
- First 3 columns are **control samples**
- Remaining columns are **experimental samples**

**To modify**: Edit the `condition` factor in the setup section.

## Installation

1. **Clone or download** this repository
2. **Install required R packages** (see Requirements section)
3. **Place your data**: Ensure `counts.csv` is in the same directory as the R Markdown file
4. **Update file path**: Modify the path in `read.delim()` if necessary

## Usage

### Option 1: R Markdown (Recommended)

Open `Pset4.Rmd` in RStudio and click "Knit" or run:

```r
rmarkdown::render("Pset4.Rmd")
```

This generates a PDF report with all analyses, plots, and results.

### Option 2: Interactive R Session

Run chunks interactively in RStudio:
1. Open the .Rmd file
2. Click "Run" on each code chunk sequentially
3. View results in the console and plots pane

### Option 3: R Script

Extract code chunks to a .R script and run:
```r
source("deseq2_analysis.R")
```

## Analysis Workflow

### 1. Data Loading and Filtering

```r
Counts <- read.delim("counts.csv", header = TRUE, row.names = 1, sep = ",")
Counts <- Counts[which(rowSums(Counts) > 50), ]
```

**Purpose**: 
- Loads count data
- Filters out lowly expressed genes (total counts ≤ 50)
- Removes noise and improves statistical power

### 2. Experimental Design Setup

```r
condition <- factor(c(rep("C", 3), rep("S", ncol(Counts) - 3)))
coldata <- data.frame(row.names = colnames(Counts), condition)
```

**Customization**:
For specific timepoints or conditions:
```r
condition <- factor(c("Control", "Control", "Control",
                      "Treatment_6h", "Treatment_6h",
                      "Treatment_24h", "Treatment_24h"))
```

### 3. DESeq2 Analysis

```r
dds <- DESeqDataSetFromMatrix(countData = Counts, 
                               colData = coldata, 
                               design = ~condition)
dds <- DESeq(dds)
```

**What happens**:
- Normalization for library size
- Dispersion estimation
- Negative binomial statistical testing

### 4. Quality Control

#### PCA Plot
```r
plotPCA(vsdata, intgroup = "condition")
```
**Interpretation**:
- Controls should cluster together
- Samples should cluster together
- Clear separation = good experimental design

#### Dispersion Plot
```r
plotDispEsts(dds)
```
**Interpretation**:
- Points should follow the fitted curve
- Validates statistical model assumptions

### 5. Results Extraction

```r
res <- results(dds, contrast = c("condition", "S", "C"))
```

**Output Columns**:
- `baseMean`: Average normalized expression
- `log2FoldChange`: Log2(Sample/Control)
- `lfcSE`: Standard error
- `stat`: Wald test statistic
- `pvalue`: Raw p-value
- `padj`: Adjusted p-value (Benjamini-Hochberg)

### 6. Filtering Significant Genes

#### Basic Filter (padj < 0.05)
```r
sigs <- na.omit(res)
sigs <- sigs[sigs$padj < 0.05, ]
```

#### Stringent Filter (padj < 0.05 AND |log2FC| > 1)
```r
sigs_stringent <- sigs[abs(sigs$log2FoldChange) > 1, ]
```

**Thresholds**:
- `padj < 0.05`: 5% false discovery rate
- `|log2FC| > 1`: 2-fold change minimum

## Output Files

### Automatically Generated Files:

1. **Pset4.pdf** (or Pset4.html)
   - Complete analysis report
   - All plots and tables
   - Generated when knitting the R Markdown

2. **deseq2_all_results.csv**
   - All genes tested
   - All statistics included
   - For comprehensive downstream analysis

3. **deseq2_significant_genes.csv**
   - Only genes with padj < 0.05
   - Quick reference for significant hits

4. **deseq2_significant_genes_stringent.csv**
   - Genes with padj < 0.05 AND |log2FC| > 1
   - High-confidence differentially expressed genes

### File Locations:
All output files are saved in the same directory as the R Markdown file.

## Interpreting Results

### Log2 Fold Change
- **Positive values**: Gene is upregulated in samples vs controls
- **Negative values**: Gene is downregulated in samples vs controls
- **log2FC = 1**: 2-fold increase
- **log2FC = -1**: 2-fold decrease
- **log2FC = 2**: 4-fold increase

### Statistical Significance
- **padj < 0.05**: Statistically significant (5% FDR)
- **padj < 0.01**: Highly significant (1% FDR)
- **padj < 0.001**: Very highly significant (0.1% FDR)

### Biological Significance
Combine statistical and biological criteria:
- **padj < 0.05** AND **|log2FC| > 1**: Moderate effect
- **padj < 0.01** AND **|log2FC| > 2**: Strong effect

### Volcano Plot
- **X-axis**: log2 Fold Change (effect size)
- **Y-axis**: -log10(p-value) (significance)
- **Red points**: Significantly upregulated
- **Blue points**: Significantly downregulated
- **Gray points**: Not significant

### MA Plot
- **X-axis**: Average expression level
- **Y-axis**: log2 Fold Change
- Shows relationship between expression and fold change
- Blue points are significant (padj < 0.1)

## Customization

### Adjust P-value Threshold

```r
# More stringent
sigs <- sigs[sigs$padj < 0.01, ]

# Less stringent
sigs <- sigs[sigs$padj < 0.10, ]
```

### Adjust Fold Change Threshold

```r
# Require 4-fold change
sigs_stringent <- sigs[abs(sigs$log2FoldChange) > 2, ]

# Require 1.5-fold change
sigs_stringent <- sigs[abs(sigs$log2FoldChange) > 0.585, ]
```

### Multiple Comparisons

For multiple experimental groups:

```r
# Define more complex conditions
condition <- factor(c("Control", "Control", "Control",
                      "Treatment_A", "Treatment_A",
                      "Treatment_B", "Treatment_B"))

# Compare specific contrasts
res_A_vs_C <- results(dds, contrast = c("condition", "Treatment_A", "Control"))
res_B_vs_C <- results(dds, contrast = c("condition", "Treatment_B", "Control"))
res_A_vs_B <- results(dds, contrast = c("condition", "Treatment_A", "Treatment_B"))
```

### Change Plot Appearance

```r
# Adjust volcano plot colors
with(subset(res, padj < 0.05 & log2FoldChange > 1), 
     points(log2FoldChange, -log10(pvalue), pch = 20, col = "darkred"))

# Change MA plot limits
plotMA(res, ylim = c(-10, 10))
```

## Troubleshooting

### Issue: "Error in read.delim: cannot open file"
**Solution**: 
- Check file path is correct
- Ensure `counts.csv` is in the working directory
- Use `getwd()` to check current directory
- Use absolute path: `/full/path/to/counts.csv`

### Issue: "All gene-wise dispersion estimates are within 2 orders of magnitude"
**Solution**: This is usually fine, just a warning. Proceed normally.

### Issue: Too few significant genes (< 10)
**Possible causes**:
- Weak treatment effect
- High biological variability
- Small sample size
- Incorrect experimental design specification

**Solutions**:
- Verify sample labels are correct
- Check PCA plot for separation
- Consider less stringent thresholds (padj < 0.1)
- Increase sample size if possible

### Issue: Too many significant genes (> 5000)
**Possible causes**:
- Very strong treatment effect
- Technical artifacts or batch effects
- Incorrect normalization

**Solutions**:
- Check PCA for outliers
- Apply more stringent filtering (add fold change threshold)
- Verify samples are labeled correctly

### Issue: "NA values in results"
**Cause**: 
- Genes with all zero counts
- Extreme outliers
- Cook's distance filtering

**Solution**: Use `na.omit()` as shown in code (already implemented)

### Issue: Poor PCA separation
**Cause**:
- Weak treatment effect
- High variability within groups
- Batch effects
- Mislabeled samples

**Solutions**:
- Verify sample labels
- Check for batch effects (add batch to design)
- Consider if samples are truly different

## Best Practices

1. **Always check quality control plots** before interpreting results
2. **Use adjusted p-values** (padj), never raw p-values
3. **Combine statistical and biological significance** (padj + fold change)
4. **Document your filtering criteria** clearly
5. **Keep raw count data** - never use normalized counts as input
6. **Report number of replicates** for reproducibility
7. **Save your R session**: `save.image("analysis.RData")`

## Common Analysis Extensions

### Gene Set Enrichment Analysis
After identifying significant genes, perform pathway analysis:
```r
# Export gene list
write.table(rownames(sigs), file = "significant_genes.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
```
Use with tools like:
- DAVID (https://david.ncifcrf.gov/)
- Enrichr (https://maayanlab.cloud/Enrichr/)
- Gene Ontology (http://geneontology.org/)

### Heatmap of Top Genes
```r
# Get top 50 most significant genes
top_genes <- head(order(res$padj), 50)

# Extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)[top_genes, ]

# Create heatmap
library(pheatmap)
pheatmap(log2(norm_counts + 1), 
         scale = "row",
         show_rownames = TRUE)
```

### Time Course Analysis
For multiple timepoints, use likelihood ratio test:
```r
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res_lrt <- results(dds_lrt)
```

## References

### Primary Citation
Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15:550.
https://doi.org/10.1186/s13059-014-0550-8

### Documentation
- DESeq2 Vignette: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- Bioconductor Page: https://bioconductor.org/packages/DESeq2/
- User Support: https://support.bioconductor.org/

### Additional Resources
- RNA-seq analysis guide: https://www.bioconductor.org/help/course-materials/
- DESeq2 workflow: https://www.bioconductor.org/packages/release/workflows/html/rnaseqGene.html

## License

This analysis pipeline is provided for educational and research purposes.

## Citation

If you use this pipeline in your research, please cite:
- The DESeq2 paper (see References)
- This repository (if publicly shared)

## Support

For issues with:
- **DESeq2 package**: Post on Bioconductor support forum
- **This pipeline**: Check Troubleshooting section above
- **Statistical interpretation**: Consult with a biostatistician

## Version History

- **v1.0** (2024-04-13): Initial version with basic differential expression analysis
- **v1.1** (Current): Added volcano plots, MA plots, enhanced filtering, and automated export

## Author Notes

This pipeline was developed for Problem Set 4 (Pset4) and represents a standard differential expression workflow. Modify the experimental design section to match your specific study design.
