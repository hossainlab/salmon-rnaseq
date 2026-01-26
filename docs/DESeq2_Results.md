# Interpretation of Differential Expression Result Columns (DESeq2)

The following columns are standard outputs from **DESeq2** differential expression analysis. They quantify gene expression abundance, effect size, statistical uncertainty, and significance.

## 1. `baseMean`
The average of **normalized read counts** for a gene across all samples.


**Interpretation**

- Reflects overall gene expression level
- Independent of experimental condition
- Low values indicate lowly expressed genes

**Practical use**

- Used for filtering unreliable low-count genes  
- Common cutoff: `baseMean > 10`

## 2. `log2FoldChange`

The log2-transformed ratio of expression between two conditions.

**Interpretation**

- `+1` → 2-fold upregulation
- `−1` → 2-fold downregulation
- `0` → no change

**Notes**
- In DESeq2, fold changes are **shrunken** using empirical Bayes methods
- Shrinkage stabilizes estimates for low-count genes

## 3. `lfcSE` (Log Fold Change Standard Error)
The standard error associated with the estimated `log2FoldChange`.

**Interpretation**

- Smaller values indicate more precise estimates
- Larger values indicate higher uncertainty

## 4. `stat`
The **Wald test statistic** used to test the null hypothesis:

**Interpretation**
- Large absolute values indicate strong evidence for differential expression
- Sign reflects direction of regulation

## 5. `pvalue`
The unadjusted p-value derived from the Wald test.


**Interpretation**

- Probability of observing the test statistic under the null hypothesis
- Not corrected for multiple testing

**Limitation**

- Inflated false positives when testing thousands of genes

## 6. `padj` (Adjusted p-value)

The p-value adjusted for multiple testing using the  
**Benjamini–Hochberg False Discovery Rate (FDR)** method.

**Interpretation**

- Controls expected proportion of false discoveries
- Common significance threshold: `padj < 0.05`
