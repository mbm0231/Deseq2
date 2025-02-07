# Deseq2
### This is an R script
Differential gene expression
## Create a proper count matrix with unique gene entries. First, let's get unique genes and their mean counts per condition.
```
cleaned_df <- aggregate(cbind(Control1, Control2, Treatment1) ~ Gene, 
                       data = merged_df, 
                       FUN = function(x) mean(x, na.rm = TRUE))
```
## Remove rows where any condition has NA
```
cleaned_df <- na.omit(cleaned_df)
```
## Create the count matrix
counts_matrix <- as.matrix(cleaned_df[, -1])
rownames(counts_matrix) <- cleaned_df$Gene

## Create the column data
```
colData <- data.frame(
    condition = factor(c("Control", "Control", "Treatment"))
)
rownames(colData) <- colnames(counts_matrix)
```
## Run DESeq2 analysis
```
dds <- DESeqDataSetFromMatrix(
    countData = round(counts_matrix),
    colData = colData,
    design = ~ condition
)
```
## Run the analysis
```
dds <- DESeq(dds)
```
## Get results
```
res <- results(dds)
res_ordered <- res[order(res$padj), ]
```
## Save to CSV
```
write.csv(as.data.frame(res_ordered), "DESeq2_results_cleaned.csv")
```
## Show summary
```
print("Summary of DESeq2 results:")
summary(res)
print("\
Top 10 differentially expressed genes:")
head(res_ordered, 10)
```
