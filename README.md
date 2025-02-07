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
## Save the cleaned DESeq2 results to a CSV file
```
write.csv(as.data.frame(res_ordered), file = "DESeq2_results_cleaned.csv")
```
## Verify the file was created
```
file.exists("DESeq2_results_cleaned.csv")

```
## Save the cleaned DESeq2 results to a CSV file
```
write.csv(as.data.frame(res_ordered), file = "DESeq2_results_cleaned.csv")
```
## Verify the file was created
```
file.exists("DESeq2_results_cleaned.csv")
```
## Load control data from the Excel files
```
library(readxl)
df1 <- read_excel("1_S27_dataframe.xlsx")
df2 <- read_excel("2_S28_dataframe.xlsx")
```
## Clean control data
```
df1_clean <- df1[, c("V1", "V3")]
df2_clean <- df2[, c("V1", "V3")]
colnames(df1_clean) <- c("Gene", "Control1")
colnames(df2_clean) <- c("Gene", "Control2")
```
## Load treatment data
```
treatment_data <- read.table("3_S29counts.txt", header = FALSE, sep = "\	", 
                           col.names = c("Gene", "Symbol", "Treatment2"))
treatment_data_clean <- treatment_data[, c("Gene", "Treatment2")]
```
## Merge all data
```
merged_df <- merge(df1_clean, df2_clean, by = "Gene", all = TRUE)
merged_df <- merge(merged_df, treatment_data_clean, by = "Gene", all = TRUE)
```
## Remove rows with missing values
```
merged_df <- na.omit(merged_df)
```
## Print dimensions of merged data
```
print("Dimensions of merged dataset:")
dim(merged_df)
```
## Show first few rows
```
print("\
First few rows of merged dataset:")
head(merged_df)

```

## Install BiocManager if not already installed
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
```
## Install DESeq2
```
BiocManager::install("DESeq2")
```
##Load library
```
library(DESeq2)
```
## Run DESeq2 analysis
```
dds <- DESeqDataSetFromMatrix(
    countData = round(counts_matrix),
    colData = colData,
    design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds)
res_ordered <- res[order(res$padj), ]
```
## Save results
```
write.csv(as.data.frame(res_ordered), "DESeq2_results_with_new_treatment.csv")
```
## Print summary
```
print("Summary of DESeq2 results:")
summary(res)

print("\
Top 10 differentially expressed genes:")
head(res_ordered, 10)
```
