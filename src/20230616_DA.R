
### to work on a cluster, need to load packages within script, installed using R interactive mode
# The code combines data loading, manipulation, and file exports.
### Run on the cluster using split_contrasts.sh from the code folder: 

## 1) connect to graham on your terminal: ssh -vvv -Y username@graham.computecanada.ca
## 2) Navigate to the folder: cd projects/def-profname-ab/yourusername/
## 3) Create a folder "test": mkdir test

#### In a new terminal tab: 
## 2) transfer the 20230616_DA.R and the 20230616_DA.sh scripts as well as the data 
## to the test folder using the scp or sftp protocols as described in notes.txt

## 3) connect to graham on your terminal: ssh -vvv -Y username@graham.computecanada.ca
## 4) Navigate to the folder: cd projects/def-profname-ab/youruser/test
## 5) Submit the job: sbatch 20230616_DA.sh




library(here)
library(tidyverse)
library(DESeq2)

# ------------------ Load data ------------------
com_path <- here('data', 'com.csv')
meta_path <- here('data', 'metadata.csv')

meta <- read.csv(meta_path, sep = ",",comment.char = "", header = T, row.names = 1)
com <- read.csv(com_path, sep = ",",comment.char = "", header = T, row.names = 1)

# ------------------ Clean data ------------------
tax <- as.data.frame(com$taxonomy)
colnames(tax) <- "taxonomy"
rownames(tax) <- rownames(com)

com <- com[, 1:ncol(com)-1]
# Get the matching column indices based on row names
matching_columns <- match(rownames(meta), colnames(com))

# Reorder the columns in 'com' based on the matching column indices
com <- com[, matching_columns]

rm(com_path)
rm(matching_columns)
rm(meta_path)

print('concluded data cleaning... Starting DA')
# ------------------ Differential analysis ------------------

dds <- DESeqDataSetFromMatrix(countData = com,
                              colData = meta,
                              design = ~ contamination)
# Perform normalization
dds <- DESeq(dds)

# Conduct differential analysis
res <- results(dds)
res_df <- as.data.frame(res)
res_df$asv_id <- rownames(res_df)

# Volcano plot
volcano <- ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(size = 2, color = ifelse(res_df$padj < 0.05, "red", "black"), alpha = 0.7) +
    xlim(c(-5, 5)) +
    ylim(c(0, 10)) +
    labs(x = "Log2 Fold Change",
         y = "-Log10 p-value",
         title = "Volcano Plot",
         subtitle = "Differential Abundance Analysis") +
    theme_light() +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 12),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none")
#volcano
print('concluded DA ... Exporting results')
# ------------------ Export results ------------------

# Save the plot as a high-resolution image (e.g., PDF)
ggsave(here("out","volcano_plot.pdf"), plot = volcano, width = 6, height = 4, units = "in", dpi = 300)
write.csv(res_df, file = here('out', 'DA_results.csv'))


print('concluded exporting results ... analysis complete')

