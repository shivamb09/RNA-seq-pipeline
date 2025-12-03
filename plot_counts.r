#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
counts_file = args[1]
out_pdf = args[2]

library(ggplot2)
library(reshape2)
library(gridExtra)
library(pheatmap)

fc = read.delim(counts_file, comment.char="#", check.names=FALSE)

std_cols = c("Geneid","Chr","Start","End","Strand","Length")
sample_cols = setdiff(colnames(fc), std_cols)

counts = as.matrix(fc[, sample_cols])
rownames(counts) = fc$Geneid

totals = colSums(counts)
tot_df = data.frame(sample=names(totals), total=as.numeric(totals))

p1 = ggplot(tot_df, aes(x=sample, y=total)) +
  geom_bar(stat="identity") +
  theme_minimal() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ylab("Total counts") + xlab("Sample") + ggtitle("Total mapped counts")

cpm = sweep(counts, 2, colSums(counts), "/") * 1e6
var_genes = head(order(apply(cpm,1,var), decreasing=TRUE), 30)
mat = log2(cpm[var_genes,] + 1)

pdf(out_pdf, width=10, height=8)
grid.arrange(p1)
pheatmap(mat, cluster_rows=TRUE, cluster_cols=TRUE,
         main="Top 30 variable genes (log2 CPM+1)")
dev.off()
