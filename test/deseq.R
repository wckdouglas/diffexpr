#!/usr/bin/env Rscript

library(DESeq2)

df <- read.table('data/ercc.tsv', header=T) 

sample_df <- data.frame(samplename = names(df)[2:ncol(df)])
sample_df$sample <- substr(sample_df$samplename, 1,1)
sample_df$batch <- substr(sample_df$samplename, 3,3)
row.names(sample_df) <- sample_df$samplename
sample_df$sample <- factor(sample_df$sample)


count_matrix <- subset(df, , -id)
row.names(count_matrix) <- df$id
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                       colData = sample_df, 
                       design = ~ batch + sample)
dds <- DESeq(dds)
res <- results(dds)
res$id = row.names(res)
write.table(data.frame(res), 
            file = 'data/R_deseq.tsv', 
            sep='\t', 
            row.names=F,
            quote=F)

dds <- DESeq(dds, test="LRT", reduced=~ batch)
res <- results(dds)
res$id = row.names(res)
write.table(data.frame(res), 
            file = 'data/R_deseq_reduced.tsv', 
            sep='\t', 
            row.names=F,
            quote=F)