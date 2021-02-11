

# Install Packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximportData")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")

library(tximportData)
library(tximport)
library(DESeq2)

samplesheet <- read.csv("con_ddt_results/pipeline_info/samplesheet.valid.csv")
samples <- strsplit(samplesheet$sample, split = "_T1")

files <- file.path("con_ddt_results/star_rsem", paste0(samples, ".genes.results"))
file.exists(files)

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

txi.rsem$length[txi.rsem$length == 0] <- 1 # https://support.bioconductor.org/p/92763/
head(txi.rsem$counts)

dim(txi.rsem$counts)

colnames(txi.rsem$counts) <- samples
head(txi.rsem$counts)

sampleTable <- data.frame(condition = c(rep("Control", 29), rep("DDT", 34)))

rownames(sampleTable) <- colnames(txi.rsem$counts)

dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
dds <- DESeq(dds)
res <- results(dds)
res

summary(res)

write.csv(as.data.frame(res),file="con_ddt_results/star_rsem/con_ddt_deseq2_diffexp.csv")