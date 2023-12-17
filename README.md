# GEA
#Gene Expression Analysis
# Download Data from cbioportal
# Change this to your own directory.
path = "/Users/sallyzamuhammad/Downloads"

# untar folder.
folder_name = "brca_tcga_pan_can_atlas_2018.tar.gz"
folder = paste(path, folder_name, sep = "/")
untar(folder)

# Go to new path
new_dir = paste(getwd(),"brca_tcga_pan_can_atlas_2018", sep = "/" )
setwd(new_dir)
new_dir

#(3)Read RNASeq File
rna = read.delim("data_mrna_seq_v2_rsem.txt")

# in this assignment we will delete the genes for which there's more than one Hugo Symbol
# These are typically genes with no Hugo Symbol ("" as an entry) or pseudogenes.
# This is more for simplicity.If you keep your analysis would still be correct so no worries.
keep = !duplicated(rna[,1])
rna = rna[keep,]
View(rna)
# set rownames of rnaseq to hugo symbols
rownames(rna)  = rna[,1]

#Read Copy Number Aberrations Data 
cnadata = read.delim("data_cna.txt")
# find ERBB2 in cna
erbb2_indx = which(cnadata[,1] == 'ERBB2')
# Plot histogram to visualize explore the data.
hist(as.numeric(cnadata[erbb2_indx,-c(1,2)]))
View(cnadata)

##Match cna data and rna data 
rna_cna_id = which(is.element(colnames(rna[,-c(1,2)]), colnames(cnadata[,-c(1,2)])))

# select only the rna cases which have cna data.
rna_cna_sub = rna[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna
no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(rna[,2+rna_cna_id]), colnames(cnadata[,-c(1,2)]))) 

# sanity check.This will print an error if the result is not the same.
sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]

# Pre-allocate memory for ERBB2
meta_erbb2 = matrix(0,length(rna_cna_id),1)

for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cnadata)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cnadata[erbb2_indx,col_cna]>0)
  
}

# This are some checks you can do to make sure your code worked.
# There's some more systematic checks you can do. See unit testing.
# simple checks to make sure. 

col_i = colnames(rna_cna_sub)[1]
col_cna = which(colnames(cnadata)==col_i)

# sanity check
(cnadata[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.
pos_example = which(meta_erbb2==1)[1]
col_i = colnames(rna_cna_sub)[pos_example]
col_cna = which(colnames(cnadata)==col_i)


# sanity check
(cnadata[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]


# botch checks should print true.
# We will add a title to the metadata.
colnames(meta_erbb2) = "ERBB2Amp"
# transform into intege


BiocManager::install("ReactomePA")
BiocManager::install("pathview")
BiocManager::install("msigdbr")
BiocManager::install("ggnewscale")

library(ggnewscale)
library(pathview)
library(dplyr)
library(tidyverse)
library (Geoquery)


# BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Install DeSeq2
BiocManager::install("DESeq2")
library(DESeq2)

# Build DESeq Object

dds <- DESeqDataSetFromMatrix(countData = round(rna_cna_sub),
                              colData = meta_erbb2,
                              design = ~ ERBB2Amp)

# Filter
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]


# Normalize
dds <- DESeq(dds)

# Transform the data to visualize
rld <- vst(dds, blind=FALSE)

# Do Principal Components Analysis
pc = prcomp(assay(rld))
head(pc)

# Plot 
plot(pc$rotation[,1], pc$rotation[,2], col = 1+(meta_erbb2), pch = 19)

# Get Results - perform statistical analyst to identify differentially expressed genes.
res <- results(dds)
res

# Summary

summary(res)
rownames(res) = rna[keep,1]
rownames(res)

# Significantly Differentially Expressed
signif = which(res$padj<0.05)
deg = res[signif,]
deg


# Separate them 
dup = deg[deg[,2]>0.,]
head(dup)
ddown = deg[deg[,2]<0.,]
head(ddown)

# For Pathway Enrichment we need Entrez IDs

entrez_ids = rna[keep,2]
entrez_all = entrez_ids[signif]
entrez_up = entrez_all[signif[deg[,2]>0.]]
entrez_down = entrez_all[signif[deg[,2]<0.]]

# Pathway Enrichment

BiocManager::install("clusterProfiler")

library(clusterProfiler)

# Do a KEGG pathway over-representation analysis
all_paths =   enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.05)
head(all_paths)

# Optionally you can divide between up and down.
# Both options are Ok for the assignment.

up_paths = enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.05)
head(up_paths)

down_paths = enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.05)
head(down_paths)

###P Value
resOrdered <- res[order(res$pvalue),]
resOrdered
summary(res)

