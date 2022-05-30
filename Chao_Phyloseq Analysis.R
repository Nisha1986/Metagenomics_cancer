if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
a
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")

#load libraries
library("phyloseq")
library("ggplot2")
library("dplyr")
getwd()
setwd("C:/Users/91809/Desktop/OMICS World/bioinformatics-master")
otu_mat<- read.table("Chao_fecal_sample_Abundance_table.txt", header = TRUE, sep = "\t")
tax_mat<- read.table("Chao_fecal_sample_Taxonomy_table.txt", header = TRUE, sep = "\t")
samples_df <- read.table("CHOW_HFD_supplementaryTable.txt", header = TRUE, sep = "\t")
row.names(otu_mat) <- otu_mat$OTUs#define row name each data
#Remove the column which is  already used as row from the data
otu_mat <- otu_mat %>% select (-OTUs)
row.names(tax_mat) <- tax_mat$class
tax_mat <- tax_mat %>% select (-class)
row.names(samples_df) <- samples_df$Run
samples_df <- samples_df %>% select (-Run)
sampletype <- unique(row.names(samples_df))
# Transform into matrix
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#Transform matrix data as input for Phyloseq
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
View(samples)
View(OTU)
View(TAX)
chao <- phyloseq(OTU, TAX, samples)
View(chao)
#Visualization: bar chart
plot_bar(chao, fill = "Phylum")
+ theme(legend.key.size = unit(0.3, "cm"), legend.text = element_text( color="Black", size=7.5), axis.text.x = element_text( color="Black", size=7.5))
#NMDS - non-metric multidimensional scaling ordination
chao.ord <- ordinate(chao, "NMDS", "bray")
sample_variables(chao)
#Richness plot
pp=plot_richness(chao,color="LibraryName",measures=c("Chao1", "Shannon"))
pp+theme(legend.key.size = unit(0.3, "cm"), legend.text = element_text( color="Black", size=7.5), axis.text.x = element_text( color="Black", size=6.0))

s1 <- cbind(samples,samples)
# Rename column names
colnames(s1) <- c("LibraryName","copy")
plot_ordination(chao, chao.ord, type="s1", color="LibraryName")
#Split by groups
plot_ordination(chao, chao.ord, type="split", color="LibraryName")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")
a
library(dada2)







