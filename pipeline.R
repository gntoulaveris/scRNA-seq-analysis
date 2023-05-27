
# This is an automated scRNA-seq data analysis pipeline.

# Packages ----
#install.packages("Seurat")
#install.packages("tidyverse")
#install.packages("reshape2")
#install.packages("mclust")
#install.packages('installr')
#install.packages("rtools")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("Matrix")
#install.packages("factoextra")
#install.packages("knitr")
#install.packages("GGally")
#install.packages("viridis")
#install.packages("hrbrthemes")

library(rlang)
library(htmltools)
library(Seurat)
library(dplyr)
library(ggplot2)
library(knitr)
library(utils)
library(reshape2)
library(mclust, quietly = TRUE)
library(Matrix)
library(installr)
library(factoextra)
library(knitr)
library(tidyr)
library(hrbrthemes)
library(GGally)
library(viridis)


# Load dataset ----

load_dataset <- function(zip_file, csv_file) {
  
  # This function loads the dataset from a csv file.
  
  data <- read.csv(unz(zip_file, csv_file), sep = ",")
  print(head(data))
  
  return(data)
}


# Data exploration ----

# 
melt_dataset <- function(data) {
  
  # This function melts a dataset data frame into long format.
  
  melted_data <- melt(data, value.name = 'expression')
  melted_data <- melted_data %>%
    rename(cells = "X", genes = "variable")
  
  return(melted_data)
  
}

#
get_library_size <- function(melted_data){
  
  # This function returns the library size for all cells in the dataset.
  
  library_size <- melted_data %>%
    group_by(cells) %>%
    summarise(library_size = sum(expression))
  
  print(library_size)
  
  return(library_size)
  
}
  
#
plot_GeneExpression_heatmap <- function(melted_data) {
  
  # This function creates a gene expression heatmap plot
  
  heatmap_plot <- ggplot(melted_data, aes(x = genes, y = cells, fill = expression)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +  
    labs(x = "Genes", y = "Cells", title = "Gene Expression Heatmap") +
    theme(axis.text = element_blank())
  
  return(heatmap_plot)
  
}



# Data preprocessing ----

# Creation of Seurat object from the expression matrix data

#
build_seurat_object <- function(data, project_name){
  
  # This function builds a Seurat object from a scRNA seq expression matrix dataset.
  # Arguments:
  #     - data (matrix): gene expression matrix
  #     - project_name (str): desired naming of the project
  
  # Remove cell ID column from the expression matrix and make it index
  rownames(data) <- as.character(data[, 1])
  data <- data[, -1]
  
  # Transpose data matrix to be in correct format to be fit in a Seurat object
  data_transpose <- t(data)
  
  # Initialize the Seurat object
  data_seurat <- CreateSeuratObject(counts = data_transpose, project = project_name)
  
  return(data_seurat)
  
}


# Quality Control ----

# Step to identify low quality cells from the downstream analysis.

#
find_mitochondrial_genes <- function(seurat_object){
  
  # This function finds the mitochondrial genes in a Seurat object,
  # based on the "MT" label they have.
  seurat_object[["percent_mito"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  return(seurat_object)
  
}


#
plot_Features_violin <- function(seurat_object){
  
  # This function creates a violin plot of the specified features.
  
  vln_plot <- VlnPlot(seurat_object, features = c("nCount_RNA", "nFeature_RNA", "percent_mito"), ncol = 3)
  
  return(vln_plot)
  
}


#
plot_Features_scatter <- function(seurat_object){
  
  # This function creates a scatter plot of the specified features
  
  scatter_plot <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = 'lm')
  
  return(scatter_plot)
  
}


# Filtering ----

# Step to remove low quality cells from the downstream analysis.

# 
filter_cells <- function(seurat_object, nFeature_RNA_min, nFeature_RNA_max, percent_mito_max){
  
  # This function performs filtering to remove low quality cells, according to 
  # the specified filtering criteria.
  
  filtered_seurat <- subset(seurat_object, 
                            subset = nFeature_RNA > nFeature_RNA_min & 
                              nFeature_RNA < nFeature_RNA_max & 
                              percent_mito < percent_mito_max)
  
  return(filtered_seurat)
  
}


# Data Normalization ----

# The function NormalizeData from the Seurat package is called for this purpose. 


# Highly variable gene identification ----

#
FindPlot_variable_genes <- function(seurat_object){
  
  # This function finds the variable genes of the Seurat object and 
  # creates a plot showing the 10 most variable genes.
  
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
  
  top10_genes <- head(VariableFeatures(seurat_object), 10)
  
  plot_variable_genes <- VariableFeaturePlot(seurat_object)
  plot_variable_genes <- LabelPoints(plot = plot_variable_genes, points = top10_genes, repel = TRUE)
  
  print(plot_variable_genes)
  
  return(seurat_object)
  
}






















