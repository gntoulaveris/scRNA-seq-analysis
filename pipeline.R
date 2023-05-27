
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
create_GeneExpression_heatmap <- function(melted_data) {
  
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

# Step to remove high quality cells from the downstream analysis.




























