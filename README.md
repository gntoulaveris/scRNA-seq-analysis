# scRNA-seq-analysis
This an unsupervised learning project analyzing single cell RNA seq data. The data are model generated expression matrices. Each dataset consists of the expression matrix of 200 cells and 200 genes. The pipeline code was written in R and implemented for each dataset in R markdown files. 

Data preprocessing and quality control were performed using the Seurat package, which is one of the most famous packages for scRNA seq data analysis. Using the same package the data were filtered to remove low-quality cells from consideration and then normalized in log space. The most highly variable genes were discovered and then the data were scaled.

Dimensionality reduction was performed in three different ways: PCA, uMAP, tSNE. Each technique was implemented using the Seurat package. Each technique led to different clustering results.

Clustering analysis for all dimensionality reduction methods was performed with Gaussian Mixture Models (GMM). The optimal GMM model in terms of the number of components and covariance matrix was discovered by using the BIC criterion. For clustering the package mclust was utilized, which is a contributed R package for model-based clustering, classification, and density estimation based on finite normal mixture modelling. It provides functions for parameter estimation via the EM algorithm for normal mixture models with a variety of covariance structures, and functions for simulation from these models. It also utilized the BIC criterion for model selection.

Visualization of the clustering results was performed either by utilizing premade mclust functions or using ggplot2, a famous R package for data visualization.


## Files
* scRNA_pipeline: R source file containing the functions needed to fully analyze and visualize a scRNA seq dataset.
* dataset#_analysis: the R markdown files containing the analysis of each dataset, by imlementing the scRNA_pipeline.
* Assignment_3_2023_datasets.zip: the datasets that were analyzed.


## Dependencies
In order to successfully implement the analysis of a scRNA seq dataset with the scRNA_pipeline code the following dependencies must be met. For ease of implementation a dataset#_analysis markdown file can be utilized with the user only having to pass on the following dependencies.

* Pipeline R file, implementation (analysis) markdown file(s) and datasets' folder need to be in the same working directory.
* The working directory must be set by the user before the code is run.
* The datasets need to be inside a zip folder and that zip file needs to be inside the working directory.
* The datasets need to be in .csv format.

* The scRNA_pipeline needs to be called inside the R markdown file (command: source("scRNA_pipeline.R"))
* The user needs to specify as variables the zip file name and the dataset .csv name and pass them as arguments to the load_dataset function.
* The user needs to call the get_dataset_name function and pass as argument the previously created variable for the dataset .csv file. The output of the function should be stored in a variable called dataset_name.

After the above steps the rest of the markdown file runs automatically. If an example markdown file is not used, then the user can call the functions stored in scRNA_pipeline in order of appearence and simply pass to them as arguments the above listed variables.
