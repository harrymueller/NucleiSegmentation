print("Installing R packages...")

repo="https://cran.curtin.edu.au/"

# pkg managers
install.packages("devtools", repos=repo) 
install.packages("BiocManager", repos=repo)

# Seurat
install.packages("Seurat", version = "4.0.0", repos = repo)

# singleR v1.8.0 (report was v1.6.1)
BiocManager::install("SingleR") 

# scibet v1.0
devtools::install_github("PaulingLiu/scibet")

# clusterProfiler v4.3.0 (report was v4.0.0)
BiocManager::install("clusterProfiler")

# monocle3 v1.0.0
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

install.packages("markdown")

# gem2rds.R from BGI
install.packages("raster")
install.packages("rgdal")
install.packages("argparser")

# markers
install.packages("rJava")
install.packages("xlsx")
install.packages('dplyr', repos = 'https://cloud.r-project.org')

# cell annotations
BiocManager::install("celldex")
BiocManager::install("scRNAseq")

print("Installed R packages")