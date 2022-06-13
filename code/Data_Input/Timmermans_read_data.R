library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(Seurat)
library(kableExtra)
library(harmony)
library(rliger)
library(data.table)
library(ggsci)


### path.variables ###
## for Timmer
Timmer.path1 <- "C:/BGI/root_integration/data/Ara/Timmermans_data/GSE123818_Root_single_cell_shr_datamatrix.csv"
Timmer.path2 <- "C:/BGI/root_integration/data/Ara/Timmermans_data/GSE123818_Root_single_cell_wt_datamatrix.csv"

# load sepsis data
Timmer.counts1 <- read.csv(file= Timmer.path1, sep=",", header = T, row.names = 1)
Timmer.counts2 <- read.csv(file= Timmer.path2, sep=",", header = T, row.names = 1)

# csv合并成数据框
Timmer <-  data.frame(Timmer.counts1,Timmer.counts2)

# 数据框转换成稀疏矩阵matrix
Timmer <- as.sparse(Timmer)

# create SeuratObject
Timmer_shr <- CreateSeuratObject(counts = Timmer.counts1, min.cells = 3, min.features=200,
                             project = 'Timmermans_Shr')
Timmer_wt <- CreateSeuratObject(counts = Timmer.counts2, min.cells = 3, min.features=200,
                                 project = 'Timmermans_wt')

##merge the two dataset
Timmer <- merge(Timmer_shr, Timmer_wt)

# Save the SeuratObject as rds data
saveRDS(Timmer, file = "C:/BGI/root_integration/data/Ara/Timmermans_data/Timmer_Ara.rds")

