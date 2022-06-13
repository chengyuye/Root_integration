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
## for Wang
Wang.path <- "C:/BGI/root_integration/data/Ara/Wang_data/"

# load Wang's data
Wang.counts <- Read10X(data.dir = Wang.path)


# create SeuratObject ####14461 cells; 22234 genes
Wang <- CreateSeuratObject(counts = Wang.counts, min.cells = 3, min.features=200,
                                 project = 'Wang')


# Save the SeuratObject as rds data
saveRDS(Wang, file = "C:/BGI/root_integration/data/Ara/Wang_data/Wang_Ara.rds")

