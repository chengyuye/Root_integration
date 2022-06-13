library(Seurat)
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(umap)
library(reticulate)
library(RColorBrewer)
library(data.table)
library(monocle)

### path.variables ###
cuperus.path <- "C:/BGI/root_integration/data/Ara/Cuperus/"
metadata.path <- "C:/BGI/root_integration/data/Ara/Cuperus/GSE121619_genes.tsv.gz"

# patient IDs
patient.ID <- c("GSM3823939_control", "GSM3823940_control", "GSM3823941_control",
                "GSM3823942_diabetes", "GSM3823943_diabetes", "GSM3823944_diabetes")

### load count matrices ###
cuperus <- readRDS("C:/BGI/root_integration/data/Ara/Cuperus/GSE121619_Control_Heatshock_cds.rds")
### control data
control <- cuperus$control
count.mat <- Biobase::exprs(control)
meta.df <- as.data.frame(Biobase::pData(control))

control <- CreateSeuratObject(counts = count.mat,
                              project = "Cuperus")
                                #assay = "RNA",
                                #meta.data = meta.df)

control$Study <- 'Cuperus'
control$Sample <- 'Control'

###heatshock data
heatshock <- cuperus$heatshock
count.mat <- Biobase::exprs(heatshock)
meta.df <- as.data.frame(Biobase::pData(heatshock))

heatshock <- CreateSeuratObject(counts = count.mat,
                                project = "Cuperus")

heatshock$Study <- 'Cuperus'
heatshock$Sample <- 'Heatshock'

###merge them together 2085 cells
cuperus <- merge(control, heatshock)

### save Seurat Object ###
saveRDS(cuperus, "C:/BGI/root_integration/data/Ara/Cuperus/Cuperus_Ara.rds")


