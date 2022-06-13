library(Seurat)
library(dplyr)
library(umap)
library(reticulate)
library(ggplot2)
library(patchwork)
library(celda)

### read raw data ###
Rybel.path <- "C:/BGI/root_integration/data/Ara/Rybel/GSE141730_RAW/"
Sample.Id <- c('GSM4212550_BDR3','GSM4212551_BDR4','GSM4212552_BDR5')
simple_sampleO_id <- c('BDR3', 'BDR4','BDR5')

# load raw lung samples separately #
Rybel.counts <- Read10X_h5(Rybel.path, use.names = T)

Rybel <- CreateSeuratObject(Rybel.counts, min.cells = 3, project = "Rybel", min.features=200)

for(i in 1:length(Sample.Id)){
  assign(paste0("Rybel.", Sample.Id[i],".raw.data"), Read10X_h5(paste0(Rybel.path, Sample.Id[i], "_filtered_gene_bc_matrices.h5")))
}

# initialize seurat objects #
for(i in 1:length(Sample.Id)){
  
  # name of seurat objects
  data <- paste0("Rybel.", Sample.Id[i],".raw")
  
  # name of count data
  count.data <- paste0("Rybel.", Sample.Id[i],".raw.data")
  
  # initialize seurat object
  eval(parse(text=paste0(data," <- CreateSeuratObject(counts = ", count.data,", project = \"Rybel\", min.cells=3, min.features=200)")))
  eval(parse(text=paste0("rm(", count.data, ")")))
  
  # add metadata: study
  eval(parse(text=paste0(data,"$Study <- \"Rybel\"")))
  
  # add metadata: sample ID
  eval(parse(text=paste0(data,"$Sample <- \"", simple_sampleO_id[i], "\"")))
  
  # add meta data: calculate missing metrics
  eval(parse(text=paste0(data,"[[\"percent.mt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^ATMG\")")))
  
  # add meta data: calculate missing metrics
  eval(parse(text=paste0(data,"[[\"percent.ct\"]] <- PercentageFeatureSet(", data, ", pattern = \"^ATCG\")")))
}

### merge raw data ###  ---- 5145 cells
## merge Seurat objects
Rybel.raw <- merge(Rybel.GSM4212550_BDR3.raw, 
                   c(Rybel.GSM4212551_BDR4.raw, Rybel.GSM4212552_BDR5.raw))

### save merged raw data ###
saveRDS(Rybel.raw, file = "C:/BGI/root_integration/data/Ara/Rybel/Rybel_Ara.rds")

