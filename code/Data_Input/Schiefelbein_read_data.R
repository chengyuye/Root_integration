library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(umap)
library(reticulate)

### load Schiefelbein data ###
Schiefelbein.data.dir <- "C:/BGI/root_integration/data/Ara/Schiefelbein_data/"
Schiefelbein.data.samples <- c("Sample_gl2","Sample_rhd6", "Sample_WT-WERGFP",
                               'Sample_WT-WERGFP_2', 'Sample_WT-WERGFP_3')

Sample.Id <- c('g12','rhd6','WT1','WT2','WT3')


samples.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/samples/"
merge.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/raredon_lung.rds"
raredon_lung.path <- "/home/s1987963/ds_group/Niklas/raredon_lung/"

# load lung samples separately #
for(i in 1:length(Schiefelbein.data.samples)){
  assign(paste0(Schiefelbein.data.samples[i],'.counts'), Read10X(data.dir = paste0(Schiefelbein.data.dir, Schiefelbein.data.samples[i])))
}

# initialize Seurat objects #
for(i in 1:length(Schiefelbein.data.samples)){
  data <- Schiefelbein.data.samples[i]
  count_data <- paste0(Schiefelbein.data.samples[i],".counts")
  samp.ID <- Sample.Id[i]
  eval(parse(text=paste0(data," <- CreateSeuratObject(counts = ", count_data,", project = \"Schiefelbein\", min.cells=3, min.features=200)")))
  eval(parse(text=paste0("rm(", count_data, ")")))
  
  # add metadata: author_dataset
  eval(parse(text=paste0(data,"$Study <- \"Schiefelbein\"")))
  
  # add metadata: Sample ID
  eval(parse(text=paste0(data,"$Sample <- \"", samp.ID, "\"")))
  
  # add meta data: calculate missing metrics
  eval(parse(text=paste0(data,"[[\"percent.mt\"]] <- PercentageFeatureSet(", data, ", pattern = \"^MT-\")")))
}

### load data ###
# load samples samples separately #
Schiefelbein.counts <- read.table(file= Schiefelbein.path, sep="\t", header = T,fill = TRUE)
rownames(Schiefelbein.counts) <- Schiefelbein.counts$X
Schiefelbein.counts <- Schiefelbein.counts[,-1]

### create SeuratObject   #### 11030 cells
Schiefelbein <- CreateSeuratObject(counts = Schiefelbein.counts, 
                                   min.cells = 3, min.features=200, 
                                   project = 'Schiefelbein', names.delim = "_")

Schiefelbein$Sample <- Schiefelbein$orig.ident
Schiefelbein$Study <- 'Schiefelbein'

# Save the SeuratObject as rds data
saveRDS(Schiefelbein, 'C:/BGI/root_integration/data/Ara/Schiefelbein_data/Schiefelbein_Ara.rds')
