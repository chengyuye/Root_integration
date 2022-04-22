library(Seurat)
library(dplyr)
library(patchwork)
library(tidyverse)
library(Seurat)
library(kableExtra)
library(devtools)
#install_github('welch-lab/liger')
library(rliger)
library(data.table)
library(ggsci)


# read in all the root data 
# Ara
Ara_root <-  readRDS(file = "C:/BGI/root_integration/data/ninanjie_result.rds")
head(rownames(Ara_root))

# Lot
Lj_root = readRDS("C:/BGI/root_integration/data/Lj_root.rds")
head(rownames(Lj_root))

# rice
rice_root = readRDS("C:/BGI/root_integration/data/shuidao_result.rds")
head(rownames(rice_root))

# overview of samples
datasets <- c("Ara_root", "Lj_root", "rice_root")

data.frame(
  "Species" = c("Ara","Lot","Rice"),
  "Clusters" = map_int(datasets, ~ length(unique(get(.x)@active.ident)))
) %>%
  kable(caption = "Overview") %>%
  kable_styling(latex_options = "hold_position")



plot_settings <- list(
  theme_void(),
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10))
)
cowplot::plot_grid(plotlist = list(
  DimPlot(Ara_root, label = TRUE) +
    #scale_color_manual(values = Ara_root@misc$cluster_colors) +
    ggtitle("Ara") +
    plot_settings,
  DimPlot(Lj_root, label = TRUE) +
    ggtitle("Lot") +
    #scale_color_manual(values = nhp@misc$cluster_colors) +
    plot_settings,
  DimPlot(rice_root, label = TRUE) +
    ggtitle("Rice") +
    #scale_color_manual(values = nhp@misc$cluster_colors) +
    plot_settings
))

# keep metadata for future use
metadata <- map(
  datasets,
  ~ data.frame(
    "cell" = rownames(get(.x)@meta.data),
    "Cluster" = get(.x)@active.ident)
)
names(metadata) <- datasets



# Preprocessing

#Before I integrate the data, 
#I'm going to keep only the genes that are shared in all datasets.

# isolate just the raw count matrix files by species
mtx <- map(datasets, ~ get(.x)@assays$RNA@counts)
names(mtx) <- datasets


print("Reading in ortholog info from Ensembl ...")
# read in ortholog info and only keep orthologs with 1:1 mapping
## read in orthogroups result
orthogroups <- read.csv("C:/BGI/root_integration/data/Orthogroup_result.csv")
orthogroups <- orthogroups %>% filter(X4 != "")
#orthogroups<- orthogroups[orthogroups[,5]!="",] (Alternative way)


# get unique genes by species
unique_ara <- orthogroups$X1
unique_lot <- orthogroups$X3
unique_lot <- gsub('_',"-",unique_lot)
unique_rice <- orthogroups$X4
unique_rice <- gsub('_',"-",unique_rice)
gene_num <- length(unique_ara) %>% as.numeric()

## replace the lot name to ara name
## replace the rice name to ara name
## this step also keeps their other gene names 
## rather than only keep those common genes
for(i in 1:gene_num){
  rownames(mtx$Lj_root) <- gsub(unique_lot[i], unique_ara[i],rownames(mtx$Lj_root))
  rownames(mtx$rice_root) <- gsub(unique_rice[i], unique_ara[i],rownames(mtx$rice_root))
}

# rename the cell to avoid erros
colnames(mtx$Ara_root) <- paste("Ara", colnames(mtx$Ara_root), sep = "_")
colnames(mtx$Lj_root) <- paste("LOT", colnames(mtx$Lj_root), sep = "_")
colnames(mtx$rice_root) <- paste("Rice", colnames(mtx$rice_root), sep = "_")

###using liger for integration
species.liger <- createLiger(list(ara = mtx$Ara_root, 
                                  lot = mtx$Lj_root,
                                  rice = mtx$rice_root))

###Normalize the datasets. 
#The normalization is applied to the datasets in their entirety.
species.liger <- normalize(species.liger)

#Select shared, homologous genes between the two species, 
#as well as unshared, non-homologous genes
species.liger <- selectGenes(species.liger, var.thres= 0.2,unshared = TRUE, 
                             unshared.datasets = list(1,2,3), unshared.thresh= c(0.3,0.2,0.4),
                             do.plot = T)

#scale, but do not center, the data.
species.liger <- scaleNotCenter(species.liger)

#Selecting an appropriate K-Value
seed = sample(1:200, 3, replace = FALSE)
alignment_score = list()
for (iter in seed){
  K_values = c(10, 20, 30, 40, 50)
  for (i in K_values){
    osm.liger <- optimizeALS(osm.liger,  lambda = 5 , use.unshared = TRUE, max.iters = 30, thresh=1e-10, k =i, rand.seed = iter)
    osm.liger <- quantile_norm(osm.liger, ref_dataset = "rna")
    osm.liger <- louvainCluster(osm.liger)
    new_alignment = calcAlignment(osm.liger)
    names(new_alignment) = paste0("Seed:", iter, "_K:",i)
    alignment_score = append(alignment_score, new_alignment)
  }
}

####Joint Matrix Factorization
species.liger <- optimizeALS(species.liger,  lambda = 10, 
                             use.unshared = TRUE, thresh=1e-10, k =20)

#Quantile Normalization and Joint Clustering
species.liger <- quantile_norm(species.liger, ref_dataset = "ara")
species.liger <- louvainCluster(species.liger,resolution = 0.25)

##Visualizations and Downstream processing
#UMAP
species.liger <- runUMAP(species.liger)

umap_plots <- plotByDatasetAndCluster(species.liger, 
                                      axis.labels = c("UMAP1","UMAP2"), 
                                      return.plots = TRUE) 
pdf(file = 'C:/BGI/root_integration/Ara_lot_rice_integrated_atlas.pdf',
    width = 7, height = 5.5)
umap_plots[[1]] 
dev.off()

##TSNE
species.liger <- runTSNE(species.liger)

tsne_plots <- plotByDatasetAndCluster(species.liger, 
                                      axis.labels = c("TSNE1","TSNE2"), 
                                      return.plots = TRUE) 
pdf(file = 'C:/BGI/root_integration/Ara_lot_rice_integrated_atlas_TSNE.pdf',
    width = 7, height = 5.5)
tsne_plots[[1]] 
dev.off()

# only keep genes that are in the dataset
unique_ara <- unique_ara[unique_ara %in% rownames(Ara_root)]
unique_lot <- unique_lot[unique_lot %in% rownames(Lj_root)]

# keep orthologs
orthologs <- orthologs %>% filter(
  `Gene name` %in% unique_mouse & `Gene name_1` %in% unique_nhp
)
Ara_filter <- subset(x = Ara_root, features = orthofinder$X1)
length(rownames(Ara_filter))


orthofinder$X3 = gsub('_',"-",orthofinder$X3)
orthofinder$X4 = gsub('_',"-",orthofinder$X4)
Lj_filter <- subset(x = Lj_root, features = orthofinder$X3)
length(rownames(Lj_filter))

result=orthofinder[0,]
for (i in 1:length(Lj_filter@assays$RNA@counts@Dimnames[[1]])){
  results = subset(orthofinder,X3==Lj_filter@assays$RNA@counts@Dimnames[[1]][i])
  result=rbind(result,results)
}

Lj_filter@assays$RNA@counts@Dimnames[[1]] <- result$X1
Lj_filter@assays$RNA@data@Dimnames[[1]] <- result$X1
Lj_filter@assays$RNA@meta.features <- data.frame (row.names = rownames (Lj_filter@assays$RNA))

# ˮ????????
rice_root = readRDS(file = "D:/singlecell-gu/Nip_data/result/result.rds")
rice_filter <- subset(x = rice_root, features = orthofinder$X4)
length(rownames(rice_filter))

result=orthofinder[0,]
for (i in 1:length(rice_filter@assays$RNA@counts@Dimnames[[1]])){
  results = subset(orthofinder,X4==rice_filter@assays$RNA@counts@Dimnames[[1]][i])
  result=rbind(result,results)
}

rice_filter@assays$RNA@counts@Dimnames[[1]] <- result$X1
rice_filter@assays$RNA@data@Dimnames[[1]] <- result$X1
rice_filter@assays$RNA@meta.features <- data.frame (row.names = rownames (rice_filter@assays$RNA))



# RenameGenesSeurat <- function(obj, newnames ) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
#   print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
#   RNA <- obj@assays$RNA
#   
#   if (nrow(RNA) == length(newnames)) {
#     if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
#     if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
#   } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
#   obj@assays$RNA <- RNA
#   return(obj)
# }
# RenameGenesSeurat(obj = rice_filter, newnames = result$Arabidopsis_thaliana.TAIR10.pep.all)


merge1 = SplitObject(Ara_filter, split.by = "orig.ident")
merge2 = SplitObject(Lj_filter, split.by = "orig.ident")
merge3 = SplitObject(rice_filter, split.by = "orig.ident")
merge.list <- c(merge1$tair_root, merge2$`Cro-root`,merge3$Nip_root)


for (i in 1:length(merge.list)) {
  merge.list[[i]] <- NormalizeData(merge.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  merge.list[[i]] <- FindVariableFeatures(merge.list[[i]], selection.method = "vst",
                                             nfeatures = 1000, verbose = FALSE)
}

features = SelectIntegrationFeatures(object.list = merge.list,nfeatures=500)
#?FindIntegrationAnchors
pancreas.anchors <- FindIntegrationAnchors(object.list = merge.list , anchor.features = features, dims = 1:30 )  ## dims = 1:20,k.anchor = 5,k.filter = 30
#?IntegrateData
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:20)  ## dims = 1:20

DefaultAssay(pancreas.integrated) <- "integrated"        
pancreas.integrated <- ScaleData(pancreas.integrated,features=VariableFeatures(pancreas.integrated)) 
pancreas.integrated <- RunPCA(pancreas.integrated, features = VariableFeatures(object = pancreas.integrated))
ElbowPlot(pancreas.integrated)
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:50)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:50)

# pancreas.integrated <- RunTSNE(pancreas.integrated, dims = 1:16)
DimPlot(pancreas.integrated, reduction = "umap",label = F,group.by="orig.ident")
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident") 
DimPlot(pancreas.integrated, reduction = "umap") 

Idents(pancreas.integrated) = "celltype"
# pancreas.integrated@meta.data[pancreas.integrated$celltype == "Metaxylem",]$celltype = "xylem"
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident") 
DimPlot(pancreas.integrated, reduction = "umap",label = T, split.by  = "orig.ident") 

DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = c("Epidermis","Epidermis (near root hair)")))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = "xylem"))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = c("Root_cap","Root_hair","Cap&Root hair")))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = "Stele"))
DimPlot(pancreas.integrated, reduction = "umap", split.by  = "orig.ident",
        cells.highlight = WhichCells(pancreas.integrated, idents = c("Root_cap")))



saveRDS(pancreas.integrated, file = "D:/singlecell-liang/write/cro_ara_merge_result/ara_Lj_merge_result.rds")


### read in the merged data
root_merge <-  readRDS(file = "D:/singlecell-liang/write/cro_ara_merge_result/ara_Lj_merge_result.rds")

pbmc.markers <- FindAllMarkers(pancreas.integrated, 
                               only.pos = TRUE, 
                               min.pct = 0.5, 
                               min.diff.pct = 0.3,
                               logfc.threshold = 0.25) 


DefaultAssay(pancreas.integrated) <- "RNA" 
FeaturePlot(pancreas.integrated, features = c("AT5G05500", "AT1G30870", "AT1G12560"),max.cutoff = 'q80',keep.scale = "all")

################ ??Ⱥ??????????ͼ
pancreas.integrated = readRDS(file = "D:/singlecell-liang/write/cro_ara_merge_result/ara_Lj_merge_result.rds")
orthofinder = read.csv('D:/singlecell-liang/write/Orthogroup_result.csv')
orthofinder = orthofinder[orthofinder[,5]!="",]
list = list(orthofinder[1:50,"X1"])
object <- AddModuleScore(object = pancreas.integrated, features = list, name = "gene_type")
FeaturePlot(object = object, features = "gene_type1")
### top 20% of the cells with high expression levels           
object <- MetaFeature(object = pancreas.integrated, features = orthofinder[1:50,"X1"], meta.name = "Aggregate_Feature") 
data = object@meta.data$Aggregate_Feature
threshold = sort(data, decreasing = T)[0.2*length(data)]
object@meta.data$Aggregate_Feature = replace(data,data<threshold,0)
FeaturePlot(object = object, features = "Aggregate_Feature",max.cutoff = "q80")
################
