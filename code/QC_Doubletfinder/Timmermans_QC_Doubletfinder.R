library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)

## Timmermans path
### path.variables ###
timmermans.path <- "C:/BGI/root_integration/data/Ara/Timmermans_data/Timmer_Ara.rds"

## read in data 
### load data ###
Timmermans <- readRDS(timmermans.path) # 5826 cells

## add study information
Timmermans$Study <- 'Timmermans'
Timmermans$Sample <- Timmermans$orig.ident

##提取线粒体基因
Timmermans[["percent.mt"]] <- PercentageFeatureSet(Timmermans, pattern='^ATMG')
Timmermans[["percent.ct"]] <- PercentageFeatureSet(Timmermans, pattern='^ATCG')

violin <- VlnPlot(Timmermans,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
        pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
        ncol = 4)

ggsave("C:/BGI/root_integration/figures/Ara/QC/Timmermans-vlnplot-before-qc.pdf", plot = violin, width = 15, height = 6) 
ggsave("QC/vlnplot-before-qc.png", plot = violin, width = 15, height = 6) 
plot1 <- FeatureScatter(wang, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wang, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("QC/pearplot-before-qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot-before-qc.png", plot = pearplot, width = 12, height = 5)

##set QC critiria  5721 cells
Timmermans <- subset(Timmermans,subset=nFeature_RNA > 500 & nFeature_RNA < 12000
               & percent.mt < 0.5 & percent.ct < 1)

dim(Timmermans)

## post-QC
violin <- VlnPlot(Timmermans,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
        pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
        ncol = 4)
ggsave("C:/BGI/root_integration/figures/Ara/QC/Timmermans-vlnplot-after-qc.pdf", 
       plot = violin, width = 15, height = 6) 

## Doutbletfinder
library(DoubletFinder)

#### Normalization ####
Timmermans <- NormalizeData(Timmermans, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection: HVGs ### -----method 2
Timmermans <- FindVariableFeatures(Timmermans, selection.method = "vst", nfeatures = 2000)

### scale data ###
Timmermans <- ScaleData(Timmermans) # uncorrected

### dimensionality reduction: PCA ###
Timmermans <- RunPCA(Timmermans)

## elbow plot - CKN ##
pca.elbow.plot <- ElbowPlot(Timmermans, ndims = 50, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot_ckn.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()


## run umap
Timmermans <- RunUMAP(Timmermans, reduction = "pca", dims = 1:50)

## FindNeighbors
Timmermans <- FindNeighbors(Timmermans, reduction = "pca", dims = 1:50)

## FindClusters
Timmermans <- FindClusters(Timmermans, resolution = 0.8)

### splity sample
Shr <- subset(Timmermans, subset = Sample == 'Timmermans_Shr')
Wt <- subset(Timmermans, subset = Sample == 'Timmermans_wt')

## Shr
sweep.res.list <- paramSweep_v3(Shr, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(Shr)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(Shr$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(Shr)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
Shr <- doubletFinder_v3(Shr, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(Shr, reduction = "umap", group.by = "DF.classifications_0.25_0.3_7")
ggsave('Shr_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ---- 1018 cells remained
Shr <- subset(Shr, subset = DF.classifications_0.25_0.3_7 == 'Singlet')


## Wt
sweep.res.list <- paramSweep_v3(Wt, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(Wt)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(Wt$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(Wt)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
Wt <- doubletFinder_v3(Wt, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(Wt, reduction = "umap", group.by = "DF.classifications_0.25_0.26_163")
ggsave('Wt_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ---- 5476 cells remained
Wt <- subset(Wt, subset = DF.classifications_0.25_0.26_163 == 'Singlet')


###merge them together
Timmermans <- merge(Shr, Wt)

### save data ###   5551 cells final
saveRDS(Timmermans, "C:/BGI/root_integration/data/Ara/Timmermans_data/Timmermans_doublet_filtered.rds")

