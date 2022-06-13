library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)

## Wang path
Wang.path <- "C:/BGI/root_integration/data/Ara/Wang_data/Wang_Ara.rds"

# load Wang's data
Wang <- readRDS(Wang.path)

## add study information
Wang$Study <- 'Wang'

##提取线粒体基因
Wang[["percent.mt"]] <- PercentageFeatureSet(Wang, pattern='^ATMG')
Wang[["percent.ct"]] <- PercentageFeatureSet(Wang, pattern='^ATCG')
VlnPlot(Wang,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
        pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
        ncol = 4)

ggsave("QC/vlnplot-before-qc.pdf", plot = violin, width = 15, height = 6) 
ggsave("QC/vlnplot-before-qc.png", plot = violin, width = 15, height = 6) 
plot1 <- FeatureScatter(wang, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wang, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("QC/pearplot-before-qc.pdf", plot = pearplot, width = 12, height = 5) 
ggsave("QC/pearplot-before-qc.png", plot = pearplot, width = 12, height = 5)

##set QC critiria
Wang <- subset(Wang,subset=nFeature_RNA > 500 & nFeature_RNA < 5000
             & percent.mt < 1 & percent.ct < 3)

dim(Wang)

## post-QC
VlnPlot(Wang,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
        pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
        ncol = 4)
ggsave("QC/vlnplot-after-qc.pdf", plot = violin, width = 15, height = 6) 
ggsave("QC/vlnplot-after-qc.png", plot = violin, width = 15, height = 6)
## 基因表达量标准化
## 它的作用是让测序数据量不同的细胞的基因表达量具有可比性。计算公式如下：
## 标准化后基因表达量 = log1p（10000*基因counts/细胞总counts）
wang <- NormalizeData(wang, normalization.method = "LogNormalize", scale.factor = 10000)

## Doutbletfinder
library(DoubletFinder)

#### Normalization ####
Wang <- NormalizeData(Wang, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection: HVGs ### -----method 2
Wang <- FindVariableFeatures(Wang,mean.cutoff=c(0.0125,3),dispersion.cutoff =c(1.5,Inf))

### scale data ###
Wang <- ScaleData(Wang) # uncorrected

### dimensionality reduction: PCA ###
Wang <- RunPCA(Wang,npcs = 100)

## elbow plot - CKN ##
pca.elbow.plot <- ElbowPlot(Wang, ndims = 100, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot_ckn.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()


## run umap
Wang <- RunUMAP(Wang,n.neighbors = 30,metric = 'correlation', min.dist = 0.3, dims = 1:40)

## FindNeighbors
Wang <- FindNeighbors(Wang, reduction = "pca", dims = 1:100)

## FindClusters
Wang <- FindClusters(Wang, resolution = 0.5)


## Ѱ??????pKֵ CKN
sweep.res.list <- paramSweep_v3(Wang, PCs = 1:100)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(Wang)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(Wang$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(Wang)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
Wang <- doubletFinder_v3(Wang, PCs = 1:100, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(Wang, reduction = "umap", group.by = "DF.classifications_0.25_0.28_1153")
ggsave('Wang_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----11337 cells remained
Wang <- subset(Wang, subset = DF.classifications_0.25_0.28_1153 == 'Singlet')

### save data ###
saveRDS(Wang, "C:/BGI/root_integration/data/Ara/Wang_data/Wang_doublet_filtered.rds")

## re-run finderallmarkers 
## compute cluster marker genes ###
CK.markers <- FindAllMarkers(CK.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
CK.top10.markers <- CK.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(CK.markers, file = "C:/BGI/mianhua/CK_harmony/dim50_annotation/CK_Marker_genes_afterdoublet.csv")
write.csv(CK.top10.markers, file = "C:/BGI/mianhua/CK_harmony/dim50_annotation/top10_marker_genes_afterdoublet.csv")
