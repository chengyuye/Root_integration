## library
library(DoubletFinder)

### path.variables ###
Rybel.path <- "C:/BGI/root_integration/data/Ara/Rybel/Rybel_Ara.rds"

Sample.Id <- c('BDR3','BDR4','BDR5')

## read in data 
### load data ###
Rybel <- readRDS(Rybel.path) # 5145 cells

##QC
violin <- VlnPlot(Rybel,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
                  pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
                  ncol = 4)

ggsave("C:/BGI/root_integration/figures/Ara/QC/Rybel-vlnplot-before-qc.pdf", plot = violin, width = 15, height = 6) 

plot1 <- FeatureScatter(Rybel, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Rybel, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("C:/BGI/root_integration/figures/Ara/QC/Rybel_pearplot-before-qc.pdf", plot = pearplot, width = 12, height = 5) 

##set QC critiria    5133 cells remained
Rybel <- subset(Rybel, subset=nFeature_RNA > 500 & nFeature_RNA < 14000
                       & percent.mt < 0.02 )

## post-QC
violin <- VlnPlot(Rybel,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
                  pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
                  ncol = 4)
ggsave("C:/BGI/root_integration/figures/Ara/QC/Rybel-vlnplot-after-qc.pdf", 
       plot = violin, width = 15, height = 6) 

plot1 <- FeatureScatter(Rybel, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Rybel, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("C:/BGI/root_integration/figures/Ara/QC/Schiefelbein_pearplot-after-qc.pdf", plot = pearplot, width = 12, height = 5) 

#### Normalization ####
Rybel <- NormalizeData(Rybel, normalization.method = "LogNormalize", scale.factor = 10000)


### feature selection: HVGs ### -----method 2
Rybel <- FindVariableFeatures(Rybel, selection.method = "vst", nfeatures = 2000)


### scale data ###
Rybel <- ScaleData(Rybel) # uncorrected

### dimensionality reduction: PCA ###
Rybel <- RunPCA(Rybel)


## elbow plot - Rybel ##
pca.elbow.plot <- ElbowPlot(Rybel, ndims = 50, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot_ckn.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## run umap
Rybel <- RunUMAP(Rybel, reduction = "pca", dims = 1:50)

## FindNeighbors
Rybel <- FindNeighbors(Rybel, reduction = "pca", dims = 1:50)

## FindClusters
Rybel <- FindClusters(Rybel, resolution = 0.5)

##split samples and run doubletfiner for each samples 
Rybel.list <- SplitObject(Rybel, split.by = "Sample")
BDR3 <- Rybel.list[[1]]
BDR4 <- Rybel.list[[2]]
BDR5 <- Rybel.list[[3]]

## BDR3 pKֵ CKN
sweep.res.list <- paramSweep_v3(BDR3, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(BDR3)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(BDR3$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(BDR3)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
BDR3 <- doubletFinder_v3(BDR3, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(BDR3, reduction = "umap", group.by = "DF.classifications_0.25_0.3_29")
ggsave('BDR3_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----1942 cells remained
BDR3 <- subset(BDR3, subset = DF.classifications_0.25_0.3_29 == 'Singlet')
dim(BDR3)

##  BDR4  pKֵ CKN
sweep.res.list <- paramSweep_v3(BDR4, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(BDR4)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(BDR4$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(BDR4)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
BDR4 <- doubletFinder_v3(BDR4, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(BDR4, reduction = "umap", group.by = "DF.classifications_0.25_0.05_16")
ggsave('BDR4_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----1452 cells remained
BDR4 <- subset(BDR4, subset = DF.classifications_0.25_0.05_16 == 'Singlet')
dim(BDR4)

## Ѱ??????pKֵ CKN
sweep.res.list <- paramSweep_v3(BDR5, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(BDR5)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(BDR5$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(BDR5)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
BDR5 <- doubletFinder_v3(BDR5, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(BDR5, reduction = "umap", group.by = "DF.classifications_0.25_0.19_21")
ggsave('BDR5_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----1673 cells remained
BDR5 <- subset(BDR5, subset = DF.classifications_0.25_0.19_21 == 'Singlet')
dim(BDR5)

### merge all the datasets after doubletfinder------ 5067 final cells
Rybel <- merge(BDR3, c(BDR4,BDR5))
dim(Rybel)

### save data ###  
saveRDS(Rybel, "C:/BGI/root_integration/data/Ara/Rybel/Rybel_doublet_filtered.rds")

