## library
library(DoubletFinder)

### path.variables ###
Schiefelbein.path <- "C:/BGI/root_integration/data/Ara/Schiefelbein_data/Schiefelbein_Ara.rds"

Sample.Id <- c('gl2','rhd6','WT1','WT2','WT3')
## read in data 
### load data ###
Schiefelbein <- readRDS(Schiefelbein.path) # 11030 cells


##提取线粒体基因
Schiefelbein[["percent.mt"]] <- PercentageFeatureSet(Schiefelbein, pattern='^ATMG')
Schiefelbein[["percent.ct"]] <- PercentageFeatureSet(Schiefelbein, pattern='^ATCG')

violin <- VlnPlot(Schiefelbein,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
                  pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
                  ncol = 4)

ggsave("C:/BGI/root_integration/figures/Ara/QC/Schiefelbein-vlnplot-before-qc.pdf", plot = violin, width = 15, height = 6) 

plot1 <- FeatureScatter(Schiefelbein, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Schiefelbein, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("C:/BGI/root_integration/figures/Ara/QC/Schiefelbein_pearplot-before-qc.pdf", plot = pearplot, width = 12, height = 5) 

##set QC critiria   10996 cells remained
Schiefelbein <- subset(Schiefelbein,subset=nFeature_RNA > 500 & nFeature_RNA < 11000
                     & percent.mt < 0.03 )

dim(10996)

## post-QC
violin <- VlnPlot(Schiefelbein,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
                  pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
                  ncol = 4)
ggsave("C:/BGI/root_integration/figures/Ara/QC/Schiefelbein-vlnplot-after-qc.pdf", 
       plot = violin, width = 15, height = 6) 

plot1 <- FeatureScatter(Schiefelbein, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Schiefelbein, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("C:/BGI/root_integration/figures/Ara/QC/Schiefelbein_pearplot-after-qc.pdf", plot = pearplot, width = 12, height = 5) 

#### Normalization ####
Schiefelbein <- NormalizeData(Schiefelbein, normalization.method = "LogNormalize", scale.factor = 10000)


### feature selection: HVGs ### -----method 2
Schiefelbein <- FindVariableFeatures(Schiefelbein, selection.method = "vst", nfeatures = 1000)


### scale data ###
Schiefelbein <- ScaleData(Schiefelbein) # uncorrected

### dimensionality reduction: PCA ###
Schiefelbein <- RunPCA(Schiefelbein)


## elbow plot - Schiefelbein ##
pca.elbow.plot <- ElbowPlot(Schiefelbein, ndims = 50, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot_ckn.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## run umap
Schiefelbein <- RunUMAP(Schiefelbein, reduction = "pca", dims = 1:50)

## FindNeighbors
Schiefelbein <- FindNeighbors(Schiefelbein, reduction = "pca", dims = 1:50)

## FindClusters
Schiefelbein <- FindClusters(Schiefelbein, resolution = 0.2)

##split samples and run doubletfiner for each samples 
Schiefelbein.list <- SplitObject(Schiefelbein, split.by = "Sample")
WT1 <- Schiefelbein.list[[1]]
WT2 <- Schiefelbein.list[[2]]
WT3 <- Schiefelbein.list[[3]]
rhd6 <- Schiefelbein.list[[4]]
gl2 <- Schiefelbein.list[[5]]
data.list <- c(WT1,WT2,WT3,rhd6,gl2)

## Ѱ??????pKֵ CKN
sweep.res.list <- paramSweep_v3(WT1, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(WT1)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(WT1$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(WT1)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
WT1 <- doubletFinder_v3(WT1, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                                 nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(WT1, reduction = "umap", group.by = "DF.classifications_0.25_0.3_134")
ggsave('WT1_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----4260 cells remained
WT1 <- subset(WT1, subset = DF.classifications_0.25_0.3_134 == 'Singlet')
dim(WT1)

## Ѱ??????pKֵ CKN
sweep.res.list <- paramSweep_v3(WT2, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(WT2)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(WT2$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(WT2)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
WT2 <- doubletFinder_v3(WT2, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(WT2, reduction = "umap", group.by = "DF.classifications_0.25_0.16_15")
ggsave('WT2_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----1457 cells remained
WT2 <- subset(WT2, subset = DF.classifications_0.25_0.16_15 == 'Singlet')
dim(WT2)

## Ѱ??????pKֵ CKN
sweep.res.list <- paramSweep_v3(WT3, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(WT3)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(WT3$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(WT3)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
WT3 <- doubletFinder_v3(WT3, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(WT3, reduction = "umap", group.by = "DF.classifications_0.25_0.01_19")
ggsave('WT3_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----1621 cells remained
WT3 <- subset(WT3, subset = DF.classifications_0.25_0.01_19 == 'Singlet')

## Ѱ??????pKֵ CKN
sweep.res.list <- paramSweep_v3(rhd6, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(rhd6)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(rhd6$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(rhd6)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
rhd6 <- doubletFinder_v3(rhd6, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(rhd6, reduction = "umap", group.by = "DF.classifications_0.25_0.09_44")
ggsave('rhd6_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----2455 cells remained
rhd6 <- subset(rhd6, subset = DF.classifications_0.25_0.09_44 == 'Singlet')

## Ѱ??????pKֵ CKN
sweep.res.list <- paramSweep_v3(gl2, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(gl2)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(gl2$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(gl2)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
gl2 <- doubletFinder_v3(gl2, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                        nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(gl2, reduction = "umap", group.by = "DF.classifications_0.25_0.1_6")
ggsave('gl2_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----9633 cells remained
gl2 <- subset(gl2, subset = DF.classifications_0.25_0.1_6 == 'Singlet')


### merge all the datasets after doubletfinder------10778 final cells
Schiefelbein <- merge(WT1,c(WT2,WT3, rhd6,gl2))

### save data ###  
saveRDS(Schiefelbein, "C:/BGI/root_integration/data/Ara/Schiefelbein_data/Schiefelbein_doublet_filtered.rds")

