## library
library(DoubletFinder)

### path.variables ###
Cuperus.path <- "C:/BGI/root_integration/data/Ara/Cuperus/Cuperus_Ara.rds"

Sample.Id <- c('Control','Heatshock')

## read in data 
### load data ###
Cuperus <- readRDS(Cuperus.path) # 2085 cells

##QC
##提取线粒体基因
Cuperus[["percent.mt"]] <- PercentageFeatureSet(Cuperus, pattern='^ATMG')
Cuperus[["percent.ct"]] <- PercentageFeatureSet(Cuperus, pattern='^ATCG')

violin <- VlnPlot(Cuperus,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
                  pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
                  ncol = 4)

ggsave("C:/BGI/root_integration/figures/Ara/QC/Cuperus-vlnplot-before-qc.pdf", plot = violin, width = 15, height = 6) 

plot1 <- FeatureScatter(Cuperus, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cuperus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("C:/BGI/root_integration/figures/Ara/QC/Rybel_pearplot-before-qc.pdf", plot = pearplot, width = 12, height = 5) 

##set QC critiria    2028 cells remained
Cuperus <- subset(Cuperus, subset=nFeature_RNA > 500 & nFeature_RNA < 9500
                & percent.mt < 1 &  percent.ct < 0.3)

## post-QC
violin <- VlnPlot(Cuperus,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ct"),
                  pt.size = 0.1, #不需要显示点，可以设置pt.size = 0
                  ncol = 4)
ggsave("C:/BGI/root_integration/figures/Ara/QC/Cuperus-vlnplot-after-qc.pdf", 
       plot = violin, width = 15, height = 6) 

plot1 <- FeatureScatter(Cuperus, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cuperus, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pearplot <- CombinePlots(plots = list(plot1, plot2), nrow=1, legend="none") 
ggsave("C:/BGI/root_integration/figures/Ara/QC/Cuperus_pearplot-after-qc.pdf", plot = pearplot, width = 12, height = 5) 

#### Normalization ####
Cuperus <- NormalizeData(Cuperus, normalization.method = "LogNormalize", scale.factor = 10000)


### feature selection: HVGs ### -----method 2
Cuperus <- FindVariableFeatures(Cuperus, selection.method = "vst", nfeatures = 2000)


### scale data ###
Cuperus <- ScaleData(Cuperus) # uncorrected

### dimensionality reduction: PCA ###
Cuperus <- RunPCA(Cuperus)


## elbow plot - Rybel ##
pca.elbow.plot <- ElbowPlot(Cuperus, ndims = 50, reduction = "pca")
png(paste0(CK.path,"pca.elbow.plot_ckn.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

## run umap
Cuperus <- RunUMAP(Cuperus, reduction = "pca", dims = 1:50)

## FindNeighbors
Cuperus <- FindNeighbors(Cuperus, reduction = "pca", dims = 1:50)

## FindClusters
Cuperus <- FindClusters(Cuperus, resolution = 0.5)

##split samples and run doubletfiner for each samples 
Cuperus.list <- SplitObject(Cuperus, split.by = "Sample")
control <- Cuperus.list[[1]]
heatshock <- Cuperus.list[[2]]


## control pKֵ CKN
sweep.res.list <- paramSweep_v3(control, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(control)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(control$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(control)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
control <- doubletFinder_v3(control, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(control, reduction = "umap", group.by = "DF.classifications_0.25_0.01_8")
ggsave('control_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----1048 cells remained
control <- subset(control, subset = DF.classifications_0.25_0.01_8 == 'Singlet')
dim(control)

##  heatshock  pKֵ CKN
sweep.res.list <- paramSweep_v3(heatshock, PCs = 1:50)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## ?ų????ܼ?????ͬԴdoublets???Ż???????doublets??��
DoubletRate = ncol(heatshock)*8*1e-6                    # 5000ϸ????Ӧ??doublets rate??3.9%
homotypic.prop <- modelHomotypic(heatshock$seurat_clusters)   # ?????ṩcelltype
nExp_poi <- round(DoubletRate*ncol(heatshock)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ʹ??ȷ???õĲ???????doublets
heatshock <- doubletFinder_v3(heatshock, PCs = 1:50, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F)

## ????չʾ????????????pbmc@meta.data??
p <- DimPlot(heatshock, reduction = "umap", group.by = "DF.classifications_0.25_0.05_6")
ggsave('heatshock_doublet_umap.pdf', p, path = "C:/BGI/root_integration/figures/Ara/QC/",
       width = 7, height = 5)

## remove doublets ----966 cells remained
heatshock <- subset(heatshock, subset = DF.classifications_0.25_0.05_6 == 'Singlet')
dim(heatshock)


### merge all the datasets after doubletfinder------ 2014 final cells
Cuperus <- merge(control, heatshock)
dim(Cuperus)

### save data ###  
saveRDS(Cuperus, "C:/BGI/root_integration/data/Ara/Cuperus/Cuperus_doublet_filtered.rds")

