library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(data.table)
library(readxl)
library(ggsci)

### path variables ###
Schiefelbein.path <- "C:/BGI/root_integration/data/Ara/Schiefelbein_data/Schiefelbein_doublet_filtered.rds"
Wang.path <- "C:/BGI/root_integration/data/Ara/Wang_data/Wang_doublet_filtered.rds"
timmermans.path <- "C:/BGI/root_integration/data/Ara/Timmermans_data/Timmermans_doublet_filtered.rds"
Cuperus.path <- "C:/BGI/root_integration/data/Ara/Cuperus/Cuperus_doublet_filtered.rds"
Rybel.path <- "C:/BGI/root_integration/data/Ara/Rybel/Rybel_doublet_filtered.rds"

Ara.path <- "C:/BGI/root_integration/figures/Ara/Integration/"
harmony.samples.path <- "C:/BGI/root_integration/figures/Ara/Integration/harmony/"


### load data ###
Schiefelbein <- readRDS(Schiefelbein.path) # 10778 cells
Wang <- readRDS(Wang.path) # 11337 cells
Timmermans <- readRDS(timmermans.path) # 5551 cells
Cuperus <- readRDS(Cuperus.path)
Rybel <- readRDS(Rybel.path)

### add metadata to Wang data
Wang$Sample <- 'Wang'

### Merge ara together ###   34747 cells in total
Ara_root <- merge(Schiefelbein, c(Wang,Timmermans,Cuperus,Rybel))

#### Normalization ####
Ara_root <- NormalizeData(Ara_root, normalization.method = "LogNormalize", scale.factor = 10000)

### feature selection ###
# select the top 4000 HVGs of each study and only consider genes that are HVGs in all studies
Schiefelbein <- FindVariableFeatures(Schiefelbein, selection.method = "vst", nfeatures = 5000)
Wang <- FindVariableFeatures(Wang, selection.method = "vst", nfeatures = 5000)
Timmermans <- FindVariableFeatures(Timmermans, selection.method = "vst", nfeatures = 5000)
Cuperus <- FindVariableFeatures(Cuperus, selection.method = "vst", nfeatures = 5000)
Rybel <- FindVariableFeatures(Rybel, selection.method = "vst", nfeatures = 5000)

# option 2: only consider genes that are HVGs in both studies
ara.HVGs <- Reduce(intersect,list(VariableFeatures(Schiefelbein), 
                                  VariableFeatures(Wang),
                                  VariableFeatures(Timmermans),
                                  VariableFeatures(Cuperus),
                                  VariableFeatures(Rybel))) # 691 HVGs
# set HVGs
VariableFeatures(Ara_root) <- ara.HVGs


### feature selection: HVGs ### -----method 2
#Ara_root <- FindVariableFeatures(Ara_root, selection.method = "vst", nfeatures = 2000)

### feature selection: HVGs ### ----- method 3
dataset.list <- SplitObject(Ara_root, split.by = "orig.ident")
featuresdataset <- SelectIntegrationFeatures(object.list = dataset.list)
VariableFeatures(Ara_root) <- featuresdataset

### scale data ###
Ara_root <- ScaleData(Ara_root) # uncorrected

### dimensionality reduction: PCA ###
Ara_root <- RunPCA(Ara_root, features = VariableFeatures(object = Ara_root))


## elbow plot  ##
pca.elbow.plot <- ElbowPlot(Ara_root, ndims = 50, reduction = "pca")
png(paste0(Ara.path,"pca.elbow.plot.png"), width=1000,height=600,units="px")
print(pca.elbow.plot)
dev.off()

# visualize PCs
PC1_2.sample.plot <- DimPlot(object = Ara_root, reduction = "pca", pt.size = .1, group.by = "Sample")
png(paste0(Ara.path,"PC1_2.sample.png"), width=1000,height=1000,units="px")
print(PC1_2.sample.plot)
dev.off()

PC1_2.studies.plot <- DimPlot(object = Ara_root, reduction = "pca", pt.size = .1, group.by = "Study")
png(paste0(Ara.path,"PC1_2.studies.png"), width=1000,height=1000,units="px")
print(PC1_2.studies.plot)
dev.off()


# save genes making up the PCs  
sink(paste0(Ara.path, "PC_genes.txt"))
print(Ara_root[["pca"]], dims = 1:50, nfeatures = 20)
sink()


### integration with harmony ###
# harmonize samples
Ara.harmony <- Ara_root %>% RunHarmony("Sample", theta = 2, reduction.save = "harmony_theta2", plot_convergence = TRUE) 

# # ## harmony elbow plot ##
harmony.elbow.plot <- ElbowPlot(Ara.harmony, ndims = 50, reduction = "harmony_theta2")
png(paste0(harmony.samples.path,"harmony_theta2.elbow.plot.png"), width=1000,height=600,units="px")
print(harmony.elbow.plot)
dev.off()
# #
## explore harmony coordinates ##
## visualize PCs ##

harmony.PC1_2.sample.plot <- DimPlot(object = Ara.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "Sample")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.sample.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.sample.plot)
dev.off()

harmony.PC1_2.studies.plot <- DimPlot(object = Ara.harmony, reduction = "harmony_theta2", pt.size = .1, group.by = "Study")
png(paste0(harmony.samples.path,"harmony_theta2.PC1_2.studies.png"), width=1000,height=1000,units="px")
print(harmony.PC1_2.studies.plot)
dev.off()

# # ### save genes making up the PCs ###
sink(paste0(harmony.samples.path, "harmony_PC_genes.txt"))
print(Ara.harmony[["harmony_theta2"]], dims = 1:50, nfeatures = 20)
sink()
# 
## save data ###
saveRDS(Ara.harmony, paste0(harmony.samples.path, "all_Ara_harmony.rds"))

### explore different numbers of harmony PCs ###

# visualize marker gene expression
### explore different numbers of harmony PCs ###
# root cap cell (RCC) genes
RCC.genes <- c("AT1G33280","AT1G79580")

# root meristematic cell (RMC) genes 
RMC.genes <- c("AT1G73590",'AT3G20840','AT2G04025','AT3G25980','AT5G13840')

# QC
QC.genes <- c('AT1G02870')

# epidermal cell/root hair (RH) genes
RH.genes <- c('AT5G49270','AT4G22080','AT2G03720')

# epidermal cell/non root hair (NRH) genes
NRH.gene <- c('AT3G02240','AT2G37260','AT1G65310','AT1G79840')

# cortex (Co) genes
RC.genes <- c("AT1G62510",'AT5G02000','AT5G64620')

# endodermis (CS+) genes
RE.genes <- c('AT5G57620','AT4G20140','AT5G06200','AT3G11550','RALF1')

# phloem (P)
phloem.genes <- c('AT4G19840','AT1G05760','AT2G15310','AT1G79430')

# stele/vacular cell genes 
VC.genes <- c('AT3G25710','AT5G48070','AT2G31083',
              'AT1G22710','AT3G23430','AT1G32450')

# xylem (X)
xylem.genes <- c('AT1G68810','AT1G20850',
                 'AT5G44030','AT3G16920','AT2G37090')

#trichoblast
trichoblast.gemes <- c('AT5G49270','AT1G63450','AT2G39690','AT1G07795')

##procabium
procabium.genes <- c('AT5G57130', 'AT4G32880')

###pericyle
pericyle.genes <- c('AT1G32450','AT5G01740')

markers <- c(RCC.genes, RMC.genes,  RH.genes, 
             NRH.gene, RC.genes, RE.genes, phloem.genes,
             VC.genes, xylem.genes, QC.genes)

marker_symbol <- c(
                   'GLV4',
                   'GL2',
                   'ATC/VIF2',
                   'MYB36',
                   'SGN3',
                   'CASP4',
                   'CASP2',
                   'AtPP2-A1',
                   'RTM1',
                   'ARFB1A',
                   'APL',
                   'TMO5',
                   'AtXTH20',
                   'CLE5',
                   'SUC2',
                   'PHO1',
                   'NPF7.3',
                   'T5L1',
                   'XCP2',
                   'CESA4',
                   'AtCTL2',
                   'IRX9')

library(MySeuratWrappers)  


#??ɫ????UMAPͼһ??
col <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#??ɫ????  
p <- VlnPlot(Ara.harmony, features = markers,  
             stack=T, pt.size=0,  
             #cols = col,#??ɫ  
             direction = "horizontal", #ˮƽ??ͼ
             x.lab = '', y.lab = '')+#?????᲻?????κζ???  
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#????ʾ?????̶?
p <- p + scale_fill_d3('category20')


ggsave('Ara_stacked_markers.pdf', p, path =  "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 10, height = 5)
## colors 
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)
p <- StackedVlnPlot(Ara.harmony, markers, pt.size=0, cols=my36colors)


###dot plot
p <- DotPlot(Ara.harmony, features = markers) +
    coord_flip() +
    theme_bw()+ #  
    theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+ 
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3)) 

ggsave('Ara_dotplot_markers.pdf', p, path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 6, height = 8)
### For ara ###
dims <- c(40,50)
for(d in dims){
  
  # create folder
  dir.create(paste0(harmony.samples.path, "dim", d, "_annotation"))
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # run UMAP
  Ara.harmony <- RunUMAP(Ara.harmony, reduction = "harmony_theta2", dims = 1:d, seed.use=1)
  
  ## batch plots ##
  # no.1 overview - by patients | by study
  sample.batch.plot <- DimPlot(Ara.harmony, reduction = "umap", group.by = "Sample", pt.size = 0.01)
  study.batch.plot <- DimPlot(Ara.harmony, reduction = "umap", group.by = "Study", pt.size = 0.1)
  batch1.plot <- sample.batch.plot + study.batch.plot
  png(paste0(dim.path, "UMAP_dim", d, ".batch.png"), width=1500, height=600, units="px")
  print(batch1.plot)
  dev.off()
  
  # no.2 - colored by sample
  batch2.plot <- DimPlot(Ara.harmony, reduction = "umap", group.by = "Sample", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".sample.batch1.png"), width=1200, height=800, units="px")
  print(batch2.plot)
  dev.off()
  
  # no.3 - split by study
  batch3.plot <- DimPlot(Ara.harmony, reduction = "umap", group.by = "Study", split.by = "Sample", pt.size = 0.01, ncol = 4)
  png(paste0(dim.path, "UMAP_dim", d, ".sample.batch2.png"), width=1500, height=1200, units="px")
  print(batch3.plot)
  dev.off()
  
  # no.4 - colored by study
  batch4.plot <- DimPlot(Ara.harmony, reduction = "umap", group.by = "Study", pt.size = 0.01)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch1.png"), width=1400, height=1200, units="px")
  print(batch4.plot)
  dev.off()
  
  # no.5 - split by study
  batch5.plot <- DimPlot(Ara.harmony, reduction = "umap", group.by = "Study", split.by = "Study", pt.size = 0.01, ncol = 3)
  png(paste0(dim.path, "UMAP_dim", d, ".study.batch2.png"), width=1800, height=600, units="px")
  print(batch5.plot)
  dev.off()
  
  # ## compare original annotation ##
  # annotation <- DimPlot(CK.harmony, reduction = "umap", split.by = "Condition", group.by = "celltype",
  #                       label = TRUE, repel = TRUE, pt.size = 0.01) + NoLegend()
  # png(paste0(dim.path, "UMAP_dim", d, "healthy_fibrotic_2.png"), width=1800, height=900, units="px")
  # print(annotation)
  # dev.off()
  
  # ## compare original annotation in healthy vs. fibrotic data ##
  # plot <- compare / annotation
  # png(paste0(dim.path,"UMAP_dim", d, "healthy_fibrotic_3.png"), width=2000, height=2000, units="px")
  # print(plot)
  # dev.off()
  
  ## QC metrics at cluster level ##
  cluster.qc.heatmap <- FeaturePlot(Ara.harmony, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.2, ncol = 3) & 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  png(paste0(dim.path,"UMAP_dim", d, "cluster.qc.heatmap.png"), width=1500,height=500,units="px")
  print(cluster.qc.heatmap)
  dev.off()
  
  ## gene expression heatmaps ##
  # # overview feature plot
  # lineage.overview <- FeaturePlot(CK.harmony, features = lineage.genes, pt.size = 0.5, ncol = 4) & 
  #   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
  # png(paste0(dim.path,"UMAP_dim", d, "overview.markers.png"), width=1600,height=500,units="px")
  # print(lineage.overview )
  # dev.off()
  
}

### clustering ###
# only for selected number of dims (for reasons of computational efficiency!)
res <- seq(0.2,0.8,0.1)
for(d in dims){
  
  # set path variable
  dim.path <- paste0(harmony.samples.path, "dim", d, "_annotation/")
  
  # visualization with UMAP
  Ara.harmony <- RunUMAP(Ara.harmony, reduction = "harmony_theta2", dims=1:d, seed.use=1)
  
  # plot batch effect
  batch.plot <- DimPlot(Ara.harmony, reduction = "umap", group.by = "Study", pt.size = 0.1)
  
  # create kNN graph
  Ara.harmony <- FindNeighbors(Ara.harmony, reduction = "harmony_theta2", dims = 1:d)
  
  for (r in res) {
    
    Ara.harmony <- FindClusters(Ara.harmony, reduction = "harmony_theta2", resolution = r)
    umap.plot <- DimPlot(Ara.harmony, reduction = "umap", repel = T, label = TRUE, pt.size = 0.1)
    
    # create eval plot
    eval.plot <- umap.plot + batch.plot
    png(paste0(dim.path, "UMAP_dim", d, "_res", r, ".png"), width=1500, height=600, units="px")
    print(eval.plot)
    dev.off()
    
    # clustering at sample level
    sample.clustering <- DimPlot(Ara.harmony, group.by = "seurat_clusters", split.by = "Sample", pt.size = 0.1, ncol = 4)
    png(paste0(dim.path,"UMAP_dim", d, "_res", r,".sample.clustering.png"), width=1600,height=1600,units="px")
    print(sample.clustering)
    dev.off()
    
    # clustering at study level
    study.clustering <- DimPlot(Ara.harmony, group.by = "seurat_clusters", split.by = "Study", pt.size = 0.01, ncol = 3)
    png(paste0(dim.path,"UMAP_dim", d, "_res", r, ".study.clustering.png"), width=1800,height=1200,units="px")
    print(study.clustering)
    dev.off()
    
  }
}

### preliminary clustering ###
Ara.harmony <- FindNeighbors(Ara.harmony, reduction = "harmony_theta2", dims = 1:50)
Ara.harmony <- FindClusters(Ara.harmony, reduction = "harmony_theta2", resolution = 0.5)
# run UMAP
Ara.harmony <- RunUMAP(Ara.harmony, reduction = "harmony_theta2", dims=1:50, seed.use=1)

# run tsne
Ara.harmony <- RunTSNE(Ara.harmony, reduction = "harmony_theta2", dims = 1:50, seed.use = 1)
### save data ###
saveRDS(Ara.harmony, paste0(CK.path, "CK_harmony.rds"))


## compute cluster marker genes ###
Ara.markers <- FindAllMarkers(Ara.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
Ara.top10.markers <- Ara.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Ara.markers, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/CK_Marker_genes_res0.4.csv")
write.csv(Ara.top10.markers, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/top10_marker_genes_res0.4.csv")

### RUN FINALMARKERS FOR CLUSTER 6
cluster6.markers <- FindMarkers(Ara.harmony, only.pos = TRUE, ident.1 = 6, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(cluster6.markers, file = "C:/BGI/root_integration/figures/Ara/Integration/harmony/cluster6_marker_genes_res0.5.csv")

### RUN FINALMARKERS FOR CLUSTER 7
cluster7.markers <- FindMarkers(Ara.harmony, only.pos = TRUE, ident.1 = 7, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(cluster7.markers, file = "C:/BGI/root_integration/figures/Ara/Integration/harmony/cluster7_marker_genes_res0.5.csv")

### RUN FINALMARKERS FOR CLUSTER 14
cluster14.markers <- FindMarkers(Ara.harmony, only.pos = TRUE, ident.1 = 14, min.pct = 0.25, logfc.threshold = 0.25) 
write.csv(cluster14.markers, file = "C:/BGI/root_integration/figures/Ara/Integration/harmony/cluster14_marker_genes_res0.5.csv")

### cell type annotation ###
cluster.annotation <- c("Root cap cell 1", 'Stele 1', 'Non-hair root epidermal cell 1',
                        'Root hair cell 1', 'Root cortex','Stele 2', 
                        'Unknown','Quiescent Center', 'Phloem 1', 
                        'Root hair cell 2', 'Non-hair root epidermal cell 2', 'Root cap cell 2',
                        'Meristem 1', 'Xylem 1', 'Root endodermis',
                        'Meristem 2', 'Phloem 2','Root cap cell 3', 
                        'Xylem 2')


names(cluster.annotation) <- levels(Ara.harmony)
Ara.harmony <- RenameIdents(Ara.harmony, cluster.annotation)

# add cell types to meta data
cell.data <- data.table(barcode = colnames(Ara.harmony),
                        celltype = Idents(Ara.harmony))

cell.data <- data.frame(cell.data, row.names = cell.data$barcode)
cell.data$barcode <- NULL
Ara.harmony <- AddMetaData(Ara.harmony, cell.data, col.name = "celltype")

p <- DimPlot(Ara.harmony, group.by = 'celltype', label = T, label.size = 3,pt.size = 0.03, repel = T)
p <- p + scale_color_d3('category20')

ggsave('Ara_annotated_umap_res0.5.pdf', p, path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 8, height = 5)


### cell lineage annotation ###
cell.data <- data.table(barcode = colnames(Ara.harmony),
                        celltype = Idents(Ara.harmony))

lineage.annotation <- c("Root cap cell (RCC)", "Stele", "Non root hair (NRH)",
                        "Root hair (RH)", "Cortex", "Stele",
                        "Unknown", "QC","Stele(P)",
                        "Root hair (RH)","Non root hair (NRH)","Root cap cell (RCC)",
                        "Root Meristem", "Stele(X)", 'Endodermis',
                        'Root Meristem', "Stele(P)", 'Root cap cell (RCC)',
                        "Stele(X)")

# lineage.annotation <- c("Unknown","Root cap cell (RCC)","Protophloem",
#                         "Pericycle","Stele","Cortex","Epidermal cell/non root hair (NRH)",
#                         "Columella","Endodermis","QC",
#                         "Trichoblast","Root cap cell (RCC)","Stele",
#                         "Root cap cell (RCC)","Stele", "Unknown")

lineage.data <- data.table(celltype = cluster.annotation, lineage = lineage.annotation)
meta.data <- merge(cell.data, lineage.data, by = "celltype")
meta.data <- data.frame(meta.data, row.names = meta.data$barcode)
meta.data$barcode <- NULL
meta.data$celltype <- NULL
Ara.harmony <- AddMetaData(Ara.harmony, meta.data, col.name = "lineage")

p <- DimPlot(Ara.harmony, group.by = 'lineage', label = T, label.size = 3, pt.size = 0.01, repel = T)
p <- p + scale_color_d3('category20')

ggsave('Ara_annotated_umap_lineage.pdf', p, path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 8, height = 5)

### save data ### 34747 cells remained
saveRDS(Ara.harmony, "C:/BGI/root_integration/data/Ara/combined_Ara_annotated_doubletremoved.rds")



##colors
# SAVE FIGURES  ----group by study and samples
library(ggsci)
p <- DimPlot(Ara.harmony, reduction = "umap", label = T,   
             #cols= allcolour, #??????ɫ  
             pt.size = 0.01,#???õ??Ĵ?С  
             group.by = 'Study',
             repel = T)+
  theme(axis.text.y = element_blank(),   #ȥ??y???̶?ע??
        axis.ticks.y = element_blank(),    #ȥ??y???̶?
        axis.text.x = element_blank(),   #ȥ??x???̶?ע??
        axis.ticks.x = element_blank())+  #ȥ??x???̶?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


#p <- p+ scale_color_npg()

ggsave('Ara_group_by_study_umap.pdf', path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 8, height = 6)

## group by celltype but split by study
p <- DimPlot(Ara.harmony, reduction = "umap", label = F,   
             #cols= allcolour, #??????ɫ  
             pt.size = 0.01,#???õ??Ĵ?С  
             group.by = 'celltype',
             split.by = 'Study',
             ncol = 3)+
  theme(axis.text.y = element_blank(),   #ȥ??y???̶?ע??
        axis.ticks.y = element_blank(),    #ȥ??y???̶?
        axis.text.x = element_blank(),   #ȥ??x???̶?ע??
        axis.ticks.x = element_blank())+  #ȥ??x???̶?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


p <- p+ scale_color_d3('category20')

ggsave('Ara_split_by_study_umap.pdf', path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 15, height = 8)

## group by samples
p <- DimPlot(Ara.harmony, reduction = "umap", label = F,   
             #cols= allcolour, #??????ɫ  
             pt.size = 0.01,#???õ??Ĵ?С  
             group.by = 'Sample',
             repel = T)+
  theme(axis.text.y = element_blank(),   #ȥ??y???̶?ע??
        axis.ticks.y = element_blank(),    #ȥ??y???̶?
        axis.text.x = element_blank(),   #ȥ??x???̶?ע??
        axis.ticks.x = element_blank())+  #ȥ??x???̶?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


#p <- p+ scale_color_lancet()

ggsave('Ara_group_by_sample_umap.pdf', path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 8, height = 6)

## group by celltype but split by sample
p <- DimPlot(Ara.harmony, reduction = "umap", label = F,   
             #cols= allcolour, #??????ɫ  
             pt.size = 0.01,#???õ??Ĵ?С  
             group.by = 'celltype',
             split.by = 'Sample',
             ncol = 4)+
  theme(axis.text.y = element_blank(),   #ȥ??y???̶?ע??
        axis.ticks.y = element_blank(),    #ȥ??y???̶?
        axis.text.x = element_blank(),   #ȥ??x???̶?ע??
        axis.ticks.x = element_blank())+  #ȥ??x???̶?
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


p <- p+ scale_color_d3('category20')

ggsave('Ara_split_by_sample_umap.pdf', path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 15, height = 12)


### extracting umap infor for figures---- celltypes
umap = Ara.harmony@reductions$umap@cell.embeddings %>%  #??????Ϣ
  as.data.frame() %>% 
  cbind(cell_type = Ara.harmony@meta.data$celltype) # ע?ͺ???label??Ϣ ????Ϊcell_type

### extracting umap infor for figures -----cell lineage
umap = Ara.harmony@reductions$umap@cell.embeddings %>%  #??????Ϣ
  as.data.frame() %>% 
  cbind(cell_lineage = Ara.harmony@meta.data$lineage) # ע?ͺ???label??Ϣ ????Ϊcell_type

head(umap)
## colors
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
## generate figures
p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
  geom_point(size = 0.01 , alpha =1 ) + 
  scale_color_manual(values = allcolour)

p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_lineage)) +
  geom_point(size = 0.01 , alpha =1 ) + 
  scale_color_manual(values = allcolour)

p2 <- p  +
  theme(panel.grid.major = element_blank(), #????????
        panel.grid.minor = element_blank(), #????????
        panel.border = element_blank(), #?߿?
        axis.title = element_blank(),  #??????
        axis.text = element_blank(), # ?ı?
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #????ɫ
        plot.background=element_rect(fill="white"))

p3 <- p2 +         
  theme(
    legend.title = element_blank(), #ȥ??legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=12), #????legend??ǩ?Ĵ?С
    legend.key.size=unit(1,'cm') ) +  # ????legend??ǩ֮???Ĵ?С
  guides(color = guide_legend(override.aes = list(size=2))) #????legend?? ???Ĵ?С 


p4 <- p3 + 
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p4 <- p4+scale_color_d3("category20")

p4 <- p4+		scale_color_d3('category20c')

##celltype
cell_type_med <- umap %>%
  group_by(cell_type) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

## cell lineage
cell_lineage_med <- umap %>%
  group_by(cell_lineage) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

library(ggrepel)
##lineage
p4 <- p4 + geom_text_repel(aes(label= cell_lineage), fontface="bold", 
                           data = cell_lineage_med,
                           point.padding=unit(1.6, "lines"))

##cell type
p4 <- p4 + geom_text_repel(aes(label= cell_type), fontface="bold", 
                           data = cell_type_med,
                           point.padding=unit(1.6, "lines"))

##
ggsave('Ara_annotated_umap_refined.pdf', p4, path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 12, height = 8)
##
ggsave('Ara_annotated_lineage_refined.pdf', p4, path = "C:/BGI/root_integration/figures/Ara/Integration/harmony/dim50_annotation/",
       width = 12, height = 8)

##
ggsave('CK_annotated_lineage_doublet_removed_updated.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)
##
ggsave('CK_annotated_lineage_doublet_removed_labeled.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)

##
ggsave('CK_annotated_umap_doublet_removed_labeled.pdf', p4, path = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/",
       width = 12, height = 8)
### save data ###
saveRDS(CK.harmony, "C:/BGI/mianhua/CK_harmony/combined_CK_annotated.rds")


