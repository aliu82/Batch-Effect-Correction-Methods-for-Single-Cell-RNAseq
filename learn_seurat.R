install.packages("dplyr")
install.packages("CRAN")
install.packages("Seurat")
library(dplyr)
library(Seurat)
library(patchwork)

### source: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#setup-the-seurat-object-1 ###


### 1. initialize Seurat object with non-normalized data - 13209 features ###
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


### 2. explore QC metrics and select cells for further analysis ###
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")                    #percent.mt = percentage of reads that map to the mitochondrial genome
#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


### 3. normalize ###
pbmc <- NormalizeData(pbmc)


### 4. identify highly variable features - 2000 variable features ###
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) 
#top10 <- head(VariableFeatures(pbmc), 10)
#plot1 <- VariableFeaturePlot(pbmc)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2


### 5. scale ###
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


### 6. linear dimensional reduction (PCA) ###
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))    #PC1 represents the first principal direction which cells show the largest variation
#VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#DimPlot(pbmc, reduction = "pca")
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)           #default color scale: magenta-black-yellow
#DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)        #colors that are more clearly separated (heterogenity) are better


### 7. determine dimensionality ###
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)              #significant PCs have low p-values
ElbowPlot(pbmc)                               #only the first 10PCs are really representative of a robust compression of the dataset

#pbmc_start<-pbmc
pbmc<-pbmc_start

### 8. cluster ###
pbmc <- FindNeighbors(pbmc, dims = 1:10)                        #k.param = 350 for knn to be equal to louvain
pbmc <- FindClusters(pbmc, resolution = 0.5)
#pbmc_knn <- FindNeighbors(pbmc, dims = 1:10)                   #head(pbmc_knn, 5)
#pbmc_lovain <- FindClusters(pbmc_knn, resolution = 0.5)        #head(Idents(pbmc_lovain), 5)


### 9. non-linear dimensional reduction (UMAP) ###
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
#DimPlot(RunUMAP(pbmc_knn, dims = 1:10), reduction = "umap")                   #3:14
#DimPlot(RunUMAP(pbmc_lovain, dims = 1:10), reduction = "umap", label=TRUE)
#saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
#markers = FindMarkers(pbmc, ident.1='pbmc3k', genes.use=rownamespbmc_knn@scale.data)


### 10. find biomarkers that represent clusters###
#cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)                                 #head(cluster2.markers, n = 5)
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)              #find all markers distinguishing cluster 5 from 0 and 3
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)                                #reports positive markers

#VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
#VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)                                         #raw counts
#FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
#top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#DoHeatmap(pbmc, features = top10$gene) + NoLegend()                                                             #plotted top 20 markers per cluster


### 11. assign cell type identity to clusters ###
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
#new.cluster.ids <- c("NK", "CD14+ Mono", "B")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
