library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)

####################### IFNB PBMC DATA - SEURAT OBJECT #########################

# scenario 1 from 'Benchmark of batch-effect correction methods for single-cell RNA sequencing data' paper
# identical cell types, different technologies
InstallData('ifnb')
LoadData('ifnb')
ifnb.list <- SplitObject(ifnb, split.by = 'stim')


# Standard pre-processing step prior to dimensional reduction techniques
for (i in 1:length(ifnb.list)) {
  ifnb.list[[i]] <- Seurat::NormalizeData(ifnb.list[[i]], normalization.method = 'LogNormalize', scale.factor = 10000)
  # Identify variable features based on a variance stabilizing transformation (vst)
  ifnb.list[[i]] <- Seurat::FindVariableFeatures(ifnb.list[[i]], selection.method = 'vst', nfeatures = 2000)
}

####################### BASE CASE - NO BATCH CORRECTION ########################

ctrl <- ifnb.list$CTRL
stim <- ifnb.list$STIM

rcdata <- merge(x = ctrl, y = stim, add.cell.ids = c('ctrl', 'stim'), merge.data = FALSE)
rcdata <- Seurat::NormalizeData(rcdata)
var.genes <- Seurat::SelectIntegrationFeatures(SplitObject(rcdata, split.by = 'stim'), nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000, selection.method = 'vst')
rcdata <- Seurat::ScaleData(rcdata, features = var.genes)
rcdata <- Seurat::RunPCA(rcdata, features = var.genes, npcs = 40, ndims.print = 1:5, nfeatures.print = 5)
rcdata <- Seurat::FindNeighbors(rcdata, reduction = 'pca', dims = 1:20, k.param = 20)
rcdata <- Seurat::FindClusters(rcdata)      # optimal resolution = 1.1
rcdata <- Seurat::RunUMAP(rcdata, dims = 1:20, reduction = 'pca', n.neighbors = 15, min.dist = 0.5, spread = 1, metric = 'euclidean', seed.use = 1)  
#p1 <- Seurat::DimPlot(rcdata, reduction = 'umap', group.by = 'stim')
#p2 <- Seurat::DimPlot(rcdata, reduction = 'umap', group.by = 'seurat_annotations')
#table(Idents(rcdata), rcdata$stim)


# ADJUSTED RAND INDEX - BATCH MIXING
ari <- dplyr::select(rcdata[[]], stim, seurat_clusters)
ari$stim <- plyr::mapvalues(ari$stim, from = c('CTRL', 'STIM'), to = c(0, 1))
pdfCluster::adj.rand.index(as.numeric(ari$stim), as.numeric(ari$seurat_clusters))         #0.1646627


# ADJUSTED RAND INDEX - CELL PURITY
ari <- dplyr::select(rcdata[[]], seurat_annotations, seurat_clusters)
ari$seurat_annotations <- plyr::mapvalues(ari$seurat_annotations,
                                from = c('B', 'B Activated', 'CD14 Mono', 'CD16 Mono', 'CD4 Memory T', 'CD4 Naive T', 'CD8 T', 'DC', 'Eryth', 'Mk', 'NK', 'pDC', 'T activated'),
                                to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$seurat_annotations), as.numeric(ari$seurat_clusters))       #0.5933426
# optimal resolution = 1.1, ari = 0.1607096, 0.594295


# AVERAGE SILHOUETTE WIDTH (ASW)
library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = rcdata@reductions[['pca']])[, 1:30])
sil.stim <- silhouette(x = as.numeric(x = as.factor(x = rcdata$stim)), dist = dist.matrix)
rcdata$sil.stim <- sil.stim[, 3]
sil.seurat_annotations <- silhouette(x = as.numeric(x = as.factor(x = rcdata$seurat_annotations)), dist = dist.matrix)
rcdata$sil.seurat_annotations <- sil.seurat_annotations[, 3]
mean(rcdata$sil.stim)                       #0.09839579
mean(rcdata$sil.seurat_annotations)         #0.1479373

####################### SEURAT CCA #############################################

ifnb.anchors <- Seurat::FindIntegrationAnchors(object.list = ifnb.list)
ifnb.integrated <- Seurat::IntegrateData(anchorset = ifnb.anchors)


Seurat::DefaultAssay(ifnb.integrated) <- 'integrated'
ifnb.integrated <- Seurat::ScaleData(ifnb.integrated)
ifnb.integrated <- Seurat::RunPCA(ifnb.integrated, npcs = 30)
ifnb.integrated <- Seurat::RunUMAP(ifnb.integrated, reduction = 'pca', dims = 1:30)
ifnb.integrated <- FindNeighbors(ifnb.integrated, dims = 1:10)
ifnb.integrated <- FindClusters(ifnb.integrated, resolution = 0.4)        # optimal resolution = 0.4
#q1 <- Seurat::DimPlot(ifnb.integrated, reduction = 'umap', group.by = 'stim')
#q2 <- Seurat::DimPlot(ifnb.integrated, reduction = 'umap', group.by = 'seurat_annotations')


# ADJUSTED RAND INDEX - BATCH MIXING
ari <- dplyr::select(ifnb.integrated[[]], stim, seurat_clusters)
ari$stim <- plyr::mapvalues(ari$stim, from = c('CTRL', 'STIM'), to = c(0, 1))
pdfCluster::adj.rand.index(as.numeric(ari$stim), as.numeric(ari$seurat_clusters))       #0.00679777


# ADJUSTED RAND INDEX - CELL PURITY
ari <- dplyr::select(ifnb.integrated[[]], seurat_annotations, seurat_clusters)
ari$seurat_annotations <- plyr::mapvalues(ari$seurat_annotations,
                                          from = c('B', 'B Activated', 'CD14 Mono', 'CD16 Mono', 'CD4 Memory T', 'CD4 Naive T', 'CD8 T', 'DC', 'Eryth', 'Mk', 'NK', 'pDC', 'T activated'),
                                          to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$seurat_annotations), as.numeric(ari$seurat_clusters))       #0.6801388
# optimal resolution = 0.4, ari = 0.001772502, 0.8674212


# AVERAGE SILHOUETTE WIDTH (ASW)
#library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = ifnb.integrated@reductions[['pca']])[, 1:30])
sil.stim <- silhouette(x = as.numeric(x = as.factor(x = ifnb.integrated$stim)), dist = dist.matrix)
ifnb.integrated$sil.stim <- sil.stim[, 3]
sil.seurat_annotations <- silhouette(x = as.numeric(x = as.factor(x = ifnb.integrated$seurat_annotations)), dist = dist.matrix)
ifnb.integrated$sil.seurat_annotations <- sil.seurat_annotations[, 3]
mean(ifnb.integrated$sil.stim)                      #0.009264817
mean(ifnb.integrated$sil.seurat_annotations)        #0.2247273

####################### FAST MUTUAL NEAREST NEIGHBORS ##########################

mnn.seurat <- RunFastMNN(object.list = ifnb.list)
mnn.seurat <- Seurat::FindNeighbors(mnn.seurat, reduction = 'mnn', dims = 1:30)
mnn.seurat <- Seurat::FindClusters(mnn.seurat)        # optimal resolution = 0.9
mnn.seurat <- Seurat::RunUMAP(mnn.seurat, reduction = 'mnn', dims = 1:30)
#r1 <- Seurat::DimPlot(mnn.seurat, reduction = 'umap', group.by = 'stim')
#r2 <- Seurat::DimPlot(mnn.seurat, reduction = 'umap', group.by = 'seurat_annotations')


# ADJUSTED RAND INDEX - BATCH MIXING
ari <- dplyr::select(mnn.seurat[[]], stim, seurat_clusters)
ari$stim <- plyr::mapvalues(ari$stim, from = c('CTRL', 'STIM'), to = c(0, 1))
pdfCluster::adj.rand.index(as.numeric(ari$stim), as.numeric(ari$seurat_clusters))       #0.0941706


# ADJUSTED RAND INDEX - CELL PURITY
ari <- dplyr::select(mnn.seurat[[]], seurat_annotations, seurat_clusters)
ari$seurat_annotations <- plyr::mapvalues(ari$seurat_annotations,
                                          from = c('B', 'B Activated', 'CD14 Mono', 'CD16 Mono', 'CD4 Memory T', 'CD4 Naive T', 'CD8 T', 'DC', 'Eryth', 'Mk', 'NK', 'pDC', 'T activated'),
                                          to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$seurat_annotations), as.numeric(ari$seurat_clusters))       #0.6601445
# optimal resolution = 0.9, ari = 0.09400809, 0.6614772


# AVERAGE SILHOUETTE WIDTH (ASW)
#library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = mnn.seurat@reductions[['mnn']])[, 1:30])
sil.stim <- silhouette(x = as.numeric(x = as.factor(x = mnn.seurat$stim)), dist = dist.matrix)
mnn.seurat$sil.stim <- sil.stim[, 3]
sil.seurat_annotations <- silhouette(x = as.numeric(x = as.factor(x = mnn.seurat$seurat_annotations)), dist = dist.matrix)
mnn.seurat$sil.seurat_annotations <- sil.seurat_annotations[, 3]
mean(mnn.seurat$sil.stim)                     #0.02837284
mean(mnn.seurat$sil.seurat_annotations)       #0.1730848

####################### LIGER - SEURAT OBJECT ##################################

liger <- Seurat::NormalizeData(ifnb)
liger <- Seurat::FindVariableFeatures(liger)
var.liger <- Seurat::SelectIntegrationFeatures(Seurat::SplitObject(liger, split.by = 'stim'), nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000, selection.method = 'vst')
liger <- Seurat::ScaleData(liger, split.by = 'stim', features = var.liger)
liger <- RunOptimizeALS(liger, k = 30, split.by = 'stim', do.center = FALSE)
liger <- RunQuantileNorm(liger, split.by = 'stim')
liger <- Seurat::FindNeighbors(liger, reduction = 'iNMF', dims = 1:30)
liger <- Seurat:: FindClusters(liger)         # optimal resolution = 1
liger <- Seurat::RunUMAP(liger, reduction = 'iNMF', dims = 1:ncol(liger[['iNMF']]))
#s1 <- Seurat::DimPlot(liger, reduction = 'umap', group.by = 'stim')
#s2 <- Seurat::DimPlot(liger, reduction = 'umap', group.by = 'seurat_annotations')


# ADJUSTED RAND INDEX - BATCH MIXING
ari <- dplyr::select(liger[[]], stim, seurat_clusters)
ari$stim <- plyr::mapvalues(ari$stim, from = c('CTRL', 'STIM'), to = c(0, 1))
pdfCluster::adj.rand.index(as.numeric(ari$stim), as.numeric(ari$seurat_clusters))       #0.0003922591


# ADJUSTED RAND INDEX - CELL PURITY
ari <- dplyr::select(liger[[]], seurat_annotations, seurat_clusters)
ari$seurat_annotations <- plyr::mapvalues(ari$seurat_annotations,
                                          from = c('B', 'B Activated', 'CD14 Mono', 'CD16 Mono', 'CD4 Memory T', 'CD4 Naive T', 'CD8 T', 'DC', 'Eryth', 'Mk', 'NK', 'pDC', 'T activated'),
                                          to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$seurat_annotations), as.numeric(ari$seurat_clusters))       #0.6897115
# optimal resolution = 1, ari = 0.0003552363, 0.7005685


# AVERAGE SILHOUETTE WIDTH (ASW)
#library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = liger@reductions[['iNMF']])[, 1:30])
sil.stim <- silhouette(x = as.numeric(x = as.factor(x = liger$stim)), dist = dist.matrix)
liger$sil.stim <- sil.stim[, 3]
sil.seurat_annotations <- silhouette(x = as.numeric(x = as.factor(x = liger$seurat_annotations)), dist = dist.matrix)
liger$sil.seurat_annotations <- sil.seurat_annotations[, 3]
mean(liger$sil.stim)                      #0.003251635
mean(liger$sil.seurat_annotations)        #0.02560169


####################### ANALYSIS GRAPHS ########################################

ari.df <- data.frame(x= c(0.594295, 0.8674212, 0.6614772, 0.7005685), y = c(0.8392904, 0.998227498, 0.90599191, 0.9996447637), names = c('Raw', 'Seurat CCA', 'fastMNN', 'LIGER'))
asw.df <- data.frame(x = c(0.1479373, 0.1555964, 0.06459867, 0.02560169), y = c(0.90160421, 0.98472071, 0.95743892, 0.996748365), names = c('Raw', 'Seurat CCA', 'fastMNN', 'LIGER'))

library(ggplot2)
library(ggrepel)
ggplot2::ggplot(ari.df, aes(x,y)) + geom_text_repel(aes(label = names), hjust = 1.5, vjust = 1.6, direction = 'x') +
  geom_point(size = 5, color = c('darkslategray3', 'plum3', 'goldenrod2', 'coral')) +
  xlab('ARI celltype') + ylab('1 - ARI batch') +
  ggplot2::coord_cartesian(xlim = c(0.3, 1), ylim = c(0.5, 1.1)) + theme_bw() +
  theme(axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot2::ggplot(asw.df, aes(x,y)) + geom_text_repel(aes(label = names), hjust = 0, vjust = 1.5, direction = 'x') +
  geom_point(size = 5, color = c('darkslategray3', 'plum3', 'goldenrod2', 'coral')) +
  xlab('ASW celltype') + ylab('1 - ASW batch') +
  ggplot2::coord_cartesian(xlim = c(-0.1, 0.2), ylim = c(0.5, 1.1)) + theme_bw() +
  theme(axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

