library(Seurat)
library(SeuratData)
library(patchwork)

####################### PANCREAS DATA - SEURAT OBJECT ##########################

# scenario 3 from 'Benchmark of batch-effect correction methods for single-cell RNA sequencing data' paper
# multiple batches with significantly different distributions of cell types

# Load data and create Seurat objects
# nCount_RNA = total number of molecules detected within a cell
# nFeature_RNA = number of genes detected in each cell
# RNA_snn_res.0.8
InstallData('panc8')
LoadData('panc8')
panc.list <- Seurat::SplitObject(panc8, split.by = 'tech')


# Standard pre-processing step prior to dimensional reduction techniques
for (i in 1:length(panc.list)) {
  panc.list[[i]] <- Seurat::NormalizeData(panc.list[[i]])
  # Identify variable features based on a variance stabilizing transformation (vst)
  panc.list[[i]] <- Seurat::FindVariableFeatures(panc.list[[i]], selection.method = 'vst', nfeatures = 2000)
}

####################### BASE CASE - NO BATCH CORRECTION ########################

celseq <- panc.list$celseq
celseq2 <- panc.list$celseq2
fluidigmc1 <- panc.list$fluidigmc1
smartseq2 <- panc.list$smartseq2
indrop <- panc.list$indrop


gcdata <- merge(x = celseq, y = list(celseq2, fluidigmc1, smartseq2, indrop), add.cell.ids = c('celseq', 'celseq2', 'fluidigmc1', 'smartseq2', 'indrop'), merge.data = FALSE)
gcdata <- Seurat::NormalizeData(gcdata)
var.genes <- Seurat::SelectIntegrationFeatures(SplitObject(gcdata, split.by = 'tech'), nfeatures = 2000, fvf.nfeatures = 2000, selection.method = 'vst')
gcdata <- Seurat::ScaleData(gcdata, features = var.genes)
gcdata <- Seurat::RunPCA(gcdata, features = var.genes, npcs = 40, ndims.print = 1:5, nfeatures.print = 5)
#Seurat::DimPlot(gcdata, reduction = 'pca', group.by = 'tech')
gcdata <- Seurat::FindNeighbors(gcdata, reduction = 'pca', dims = 1:20, k.param = 20)
gcdata <- Seurat::FindClusters(gcdata)    #used default resolution for consistency, optimal resolution = 0.1
gcdata <- Seurat::RunUMAP(gcdata, dims = 1:20, reduction = 'pca', n.neighbors = 15, min.dist = 0.5, spread = 1, metric = 'euclidean')  
#p1 <- Seurat::DimPlot(gcdata, reduction = 'umap', group.by = 'tech')
#p2 <- Seurat::DimPlot(gcdata, reduction = 'umap', group.by = 'celltype')

# how the distribution of celltype varies per tech
Seurat::VlnPlot(subset, 'nFeature_RNA', group.by = 'celltype') + ggplot2::theme(legend.position = 'none') + ggplot2::ggtitle('smartseq2')
#subset <- subset(gcdata, subset = (tech == 'smartseq2'))


# ADJUSTED RAND INDEX - BATCH MIXING (lower value is better)
ari <- dplyr::select(gcdata[[]], tech, seurat_clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c('celseq', 'celseq2', 'fluidigmc1', 'smartseq2', 'indrop'), to = c(0, 1, 2, 3, 4))
pdfCluster::adj.rand.index(as.numeric(ari$tech), as.numeric(ari$seurat_clusters))         #0.1485275


# ADJUSTED RAND INDEX - CELL PURITY (higher value is better)
ari <- dplyr::select(gcdata[[]], celltype, seurat_clusters)
ari$celltype <- plyr::mapvalues(ari$celltype,
                                from = c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma',
                                'macrophage', 'mast', 'quiescent_stellate', 'schwann'),
                                to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$celltype), as.numeric(ari$seurat_clusters))       #0.3665525
# optimal resolution = 0.1, ari = 0.2050022, 0.5642438


# AVERAGE SILHOUETTE WIDTH (ASW)
library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = gcdata@reductions[['pca']])[, 1:30])
sil.tech <- silhouette(x = as.numeric(x = as.factor(x = gcdata$tech)), dist = dist.matrix)
gcdata$sil.tech <- sil.tech[, 3]
sil.celltype <- silhouette(x = as.numeric(x = as.factor(x = gcdata$celltype)), dist = dist.matrix)
gcdata$sil.celltype <- sil.celltype[, 3]
mean(gcdata$sil.tech)             #0.0244474
mean(gcdata$sil.celltype)         #0.2145548

####################### SEURAT CCA #############################################

panc.anchors <- Seurat::FindIntegrationAnchors(object.list = panc.list)
panc.integrated <- Seurat::IntegrateData(anchorset = panc.anchors)

# Switch to integrated assay
Seurat::DefaultAssay(panc.integrated) <- 'integrated'
# Run standard workflow for visualization and clustering
panc.integrated <- Seurat::ScaleData(panc.integrated)
panc.integrated <- Seurat::RunPCA(panc.integrated, npcs = 30)
panc.integrated <- Seurat::RunUMAP(panc.integrated, reduction = 'pca', dims = 1:30)
panc.integrated <- FindNeighbors(panc.integrated, dims = 1:10)
panc.integrated <- FindClusters(panc.integrated)      # optimal resolution = 0.15
#q1 <- Seurat::DimPlot(panc.integrated, reduction = 'umap', group.by = 'tech')
#q2 <- Seurat::DimPlot(panc.integrated, reduction = 'umap', group.by = 'celltype')
#table(Idents(panc.integrated), panc.integrated$tech)


# ADJUSTED RAND INDEX - BATCH MIXING
ari <- dplyr::select(panc.integrated[[]], tech, seurat_clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c('celseq', 'celseq2', 'fluidigmc1', 'smartseq2', 'indrop'), to = c(0, 1, 2, 3, 4))
pdfCluster::adj.rand.index(as.numeric(ari$tech), as.numeric(ari$seurat_clusters))       #0.01952463


# ADJUSTED RAND INDEX - CELL PURITY
ari <- dplyr::select(panc.integrated[[]], celltype, seurat_clusters)
ari$celltype <- plyr::mapvalues(ari$celltype,
                                from = c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma',
                                'macrophage', 'mast', 'quiescent_stellate', 'schwann'),
                                to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$celltype), as.numeric(ari$seurat_clusters))       #0.4770111
#resolution = 0.15, ari = 0.006943822, 0.9371157


# AVERAGE SILHOUETTE WIDTH (ASW)
#library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = panc.integrated@reductions[['pca']])[, 1:30])
sil.tech <- silhouette(x = as.numeric(x = as.factor(x = panc.integrated$tech)), dist = dist.matrix)
panc.integrated$sil.tech <- sil.tech[, 3]
sil.celltype <- silhouette(x = as.numeric(x = as.factor(x = panc.integrated$celltype)), dist = dist.matrix)
panc.integrated$sil.celltype <- sil.celltype[, 3]
mean(panc.integrated$sil.tech)            #-0.1074509
mean(panc.integrated$sil.celltype)        #0.2974331


# other indicators of measuring batch-effect, but couldn't get to work
#lisi <- lisi::compute_lisi(panc.integrated@assays$RNA@counts, ari, c('celltype', 'seurat_clusters'))
#kbet <- kBET::kBET(panc.integrated@assays, predictions)

####################### FAST MUTUAL NEAREST NEIGHBORS ##########################

# RunFastMNN() can be found in the 'seurat-wrappers' github for Seurat
# https://github.com/satijalab/seurat-wrappers
# there is also a batchelor::fastMNN() function that inputs a SingleCellExperiment object instead of a Seurat object
mnn.seurat <- RunFastMNN(object.list = panc.list)
mnn.seurat <- Seurat::FindNeighbors(mnn.seurat, reduction = 'mnn', dims = 1:30)
mnn.seurat <- Seurat::FindClusters(mnn.seurat)      # optimal resolution = 0.15
mnn.seurat <- Seurat::RunUMAP(mnn.seurat, reduction = 'mnn', dims = 1:30)
#r1 <- Seurat::DimPlot(mnn.seurat, reduction = 'umap', group.by = 'tech')
#r2 <- Seurat::DimPlot(mnn.seurat, reduction = 'umap', group.by = 'celltype')
#table(Idents(mnn.seurat), mnn.seurat$tech)


# ADJUSTED RAND INDEX - BATCH MIXING
ari <- dplyr::select(mnn.seurat[[]], tech, seurat_clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c('celseq', 'celseq2', 'fluidigmc1', 'smartseq2', 'indrop'), to = c(0, 1, 2, 3, 4))
pdfCluster::adj.rand.index(as.numeric(ari$tech), as.numeric(ari$seurat_clusters))       #0.0169149


# ADJUSTED RAND INDEX - CELL PURITY
ari <- dplyr::select(mnn.seurat[[]], celltype, seurat_clusters)
ari$celltype <- plyr::mapvalues(ari$celltype,
                                from = c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma',
                                'macrophage', 'mast', 'quiescent_stellate', 'schwann'),
                                to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$celltype), as.numeric(ari$seurat_clusters))       #0.6380783
#resolution = 0.15, ari = 0.006628357, 0.9465943


# AVERAGE SILHOUETTE WIDTH (ASW)
#library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = mnn.seurat@reductions[['mnn']])[, 1:30])
sil.tech <- silhouette(x = as.numeric(x = as.factor(x = mnn.seurat$tech)), dist = dist.matrix)
mnn.seurat$sil.tech <- sil.tech[, 3]
sil.celltype <- silhouette(x = as.numeric(x = as.factor(x = mnn.seurat$celltype)), dist = dist.matrix)
mnn.seurat$sil.celltype <- sil.celltype[, 3]
mean(mnn.seurat$sil.tech)           #-0.06888286
mean(mnn.seurat$sil.celltype)       #0.3192418

####################### LIGER - SEURAT OBJECT ##################################

liger <- Seurat::NormalizeData(panc8)
liger <- Seurat::FindVariableFeatures(liger)
var.liger <- Seurat::SelectIntegrationFeatures(Seurat::SplitObject(liger, split.by = 'tech'), nfeatures = 2000, fvf.nfeatures = 2000, selection.method = 'vst')
liger <- Seurat::ScaleData(liger, split.by = 'tech', features = var.liger)

# RunOptimizeALS() and RunQuantileNorm() also found in the 'seurat-wrappers' github for Seurat
liger <- RunOptimizeALS(liger, k = 30, split.by = 'tech', do.center = FALSE)
liger <- RunQuantileNorm(liger, split.by = 'tech')
liger <- Seurat::FindNeighbors(liger, reduction = 'iNMF', dims = 1:30)
liger <- Seurat::FindClusters(liger)          # optimal resolution = 0.8 (default)
liger <- Seurat::RunUMAP(liger, reduction = 'iNMF', dims = 1:ncol(liger[['iNMF']]))
#s1 <- Seurat::DimPlot(liger, reduction = 'umap', group.by = 'tech')
#s2 <- Seurat::DimPlot(liger, reduction = 'umap', group.by = 'celltype')
#table(Idents(liger), liger$tech)


# ADJUSTED RAND INDEX - BATCH MIXING
ari <- dplyr::select(liger[[]], tech, seurat_clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c('celseq', 'celseq2', 'fluidigmc1', 'smartseq2', 'indrop'), to = c(0, 1, 2, 3, 4))
pdfCluster::adj.rand.index(as.numeric(ari$tech), as.numeric(ari$seurat_clusters))       #0.00571955


# ADJUSTED RAND INDEX - CELL PURITY
ari <- dplyr::select(liger[[]], celltype, seurat_clusters)
ari$celltype <- plyr::mapvalues(ari$celltype,
                                from = c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma',
                                'macrophage', 'mast', 'quiescent_stellate', 'schwann'),
                                to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$celltype), as.numeric(ari$seurat_clusters))       #0.617431


# AVERAGE SILHOUETTE WIDTH (ASW)
#library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = liger@reductions[['iNMF']])[, 1:30])
sil.tech <- silhouette(x = as.numeric(x = as.factor(x = liger$tech)), dist = dist.matrix)
liger$sil.tech <- sil.tech[, 3]
sil.celltype <- silhouette(x = as.numeric(x = as.factor(x = liger$celltype)), dist = dist.matrix)
liger$sil.celltype <- sil.celltype[, 3]
mean(liger$sil.tech)            #-0.1722261
mean(liger$sil.celltype)        #0.1042689

####################### LIGER - SINGLE CELL OBJECT #############################

plist <- list('celseq' = panc.list$celseq, 'celseq2' = panc.list$celseq2, 'fluidigmc1' = panc.list$fluidigmc1,
              'smartseq2' = panc.list$smartseq2, 'indrop' = panc.list$indrop)
liger.panc <- rliger::createLiger(sapply(plist, function(data) data[['RNA']]@counts[, colnames(data)]))
liger.panc <- rliger::normalize(liger.panc)
liger.panc <- rliger::selectGenes(liger.panc)
liger.panc <- rliger::scaleNotCenter(liger.panc)
liger.panc <- rliger::optimizeALS(liger.panc, k=20)
liger.panc <- rliger::quantile_norm(liger.panc)             # builds a similarity graph based on shared factor neighborhoods
liger.panc <- rliger::louvainCluster(liger.panc)            # identify clusters shared and align quantiles within each cluster
liger.panc <- rliger::runUMAP(liger.panc)
all.plots <- rliger::plotByDatasetAndCluster(liger.panc, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
all.plots[[1]]        # umap graph clustered by tech
all.plots[[2]]        # umap graph clusters
#plotClusterFactors(liger.panc, use.aligned = T)


# ADJUSTED RAND INDEX - BATCH MIXING
tech <- unlist(lapply(1:length(liger.panc@H), function(x) {rep(names(liger.panc@H)[x], nrow(liger.panc@H[[x]]))}))
ari <- data.frame('tech' = tech, 'clusters' = liger.panc@clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c('celseq', 'celseq2', 'fluidigmc1', 'smartseq2', 'indrop'), to = c(0, 1, 2, 3, 4))
pdfCluster::adj.rand.index(as.numeric(ari$tech), as.numeric(ari$clusters))        #0.01202221


# ADJUSTED RAND INDEX - CELL PURITY
ari <- data.frame('celltype' = panc8$celltype, 'clusters' = liger.panc@clusters)
ari$celltype <- plyr::mapvalues(ari$celltype,
                                from = c('acinar', 'activated_stellate', 'alpha', 'beta', 'delta', 'ductal', 'endothelial', 'epsilon', 'gamma',
                                'macrophage', 'mast', 'quiescent_stellate', 'schwann'),
                                to = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
pdfCluster::adj.rand.index(as.numeric(ari$celltype), as.numeric(ari$clusters))    #0.2387594

####################### ANALYSIS GRAPHS ########################################

ari.df <- data.frame(x= c(0.5642438, 0.9371157, 0.9465943, 0.617431), y = c(0.7949978, 0.993056178, 0.993371643, 0.99428045), names = c('Raw', 'Seurat CCA', 'fastMNN', 'LIGER'))
asw.df <- data.frame(x = c(0.2145548, 0.2974331, 0.3192418, 0.1042689), y = c(0.9755526, 1.1074509, 1.06888286, 1.1722261), names = c('Raw', 'Seurat CCA', 'fastMNN', 'LIGER'))

library(ggplot2)
library(ggrepel)
ggplot2::ggplot(ari.df, aes(x,y)) + geom_text_repel(aes(label = names), hjust = 1.5, vjust = -.5, direction = 'x') +
  geom_point(size = 5, color = c('darkslategray3', 'plum3', 'goldenrod2', 'coral')) +
  xlab('ARI celltype') + ylab('1 - ARI batch') +
  ggplot2::coord_cartesian(xlim = c(0.3, 1), ylim = c(0.5, 1.1)) + theme_bw() +
  theme(axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


ggplot2::ggplot(asw.df, aes(x,y)) + geom_text_repel(aes(label = names), hjust = 1, vjust = 1, direction = 'x') +
  geom_point(size = 5, color = c('darkslategray3', 'plum3', 'goldenrod2', 'coral')) +
  xlab('ASW celltype') + ylab('1 - ASW batch') +
  ggplot2::coord_cartesian(xlim = c(-0.1, 0.4), ylim = c(0.5, 1.3)) + theme_bw() +
  theme(axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
