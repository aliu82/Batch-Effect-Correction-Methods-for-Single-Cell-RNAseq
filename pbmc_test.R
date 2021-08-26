#library(dplyr)
#library(Seurat)
#library(patchwork)
#library(cowplot)
#library(Matrix)
#library(rliger)
#library(SingleCellExperiment)
#library(batchelor)
#library(scran)
#library(devtools)

# first run test - after learning Seurat through SatijaLab, testing different approaches of methodologies

####################### PBMC DATA - SEURAT OBJECT ##############################

pbmc.data <- Seurat::Read10X(data.dir = 'pbmc/')
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, project = 'pbmc3k', min.cells = 3, min.features = 200)
#pbmc <- Seurat::VlnPlot(pbmc, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
pbmc <- Seurat::NormalizeData(pbmc)
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)#'PPBP', 'LYZ', 'S100A9', 'FTL', 'GNLY', 'IGLL5'
pbmc <- Seurat::ScaleData(pbmc, features = rownames(pbmc))
pbmc <- Seurat::RunPCA(pbmc, pc.genes = pbmc@var.genes)
pbmc <- Seurat::FindNeighbors(pbmc)
pbmc <- Seurat::RunUMAP(pbmc, dim = 1:10)
pbmc <- Seurat::FindClusters(pbmc, reduction='pca')
pbmc <- Seurat::RunTSNE(pbmc, reduction='pca')
Seurat::DimPlot(pbmc, reduction = 'tsne')                                   #VISUALIZATION WITHOUT BATCH CORRECTION

####################### PBMC DATA - SINGLE CELL OBJECT #########################

test <- Seurat::CreateSeuratObject(counts = pbmc.data, project = 'pbmc3k', min.cells = 3, min.features = 200)
test1 <- subset(x = test, cells = colnames(x = pbmc)[1:1000])
test2 <- subset(x = test, cells = colnames(x = pbmc)[1001:2000])
test1 <- Seurat::as.SingleCellExperiment(test1)
test2 <- Seurat::as.SingleCellExperiment(test2)
universe <- intersect(rownames(test1), rownames(test2))
test1 <- test1[universe,]                                     #subset the SingleCellExperiment object
test2 <- test2[universe,]


var1 <- scran::modelGeneVar(Seurat::as.SingleCellExperiment(pbmc1))
var2 <- scran::modelGeneVar(Seurat::as.SingleCellExperiment(pbmc2))
combined.var <- scran::combineVar(var1, var2)
chosen.hvgs <- combined.var$bio > 0


rescaled <- batchelor::multiBatchNorm(test1, test2)    #rescale the batch to adjust for systematic differences in SingleCellExperiment
test1 <- rescaled[[1]]
test2 <- rescaled[[2]]
test1$batch <- 'test1'
test2$batch <- 'test2'
uncorrected <- cbind(test1, test2)
set.seed(0010101010)
uncorrected <- scater::runPCA(uncorrected, subset_row = chosen.hvgs, BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
snn.uncorrected <- scran::buildSNNGraph(uncorrected, use.dimred = 'PCA')
clusters <- igraph::cluster_walktrap(snn.uncorrected)$membership
tab <- table(Cluster = clusters, Batch = uncorrected$batch)
uncorrected <- scater::runTSNE(uncorrected, dimred = 'PCA')
scater::plotTSNE(uncorrected, colour_by = 'batch')

####################### SEURAT CCA #############################################

cca.pbmc <- Seurat::RunCCA(object1 = pbmc1, object2 = pbmc2)
DimPlot(cca.pbmc, reduction = 'cca')                                #VISUALIZATION WITH CCA
#VlnPlot(cca.pbmc, features = c('PPBP', 'LYZ', 'S100A9', 'FTL', 'GNLY', 'IGLL5'))
#DimHeatmap(object = cca.pbmc, reduction = 'cca', cells = 500, dims = 1:9, balanced = TRUE)
#FeaturePlot(object = cca.pbmc, features = rownames(cca.pbmc)[1:10], reduction = 'cca')

####################### FAST MUTUAL NEAREST NEIGHBORS ##########################

mnn.pbmc <- batchelor::fastMNN(Seurat::as.SingleCellExperiment(pbmc1), Seurat::as.SingleCellExperiment(pbmc2), subset.row=chosen.hvgs, BSPARAM=BiocSingular::RandomParam(deferred=TRUE))
snn.gr <- scran::buildSNNGraph(mnn.pbmc, use.dimred = 'corrected')
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster = clusters.mnn, Batch = mnn.pbmc$batch)
norm <- SingleCellExperiment::normcounts(mnn.pbmc)
chi.prop <- colSums(tab.mnn)/sum(tab.mnn)
chi.results <- apply(tab.mnn, 1, FUN = chisq.test, p = chi.prop)
p.values <- vapply(chi.results, '[[', i = 'p.value', 0)
tab.mnn <- unclass(tab.mnn)
norm <- scuttle::normalizeCounts(tab.mnn, pseudo.count = 10)
rv <- rowVars(norm)
#DataFrame(Batch = tab.mnn, var = rv)[order(rv, decreasing = TRUE),]
mnn.pbmc <- scater::runTSNE(mnn.pbmc, dimred = 'corrected')
mnn.pbmc$batch <- factor(mnn.pbmc$batch)
scater::plotTSNE(mnn.pbmc, colour_by = 'batch')

####################### LIGER ##################################################

liger.pbmc <- rliger::seuratToLiger(as.matrix(list('pbmc1' = pbmc1, 'pbmc2' = pbmc2)))
liger.pbmc <- rliger::normalize(liger.pbmc)
liger.pbmc <- rliger::selectGenes(liger.pbmc)
liger.pbmc <- rliger::scaleNotCenter(liger.pbmc)
liger.pbmc <- rliger::optimizeALS(liger.pbmc, k=20)
liger.pbmc <- rliger::quantile_norm(liger.pbmc)             #builds a similarity graph based on shared factor neighborhoods
liger.pbmc <- rliger::louvainCluster(liger.pbmc)            #identify clusters shared and align quantiles within each cluster
liger.pbmc <- rliger::runUMAP(liger.pbmc)
#liger.pbmc <- rliger::runTSNE(liger.pbmc, use.raw = T)
#p <- rliger::plotByDatasetAndCluster(liger.pbmc, axis.labels = c('UMAP', 'UMAP2'), return.plots = T)
#print(p[[1]])
rliger::plotClusterProportions(liger.pbmc)
#rliger::plotClusterFactors(liger.pbmc, use.aligned = T)
markers <- rliger::getFactorMarkers(liger.pbmc, dataset1 = 'pbmc1', dataset2 = 'pbmc2', num.genes = 10)


liger.pbmc <- rliger::plotGeneLoadings(liger.pbmc, do.spec.plot=FALSE, return.plots=TRUE)
rliger::plotGene(liger.pbmc, gene = 'PPBP')
rliger::plotGeneViolin(liger.pbmc, gene = 'PPBP')

