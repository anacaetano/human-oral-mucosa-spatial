library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load data (example)
data_dir_A1 <- '/path/to/data/folder'
list.files(data_dir_A1)

data_dir_B1 <- '/path/to/data/folder'
list.files(data_dir_B1)

data_dir_C1 <- '/path/to/data/folder'
list.files(data_dir_C1)

data_dir_D1 <- '/path/to/data/folder'
list.files(data_dir_D1)

#Seurat workflow

A1_spatial <- Load10X_Spatial(data_dir_A1, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE)
plot1 <- VlnPlot(A1_spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(A1_spatial, features = "nCount_Spatial", pt.size.factor = 3.5) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

B1_spatial <- Load10X_Spatial(data_dir_B1, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE)
plot1 <- VlnPlot(B1_spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(B1_spatial, features = "nCount_Spatial", pt.size.factor = 2.5) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

C1_spatial <- Load10X_Spatial(data_dir_C1, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE)
plot1 <- VlnPlot(C1_spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(C1_spatial, features = "nCount_Spatial", pt.size.factor = 3) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

D1_spatial <- Load10X_Spatial(data_dir_D1, filename = "filtered_feature_bc_matrix.h5", assay = "Spatial", slice = "slice1", filter.matrix = TRUE, to.upper = FALSE)
plot1 <- VlnPlot(D1_spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(D1_spatial, features = "nCount_Spatial", pt.size.factor = 2) + theme(legend.position = "right")
wrap_plots(plot1, plot2)


A1_spatial <- SCTransform(A1_spatial, assay = "Spatial", verbose = FALSE)
B1_spatial <- SCTransform(B1_spatial, assay = "Spatial", verbose = FALSE)
C1_spatial <- SCTransform(C1_spatial, assay = "Spatial", verbose = FALSE)
D1_spatial <- SCTransform(D1_spatial, assay = "Spatial", verbose = FALSE)


#Dimensionality reduction, clustering and visualisation
A1_spatial <- RunPCA(A1_spatial, assay = "SCT", verbose = FALSE)
ElbowPlot(A1_spatial)
A1_spatial <- FindNeighbors(A1_spatial, reduction = "pca", dims = 1:10)
A1_spatial <- FindClusters(A1_spatial, verbose = FALSE)
A1_spatial <- RunUMAP(A1_spatial, reduction = "pca", dims = 1:10)

p1 <- DimPlot(A1_spatial, reduction = "umap", label = TRUE, cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink"))
p2 <- SpatialDimPlot(A1_spatial, label.size = 3, pt.size.factor = 2, cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink"))
p1 + p2

p1 <- DimPlot(A1_spatial, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(A1_spatial, label.size = 3, pt.size.factor = 2)
p1 + p2


B1_spatial <- RunPCA(B1_spatial, assay = "SCT", verbose = FALSE)
ElbowPlot(B1_spatial)
B1_spatial <- FindNeighbors(B1_spatial, reduction = "pca", dims = 1:10)
B1_spatial <- FindClusters(B1_spatial, verbose = FALSE)
B1_spatial <- RunUMAP(B1_spatial, reduction = "pca", dims = 1:10)

p1 <- DimPlot(B1_spatial, reduction = "umap", label = TRUE, cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink"))
p2 <- SpatialDimPlot(B1_spatial, label.size = 3, pt.size.factor = 3, cols = c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "pink"))
p1 + p2

C1_spatial <- RunPCA(C1_spatial, assay = "SCT", verbose = FALSE)
ElbowPlot(C1_spatial)
C1_spatial <- FindNeighbors(C1_spatial, reduction = "pca", dims = 1:7)
C1_spatial <- FindClusters(C1_spatial, verbose = FALSE)
C1_spatial <- RunUMAP(C1_spatial, reduction = "pca", dims = 1:7)

p1 <- DimPlot(C1_spatial, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(C1_spatial, label = TRUE, label.size = 3, pt.size.factor = 3.5)
p1 + p2


D1_spatial <- RunPCA(D1_spatial, assay = "SCT", verbose = FALSE)
ElbowPlot(D1_spatial)
D1_spatial <- FindNeighbors(D1_spatial, reduction = "pca", dims = 1:10)
D1_spatial <- FindClusters(D1_spatial, verbose = FALSE)
D1_spatial <- RunUMAP(D1_spatial, reduction = "pca", dims = 1:10)

p1 <- DimPlot(D1_spatial, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(D1_spatial, label = TRUE, label.size = 3)
p1 + p2


#Highlight regions
SpatialDimPlot(A1_spatial, cells.highlight = CellsByIdentities(object = A1_spatial, idents = 7), facet.highlight = TRUE, ncol = 3, pt.size.factor = 2.5,  cols.highlight = c("seagreen3", "grey50"))

SpatialDimPlot(B1_spatial, cells.highlight = CellsByIdentities(object = B1_spatial, idents = 8), facet.highlight = TRUE, ncol = 3, pt.size.factor = 3, cols.highlight = c("seagreen3", "grey50"))

SpatialDimPlot(C1_spatial, cells.highlight = CellsByIdentities(object = C1_spatial, idents = 4), facet.highlight = TRUE, ncol = 3.5, pt.size.factor = 2.5, cols.highlight = c("lightblue", "grey50"))

SpatialDimPlot(D1_spatial, cells.highlight = CellsByIdentities(object = D1_spatial, idents = 8), facet.highlight = TRUE, ncol = 3, cols.highlight = c("lightblue", "grey50"))

#Identification of spatial variable features
A1_spatial <- FindSpatiallyVariableFeatures(A1_spatial, assay = "SCT", features = VariableFeatures(A1_spatial)[1:1000], 
                                            selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(A1_spatial), selection.method = "markvariogram")
SpatialFeaturePlot(A1_spatial, features = top.features, ncol = 3, alpha = c(0.1, 1), pt.size.factor = 3.5)

head(top.features, n=20)

#Find spatial variables by region (alternative example)
##Region 0
de_markers <- FindMarkers(A1_spatial, ident.1 = 0)
SpatialFeaturePlot(object = A1_spatial, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3, pt.size.factor = 2.5)
head(de_markers, n=50)
write.table(de_markers, file = "A1-region0.csv", sep = ',', quote = F, col.names = NA)


#Integration analysis
#add single-cell rna seq data 
#The following dataset was used in the manuscript
(Caetano & Yianni et al. 2021) GSE152042

load("'/path/to/data/folder'")

perio.combined <- SCTransform(perio.combined, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:50)

Stromal<- SCTransform(Stromal, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:50)


#A1
A1_spatial <- SCTransform(A1_spatial, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
DimPlot(epithelial, group.by = "seurat_clusters", label = TRUE)

anchors <- FindTransferAnchors(reference = perio.combined, query = A1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = epithelial$seurat_clusters, prediction.assay = TRUE, 
                                  weight.reduction = A1_spatial[["pca"]], dims = 1:30, k.weight = 32)
A1_spatial[["predictions"]] <- predictions.assay

DefaultAssay(A1_spatial) <- "predictions"
SpatialFeaturePlot(A1_spatial, features = "5", ncol = 2, crop = TRUE)

A1_spatial <- FindSpatiallyVariableFeatures(A1_spatial, assay = "predictions", selection.method = "markvariogram", 
                                            features = rownames(A1_spatial), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(A1_spatial), 4)

SpatialPlot(object = A1_spatial, features = top.clusters, ncol = 2)

#B1
B1_spatial <- SCTransform(B1_spatial, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
DimPlot(perio.combined, group.by = "seurat_clusters", label = TRUE)

anchors <- FindTransferAnchors(reference = perio.combined, query = B1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Stromal$seurat_clusters, prediction.assay = TRUE, 
                                  weight.reduction = B1_spatial[["pca"]], dims = 1:30, k.weight = 32)
B1_spatial[["predictions"]] <- predictions.assay
DefaultAssay(A1_spatial) <- "predictions"
SpatialFeaturePlot(B1_spatial, features = "6", ncol = 2, crop = TRUE, pt.size.factor = 2)

B1_spatial <- FindSpatiallyVariableFeatures(B1_spatial, assay = "predictions", selection.method = "markvariogram", 
                                            features = rownames(A1_spatial), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(B1_spatial), 4)

SpatialPlot(object = B1_spatial, features = top.clusters, ncol = 2, pt.size.factor = 3)


#C1
C1_spatial <- SCTransform(C1_spatial, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
DimPlot(perio.combined, group.by = "seurat_clusters", label = TRUE)

anchors <- FindTransferAnchors(reference = perio.combined, query = C1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = perio.combined$seurat_clusters, prediction.assay = TRUE, 
                                  weight.reduction = C1_spatial[["pca"]], dims = 1:30, k.weight = 32)
C1_spatial[["predictions"]] <- predictions.assay

DefaultAssay(C1_spatial) <- "predictions"
SpatialFeaturePlot(C1_spatial, features = "15", ncol = 2, crop = TRUE, pt.size.factor = 3)

C1_spatial <- FindSpatiallyVariableFeatures(C1_spatial, assay = "predictions", selection.method = "markvariogram", 
                                            features = rownames(C1_spatial), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(C1_spatial), 4)

SpatialPlot(object = C1_spatial, features = top.clusters, ncol = 2, pt.size.factor = 3)


#D1
D1_spatial <- SCTransform(D1_spatial, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
DimPlot(perio.combined, group.by = "seurat_clusters", label = TRUE)

anchors <- FindTransferAnchors(reference = perio.combined, query = D1_spatial, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = perio.combined$seurat_clusters, prediction.assay = TRUE, 
                                  weight.reduction = D1_spatial[["pca"]], dims = 1:30, k.weight = 32)
D1_spatial[["predictions"]] <- predictions.assay

DefaultAssay(A1_spatial) <- "predictions"
SpatialFeaturePlot(D1_spatial, features = "15", ncol = 2, crop = TRUE)

D1_spatial <- FindSpatiallyVariableFeatures(D1_spatial, assay = "predictions", selection.method = "markvariogram", 
                                            features = rownames(D1_spatial), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(D1_spatial), 4)

SpatialPlot(object = D1_spatial, features = top.clusters, ncol = 2)


# single-cell data integration with merged consecutive sections (technical replicates)

#Integration analysis A1B1

A1B1.merge <- SCTransform(A1B1.merge, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)

DimPlot(perio.combined, group.by = "seurat_clusters", label = TRUE)

anchors <- FindTransferAnchors(reference = perio.combined, query = A1B1.merge, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = perio.combined$seurat_clusters, prediction.assay = TRUE, 
                                  weight.reduction = A1B1.merge[["pca"]], dims = 1:30, k.weight = 55)
A1B1.merge[["predictions"]] <- predictions.assay

DefaultAssay(A1B1.merge) <- "predictions"
SpatialFeaturePlot(A1B1.merge, features = "15", pt.size.factor = 3.5, ncol = 2, crop = TRUE)

A1B1.merge <- FindSpatiallyVariableFeatures(A1B1.merge, assay = "predictions", selection.method = "markvariogram", 
                                            features = rownames(A1B1.merge), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(A1B1.merge), 4)
SpatialPlot(object = A1B1.merge, features = top.clusters, ncol = 2, pt.size.factor = 3.5)

SpatialFeaturePlot(D1_spatial, features = c("8", "6"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))


#Integration analysis C1D1

C1D1.merge <- SCTransform(C1D1.merge, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)

DimPlot(perio.combined, group.by = "seurat_clusters", label = TRUE)

anchors <- FindTransferAnchors(reference = perio.combined, query = C1D1.merge, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = perio.combined$seurat_clusters, prediction.assay = TRUE, 
                                  weight.reduction = C1D1.merge[["pca"]], dims = 1:30, k.weight = 55)
C1D1.merge[["predictions"]] <- predictions.assay

DefaultAssay(C1D1.merge) <- "predictions"
SpatialFeaturePlot(C1D1.merge, features = "15", pt.size.factor = 3.5, ncol = 2, crop = TRUE)

C1D1.merge <- FindSpatiallyVariableFeatures(C1D1.merge, assay = "predictions", selection.method = "markvariogram", 
                                            features = rownames(C1D1.merge), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(A1B1.merge), 4)
SpatialPlot(object = A1B1.merge, features = top.clusters, ncol = 2, pt.size.factor = 3.5)

SpatialFeaturePlot(D1_spatial, features = c("8", "6"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha 



#BayesSpace
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(scuttle)
library(dplyr)

#Loading data
sce_A1 <- readVisium('/path/to/data/folder')

sce_B1 <- readVisium("/path/to/data/folder")

sce_C1 <- readVisium("/path/to/data/folder")

sce_D1 <- readVisium("/path/to/data/folder")

#Preprocessing data
set.seed(102)

sce_A1 <- spatialPreprocess(sce_A1, platform = "Visium", n.PCs = 10, n.HVGs = 200, log.normalize = TRUE)

sce_B1 <- spatialPreprocess(sce_B1, platform = "Visium", n.PCs = 10, n.HVGs = 200, log.normalize = TRUE)

sce_C1 <- spatialPreprocess(sce_C1, platform = "Visium", n.PCs = 7, n.HVGs = 200, log.normalize = TRUE)

sce_D1 <- spatialPreprocess(sce_D1, platform = "Visium", n.PCs = 7, n.HVGs = 200, log.normalize = TRUE)

#Clustering
## understand number of clusters
sce_A1 <- qTune(sce_A1, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce_A1)

sce_B1 <- qTune(sce_B1, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce_B1)

sce_C1 <- qTune(sce_C1, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce_C1)

sce_D1 <- qTune(sce_D1, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce_D1)

## clustering
set.seed(149)

sce_A1 <- spatialCluster(sce_A1, q=10, platform = "Visium", d=7, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE)
head(colData(sce_A1))


sce_B1 <- spatialCluster(sce_B1, q=8, platform = "Visium", d=7, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE)
head(colData(sce_B1))

sce_C1 <- spatialCluster(sce_C1, q=7, platform = "Visium", d=7, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE)
head(colData(sce_C1))

sce_D1 <- spatialCluster(sce_D1, q=7, platform = "Visium", d=7, init.method = "mclust", model = "t", gamma = 2, nrep = 10000, burn.in = 100, save.chain = TRUE)
head(colData(sce_D1))

#visualising spatial clusters

clusterPlot(sce_A1)

clusterPlot(sce_B1)

clusterPlot(sce_C1)

clusterPlot(sce_D1)

##customise colours
clusterPlot(sce_B1, palette=c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "darkorange", "seagreen", "seagreen3", "lightgoldenrod", "lightsteelblue"), color="black") +
  theme_bw() +
  xlab("Column") +
  ylab("Row") +
  labs(fill="BayesSpace\ncluster", title="Spatial clustering of A1")

# Enhanced resolution (example A1)
##The spatialEnhance() function will enhance the resolution of the principal components, and add these PCs as well as predicted cluster labels at subspot resolution to a new SingleCellExperiment.

sce_A1.enhanced <- spatialEnhance(sce_A1, q=10, platform="Visium", d=7,
                                  model="t", gamma=2,
                                  jitter_prior=0.3, jitter_scale=3.5,
                                  nrep=10000, burn.in=100,
                                  save.chain=TRUE)

head(colData(sce_A1.enhanced))

clusterPlot(sce_A1.enhanced)

clusterPlot(sce_A1.enhanced, palette=c("steelblue", "lightskyblue1", "grey", "navy", "darkcyan", "cadetblue", "seagreen", "seagreen3", "tan1", "lightsteelblue")) +
  theme_bw() +
  xlab("Column") +
  ylab("Row") +
  labs(fill="BayesSpace\ncluster", title="Spatial clustering of enhanced A1")

clusterPlot(sce_A1.enhanced) +
  theme_bw() +
  xlab("Column") +
  ylab("Row") +
  scale_fill_brewer('RdYlBu') +
  labs(fill="BayesSpace\ncluster", title="Spatial clustering of enhanced A1")


## Convert SCE to Seurat object and use BayesSpace cluster as identifier

sobj <- Seurat::CreateSeuratObject(counts=logcounts(sceA1.enhanced),
                                   assay='Spatial',
                                   meta.data=as.data.frame(colData(sceA1.enhanced)))

sobj <- Seurat::SetIdent(sobj, value = "spatial.cluster")

sobj@assays$Spatial@scale.data <-
  sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t



## Select top n markers from each cluster (by log fold change)
top_markers <- Seurat::FindAllMarkers(sobj, assay='Spatial', slot='data',
                                      group.by='spatial.cluster',
                                      only.pos=TRUE) %>% 
  group_by(cluster) %>% 
  top_n(5, avg_log2FC)


Seurat::DoHeatmap(sobj, features = top_markers$gene, slot='scale.data',
                  group.by = "spatial.cluster",
                  angle=0, size=4, label = FALSE, raster=FALSE) + 
  guides(scale = "none")


#Enhanced gene expression

markers <- list()
markers[["Basal keratinocyte"]] <- c("KRT5", "KRT14", "COL7A1", "CXCL14")
markers[["Suprabasal keratinocyte"]] <- c("SBSN", "SPRR1B")
markers[["Myeloid"]] <- c("LYZ")
markers[["T-cell"]] <- c("CD2", "CD3D", "CD3E", "CD3G", "CD7")
markers[["Melanocyte"]] <- c("MLANA")
markers[["Endothelial"]] <- c("TFF3", "CLDN5", "VWF")
markers[["Stromal"]] <- c("COL1A1", "COL3A1")
markers[["Regulatory T"]] <- c("FOXP3", "LAIR2")
markers[["gamma delta T"]] <- c("FXYD2", "KLRC4")
markers[["MAIT"]] <- c("SLC4A10", "IL23R", "RORC")
markers[["NK cells"]] <- c("SH2D1B", "KLRF1", "SPTSSB")
markers[["B/plasma"]] <- c("POU2AF1", "IGLC3")
markers[["Mast"]] <- c("MSA42", "P2RX1", "HPGDS")
markers[["DC"]] <- c("PPY", "XCR1", "MSR1", "FPR1")
markers[["aDC"]] <- c("LAMP3", "MERG", "HMSD")
markers[["pDC"]] <- c("PACSIN1", "SMIM5")
markers[["Langerhans"]] <- c("RASAL1", "FCEBP", "CD207")
markers[["Macrophage"]] <- c("CCL18", "FOLR2", "MMP9")

sce_A1.enhanced <- enhanceFeatures(sce_A1.enhanced, sce_A1,
                                   model="xgboost",
                                   feature_names=purrr::reduce(markers, c),
                                   nrounds=0)

sum_counts <- function(sce_A1, features) {
  if (length(features) > 1) {
    colSums(logcounts(sce_A1)[features, ])
  } else {
    logcounts(sce_A1)[features, ]
  }
}


#Giotto ligand-receptor analysis
library(Giotto)
library(janitor)

results_folder = '/path/to/data/folder'
python_path = '/path/to/data/folder'

## create instructions
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = python_path)

## provide path to visium folder
data_path = '/path/to/data/folder'


## Create Giotto object & Process Data
## directly from visium folder
A1 = createGiottoVisiumObject(visium_dir = data_path, expr_data = 'raw',
                              png_name = 'tissue_lowres_image.png',
                              gene_column_index = 2, instructions = instrs)

## update and align background image
# problem: image is not perfectly aligned
spatPlot(gobject = A1, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = '2_a_spatplot_image'))

showGiottoImageNames(A1)

A1 = updateGiottoImage(A1, image_name = 'image',
                       xmax_adj = 800, xmin_adj =600,
                       ymax_adj = 600, ymin_adj = 600)

spatPlot(gobject = A1, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = '2_b_spatplot_image_adjusted'))

## check metadata
pDataDT(A1)

## compare in tissue with provided jpg
spatPlot(gobject = A1, cell_color = 'in_tissue', point_size = 2,
         cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
         save_param = list(save_name = '2_c_in_tissue'))

## subset on spots that were covered by tissue
metadata = pDataDT(A1)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
A1 = subsetGiotto(A1, cell_ids = in_tissue_barcodes)

## filter
A1 <- filterGiotto(gobject = A1,
                   expression_threshold = 0.5,
                   gene_det_in_min_cells = 50,
                   min_det_genes_per_cell = 10,
                   expression_values = c('raw'),
                   verbose = T)

## normalize
A1 <- normalizeGiotto(gobject = A1, scalefactor = 6000, verbose = T)

## add gene & cell statistics
A1 <- addStatistics(gobject = A1)

## visualize
spatPlot2D(gobject = A1, show_image = T, point_alpha = 0.7,
           save_param = list(save_name = '2_d_spatial_locations'))

spatPlot2D(gobject = A1, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_genes', color_as_factor = F,
           save_param = list(save_name = '2_e_nr_genes'))


## Dimension Reduction
## highly variable genes (HVG)
A1 <- calculateHVG(gobject = A1,
                   save_param = list(save_name = '3_a_HVGplot'))

## run PCA on expression values (default)
A1 <- runPCA(gobject = A1, center = TRUE, scale_unit = T)
screePlot(A1, ncp = 30, save_param = list(save_name = '3_b_screeplot'))

plotPCA(gobject =A1,
        save_param = list(save_name = '3_c_PCA_reduction'))

## run UMAP and tSNE on PCA space (default)
##UMAP
A1 <- runUMAP(A1, dimensions_to_use = 1:2)
plotUMAP(gobject = A1,
         save_param = list(save_name = '3_d_UMAP_reduction'))
##tSNE
A1 <- runtSNE(A1, dimensions_to_use = 1:2)
plotTSNE(gobject = A1,
         save_param = list(save_name = '3_e_tSNE_reduction'))

##Clustering

## sNN network (default)
A1 <- createNearestNetwork(gobject = A1, dimensions_to_use = 1:7, k = 15)
## Leiden clustering
A1 <- doLeidenCluster(gobject = A1, resolution = 0.35)
plotUMAP(gobject = A1,
         cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,
         save_param = list(save_name = '4_a_UMAP_leiden'))

A1 <- doLouvainCluster(gobject = A1, resolution = 3)
plotUMAP(gobject = A1,
         cell_color = 'louvain_clus', show_NN_network = T, point_size = 2.5,
         save_param = list(save_name = '4_a_UMAP_louvain'))


##Spatial Grid
A1 <- createSpatialGrid(gobject = A1,
                        sdimx_stepsize = 400,
                        sdimy_stepsize = 400,
                        minimum_padding = 0)
spatPlot(A1, cell_color = 'cell_types', show_grid = T,
         grid_color = 'red', spatial_grid_name = 'spatial_grid', 
         save_param = c(save_name = '8_grid'))


##Spatial Network
## Delaunay network: stats + creation
plotStatDelaunayNetwork(gobject = A1, maximum_distance = 400, 
                        save_param = c(save_name = '9_a_delaunay_network'))

A1 = createSpatialNetwork(gobject = A1, minimum_k = 0)
showNetworks(A1)
spatPlot(gobject = A1, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         save_param = c(save_name = '9_b_delaunay_network'))

##Spatial Genes
## kmeans binarization
kmtest = binSpect(A1)
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = kmtest$genes[1:20], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_km'))
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = kmtest$genes[21:40], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_km2'))
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = kmtest$genes[41:60], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_km3'))

##Rank Binarisation
ranktest = binSpect(A1, bin_method = 'rank')
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = ranktest$genes[1:20], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_rank'))
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = ranktest$genes[21:40], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_rank2'))
# spatGenePlot(A1, expression_values = 'scaled',
#              genes = ranktest$genes[41:60], cow_n_col = 4, point_size = 1.5,
#              save_param = c(save_name = 'spatial_genes_rank3'))



## Cell Neighborhood analysis
cell_proximities = cellProximityEnrichment(gobject = A1,
                                           cluster_column = 'cell_types',
                                           spatial_network_name = 'Delaunay_network',
                                           adjust_method = 'fdr',
                                           number_of_simulations = 1000)

## barplot
cellProximityBarplot(gobject = A1,
                     CPscore = cell_proximities, min_orig_ints = 5, min_sim_ints = 5, 
                     save_param = c(save_name = '12_a_barplot_cell_cell_enrichment'))

## heatmap
cellProximityHeatmap(gobject = A1, CPscore = cell_proximities, order_cell_types = T, scale = T,
                     color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),
                     save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'))

## network
cellProximityNetwork(gobject = A1, CPscore = cell_proximities, remove_self_edges = T,
                     only_show_enrichment_edges = F,
                     save_param = c(save_name = '12_c_network_cell_cell_enrichment'))

## network with self-edges
cellProximityNetwork(gobject = A1, CPscore = cell_proximities,
                     remove_self_edges = F, self_loop_strength = 0.3,
                     only_show_enrichment_edges = F,
                     rescale_edge_weights = T,
                     node_size = 8,
                     edge_weight_range_depletion = c(1, 2),
                     edge_weight_range_enrichment = c(2,5),
                     save_param = c(save_name = '12_d_network_cell_cell_enrichment_self',
                                    base_height = 5, base_width = 5, save_format = 'pdf'))



## select top 25th highest expressing genes
gene_metadata = fDataDT(A1)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr_det)

quantile(gene_metadata$mean_expr_det)
high_expressed_genes = gene_metadata[mean_expr_det > 1.31]$gene_ID


# LR expression
# LR activity changes

LR_data = data.table::fread ('/path/to/data/folder')
LR_data

LR_data[, ligand_det := ifelse(ligand %in% A1@gene_ID, T, F)]
LR_data[, receptor_det := ifelse(receptor %in% A1@gene_ID, T, F)]
LR_data_det = LR_data[ligand_det == T & receptor_det == T]
select_ligands = LR_data_det$ligand
select_receptors = LR_data_det$receptor


## get statistical significance of gene pair expression changes based on expression ##
expr_only_scores = exprCellCellcom(gobject = A1,
                                   cluster_column = 'cell_types', 
                                   random_iter = 500,
                                   gene_set_1 = select_ligands,
                                   gene_set_2 = select_receptors, 
                                   verbose = FALSE)

## get statistical significance of gene pair expression changes upon cell-cell interaction
spatial_all_scores = spatCellCellcom(A1,
                                     spatial_network_name = 'Delaunay_network',
                                     cluster_column = 'cell_types', 
                                     random_iter = 500,
                                     gene_set_1 = select_ligands,
                                     gene_set_2 = select_receptors,
                                     adjust_method = 'fdr',
                                     do_parallel = T,
                                     cores = 4,
                                     verbose = 'none')


## select top LR ##
selected_spat = spatial_all_scores[p.adj <= 0.5 & abs(log2fc) > 0.1 & lig_nr >= 2 & rec_nr >= 2]
data.table::setorder(selected_spat, -PI)

top_LR_ints = unique(selected_spat[order(-abs(PI))]$LR_comb)[1:33]
top_LR_cell_ints = unique(selected_spat[order(-abs(PI))]$LR_cell_comb)[1:33]

plotCCcomDotplot(gobject = A1,
                 comScores = spatial_all_scores,
                 selected_LR = top_LR_ints,
                 selected_cell_LR = top_LR_cell_ints,
                 cluster_on = 'PI',
                 save_param = c(save_name = '14_a_communication_dotplot', save_format = 'pdf'))


