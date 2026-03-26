# Load necessary libraries
library(spacexr)
library(pheatmap)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(SpatialExperiment)
library(Seurat)
library(SeuratWrappers)
library(Banksy)
library(ggplot2)

# Load reference data
counts <- ReadMtx(
  '/mnt/bioinfo_sharing/sharing/astrid/BrCa_ref/count_matrix_sparse.mtx',
  '/mnt/bioinfo_sharing/sharing/astrid/BrCa_ref/count_matrix_barcodes.tsv',
  '/mnt/bioinfo_sharing/sharing/astrid/BrCa_ref/count_matrix_genes.tsv',
  cell.column = 1,
  feature.column = 1,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)

meta <- read.csv(file = '/mnt/bioinfo_sharing/sharing/astrid/BrCa_ref/metadata.csv', header = TRUE)

ref <- CreateSeuratObject(
  counts,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = meta,
  project = "BrCaRef"
)
#normalization
ref <- NormalizeData(ref)
ref <- ScaleData(ref)

#runPCA,findNeighbours, runUMAP
# Perform PCA on the preprocessed reference
ref <- FindVariableFeatures(ref)
ref <- RunPCA(ref, assay = "RNA", npcs = 30, verbose = FALSE)
#Check elbow plot to determine the number of significant PCs
ElbowPlot(ref, ndims = 30, reduction = 'pca')
# Perform UMAP on the PCA results
ref <- FindNeighbors(ref)
ref <- RunUMAP(ref, reduction = "pca", dims = 1:30)
#Visualize the UMAP results
umap <- DimPlot(ref, reduction = "umap", label = TRUE) + ggtitle("UMAP of Spatial Transcriptomics Data")
umap

DimPlot(ref, reduction = 'umap')
DimPlot(ref, reduction = 'umap', group.by = "celltype_major")
DimPlot(ref, reduction = 'umap', group.by = "celltype_minor", label = TRUE) + NoLegend()
DimPlot(ref, reduction = 'umap', group.by = "celltype_subset", label = TRUE) + NoLegend()


Idents(ref) <- "celltype_minor"

# Retrieve raw count matrix representing gene expression per cell
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$celltype_minor)
names(cluster) <- colnames(ref)
# Obtain total UMI counts per cell to account for sequencing depth variation
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
#build reference object
reference <- Reference(counts, cluster, nUMI)



#prepare the spatial data to use with RCTD
counts <- visium[["Visium"]]$counts
coords <- GetTissueCoordinates(visium)

#for segmented 
rownames(coords) <- coords$cell
coords <- coords[,c("x","y")]

coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))


# Create RCTD object on 8 cores
RCTD <- create.RCTD(query, reference, max_cores = 8)

#segmented
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

#add results back to visium
visium <- AddMetaData(visium, metadata = RCTD@results$weights)
dim(RCTD@spatialRNA@counts)
head(colnames(RCTD@spatialRNA@counts))
ref_cell_types <- colnames(RCTD@results$weights)

#Visualize the results for the cell types identified by RCTD
celltype_plots <- list()
for (x in ref_cell_types){
  p <- SpatialFeaturePlot(xenium_sig, features = x, image.alpha = 0.5)
  celltype_plots <- c(celltype_plots, p)
}
names(celltype_plots) <- ref_cell_types
celltype_plots



#visium visualizations
ImageDimPlot(visium, group.by = "first_type",cols = celltype_colours_simplified, alpha = 0.3)
celltype_colours_simplified <- c("B cells Memory" = "red", "B cells Naive" = "red", "Plasmablasts" = "red",
                                 "T cells CD4+" = "magenta", "T cells CD8+" = "magenta", "Cycling T-cells" = "magenta",
                                 "NK cells" = "magenta2", "NKT cells" = "magenta2",
                                 "Macrophage" = "purple", "Monocyte" = "purple", "Cycling_Myeloid" = "purple", "DCs" = "purple",
                                 "Cancer Basal SC" = "lightblue", "Cancer Cycling", "lightblue", "Cancer Her2 SC" = "lightblue",
                                 "Cancer LumA SC" = "lightblue", "Cancer LumB SC" = "lightblue",
                                 "Luminal Progenitors" = "lightgreen", "Mature Luminal" = "lightgreen", "Myoepithelial" = "lightgreen",
                                 "CAFs MSC iCAF-like" = "orange", "CAFs myCAF-like" = "orange",
                                 "Endothelial ACKR1" = "gold", "Endothelial CXCL12" = "gold", "Endothelial Lymphatic LYVE1" = "gold",
                                 "Endothelial RGS5" = "gold",
                                 "Cycling PVL" = "yellow", "PVL Differentiated" = "yellow", "PVL Immature" = "yellow")               

summary(as.factor(visium$spot_class))
VlnPlot(visium, features = "nCount_Spatial", group.by = "spot_class", pt.size = 0)




