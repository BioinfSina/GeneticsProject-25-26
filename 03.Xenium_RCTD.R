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


# Load reference data.
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
# Perform PCA on the preprocessed data
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
counts <- xenium[["Xenium"]]$counts
coords <- GetTissueCoordinates(xenium)

#for Xenium
rownames(coords) <- coords$cell
coords <- coords[,c("x","y")]

coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))


# Create RCTD object on 8 cores
RCTD <- create.RCTD(query, reference, max_cores = 8)

#segmented
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

#add results back to xenium
xenium <- AddMetaData(xenium, metadata = RCTD@results$weights)
dim(RCTD@spatialRNA@counts)
head(colnames(RCTD@spatialRNA@counts))
ref_cell_types <- colnames(RCTD@results$weights)


celltype_cols <- colnames(xenium@meta.data)[13:41]

xenium$first_type <- apply(
  xenium@meta.data[, celltype_cols],
  1,
  function(x) names(which.max(x))
)


#Visualize the results for the cell types identified by RCTD
xenium <- AddMetaData(xenium, metadata = RCTD@results$results_df)
ImageDimPlot(xenium, group.by = "first_type",cols = celltype_colours_simplified, alpha = 0.3)
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

ImageDimPlot(xenium)

#backup adding results back to xenium
results_df <- RCTD@results$results_df
results_df <- results_df[rownames(xenium@meta.data), ]
xenium@meta.data <- cbind(xenium@meta.data, results_df)
summary(as.factor(xenium$spot_class))


#visualize spot class in different subsets to justify subsetting
#subsetting is based on nfeatures(genes detected in xenium)
xenium2 <- subset(
  xenium1,
  nFeature_Xenium >= 50
)
summary(as.factor(xenium2$spot_class))
xenium3 <- subset(
  xenium2,
  nFeature_Xenium >= 100
)
summary(as.factor(xenium3$spot_class))
xenium4 <- subset(
  xenium3,
  nFeature_Xenium >= 150
)
summary(as.factor(xenium4$spot_class))
summary(as.factor(xenium1$spot_class))

spot_colours <- c( "NA" = "magenta","singlet" = "purple", 
                   "reject" = "lightblue", 
                   "doublet_uncertain" = "lightgreen", 
                   "doublet_certain" = "gold"
)               

#plots
plot_spotclass_pie <- function(obj, title){
  
  df <- as.data.frame(table(obj$spot_class, useNA = "ifany"))
  colnames(df) <- c("spot_class","count")
  
  df$spot_class <- as.character(df$spot_class)
  df$spot_class[is.na(df$spot_class)] <- "NA"
  
  df$percent <- 100 * df$count / sum(df$count)
  
  df$label <- paste0(df$spot_class,"\n",
                     df$count," (",round(df$percent,1),"%)")
  
  p <- ggplot(df, aes(x="", y=count, fill=spot_class)) +
    
    geom_bar(stat="identity", width=1, color="white") +
    
    coord_polar("y") +
    
    geom_text(
      data = subset(df, percent >= 3),
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    
    
    geom_label_repel(
      data = subset(df, percent < 3),
      aes(label = label),
      position = position_stack(vjust = 0.5),
      show.legend = FALSE
    ) +
    
    scale_fill_manual(values = spot_colours) +
    
    theme_void() +
    labs(title = title)
  
  return(p)
}
p1 <- plot_spotclass_pie(xenium1, "Xenium")
p2 <- plot_spotclass_pie(xenium2, "Xenium nFeature_Xenium >= 50")
p3 <- plot_spotclass_pie(xenium3, "Xenium nFeature_Xenium >= 100")
p4 <- plot_spotclass_pie(xenium4, "Xenium nFeature_Xenium >= 150")

(p1 | p2)/
  (p3 | p4)


#visualize by spot class
ImageDimPlot(
  xenium,
  cells = reject_cells,
  group.by = "spot_class",
  cols = "polychrome",
  dark.background = TRUE
)





