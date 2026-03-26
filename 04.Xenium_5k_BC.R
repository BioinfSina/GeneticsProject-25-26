#load necessary libraries
library(Seurat)
library(SeuratWrappers)
library(Banksy)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(arrow)
library(leidenbase)
library(Nebulosa)


xenium1 <- LoadXenium("~/practical/FFPE_5k", fov = "fov", segmentations = "cell")
xenium <- NormalizeData(xenium)
#apply the normalized factors back to data
xenium <- ScaleData(xenium)

xenium <- FindVariableFeatures(xenium)

#dimensional reduction
xenium <- RunPCA(xenium, npcs = 30, features = VariableFeatures(xenium))
#Check elbow plot to determine the number of significant PCs
#as transcripts are sparse in xenium, usually not so many PCs
# elbow point is when the variance explained sudden drops, if following PC values show relatively 
# low values, then PC is probably saturated and could be considered abolished.
ElbowPlot(xenium, ndims = 30, reduction = 'pca')
xenium <- RunUMAP(xenium, dims = 1:30)

DimPlot(xenium1, reduction = "umap", label = TRUE, cols = "brilliant_purple") +
  ggtitle("Xenium") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


#clustering
xenium <- FindNeighbors(xenium, dims = 1:30, reduction = "pca")
xenium <- FindClusters(xenium, resolution = 0.5, algorithm = 4, cluster.name="clusters_0.5")
DimPlot(xenium, reduction = "umap", group.by = "clusters_0.5", label = TRUE) +
  labs(title = "UMAP of spatial clusters_0.5")
ImageDimPlot(xenium, group.by = "clusters_0.5") +
  labs(title = "Spatial clusters on image data_0.5")


#resolution for xenium could be smaller since trasncripts are sparse.
xenium <- FindClusters(xenium, resolution = 0.3, algorithm = 4, cluster.name="clusters_0.3")
#visualization
DimPlot(xenium, reduction = "umap", group.by = "clusters_0.3", label = TRUE) +
  labs(title = "UMAP of spatial clusters_0.3")
ImageDimPlot(xenium, group.by = "clusters_0.3") +
  labs(title = "Spatial clusters on image data_0.3")

# Show only cluster 1 on the spatial image
ImageDimPlot(xenium, cells = WhichCells(xenium, idents = 7))


#cell type proportion barplot
prop_cluster_x <- prop.table(xcell_cluster, 1)
df_x <- as.data.frame(prop_cluster_x)
celltypes_x <- unique(as.character(df_x$Var2))
kelly_expanded <- colorRampPalette(pals::kelly())(30)
pals_x <- setNames(kelly_expanded, c(celltypes_x, "Other"))

df_x$fill_group <- ifelse(df_x$Freq >= 0.02, as.character(df_x$Var2), "Other")

celltypes_x <- unique(as.character(df_x$Var2))

#visualization
ggplot(df_x, aes(x = Var1, y = Freq, fill = fill_group)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.1) + 
  scale_fill_manual(values = pals_x) +
  labs(x = "Xenium BANKSY Cluster", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



#find marker genes
#logfc.threshold could be smaller-xenium is targeted panel + log-normalized-fc already small
#logfc.threshold= 0.1
#test.use = "wilcox" is stable and robust
markers <- FindAllMarkers(xenium, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, pvalue.cutoff = 0.05)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(1, avg_log2FC) %>%
  ungroup()
FeaturePlot(xenium, features = top_markers$gene, ncol = 5) +
  plot_annotation(title = "Top marker genes for each cluster")
ImageFeaturePlot(xenium, features = top_markers$gene) +
  plot_annotation(title = "Top marker genes for each cluster on image data")



#visualizing single gene/double expression distribution
ImageFeaturePlot(
  xenium,
  features = c("COL4A1"),size = 1.0, dark.background = FALSE)
ImageFeaturePlot(xenium, features = c("SOSTDC1", "AGR3"), 
                 blend = TRUE, size = 0.8, dark.background = TRUE)[c(3, 4)]

p <- ImageFeaturePlot(xenium, features = c("ESR1",  "COL4A1"), 
                      cols = c("grey90","red","black"),
                      blend = TRUE, size = 1.3, dark.background = FALSE)
p[[3]]+ coord_cartesian(
  xlim = c(3500, 6500),
  ylim = c(2000, 6500)
)  + p[[4]] + plot_layout(widths = c(4, 1))




#adding gene module score(based on marker genes selected from visium)
immune <- c("CXCL13","CXCR5","CR2","AICDA")
tumor  <- c("AGR3","VTCN1","SOSTDC1","CACNG4")

xenium <- AddModuleScore(xenium, list(immune), name = "Immune")
xenium <- AddModuleScore(xenium, list(tumor), name = "Tumor")

ic <- ImageFeaturePlot(
  xenium,
  features = c("Immune1","Tumor1"),
  blend = TRUE,
  cols = c("grey90","blue","red"),
  blend.threshold = 0.5,
  size = 1,
  dark.background = FALSE
)

ic[[3]]+ coord_cartesian(
  xlim = c(4500, 5500),
  ylim = c(6000, 7000)
) + ic[[4]] + plot_layout(widths = c(4, 1))

ic[[3]]+ ic[[4]] + plot_layout(widths = c(4, 1))




#highlight according to cell types to be compared with score plot
highlight_types <- c(
  "B cells Memory",
  "B cells Naive",
  "T cells CD4+",
  "T cells CD8+",
  "Macrophage",
  "DCs",
  "Monocyte",
  "Cancer LumB SC",
  "Cancer LumA SC",
  "Cancer Basal SC"
)

cells_use <- as.character(xenium$first_type) %in% highlight_types

xenium$highlight_type <- "Background"

xenium$highlight_type[cells_use] <- as.character(xenium$first_type)[cells_use]
highlight_pals <- pals_use
highlight_pals["Background"] <- "grey90"
highlight_pals["B cells Memory"] <- "blue"


ic2 <-ImageDimPlot(
  xenium,
  group.by = "highlight_type",
  cols = highlight_pals,
  size = 1.2,
  dark.background = FALSE
)+ coord_cartesian(
  xlim = c(4500, 5500),
  ylim = c(6000, 7000)
)

ic2 + ic[[3]]+ coord_cartesian(
  xlim = c(4500, 5500),
  ylim = c(6000, 7000)
) 


#iterated for fibroblasts
highlight_types <- c(
  "CAFs myCAF-like",
  "CAFs MSC iCAF-like",
  "T cells CD8+",
  "Endothelial Lymphatic LYVE1",
  "Cancer LumB SC",
  "Cancer LumA SC",
  "Cancer Basal SC",
  "Endothelial RGS5"
)

cells_use <- as.character(xenium$first_type) %in% highlight_types

xenium$highlight_type <- "Background"

xenium$highlight_type[cells_use] <- as.character(xenium$first_type)[cells_use]
highlight_pals <- pals_use
highlight_pals["Background"] <- "grey90"
highlight_pals["Cancer LumB SC"] <- "blue"
highlight_pals["Cancer LumA SC"] <- "lightgreen"
highlight_pals["CAFs MSC iCAF-like"] <- "yellow"
highlight_pals["CAFs myCAF-like"] <- "orange"
highlight_pals["T cells CD8+"] <- "red"


fc <-ImageDimPlot(
  xenium,
  group.by = "highlight_type",
  cols = highlight_pals,
  size = 1.5,
  dark.background = FALSE
)+ coord_cartesian(
  xlim = c(2500, 5000),
  ylim = c(3000, 6000)
)
fc

fc + p[[3]]+ coord_cartesian(
  xlim = c(4500, 5500),
  ylim = c(3000, 5000)
)  + p[[4]] 




#spatial info awared analyses
#lambda: weight of the spatial information matrix on the gene expression matrix
#lambda: 0-1 larger more spatial
#lambda smaller would be more about cell types, larger would be more of spatial niches
#features=variable: only use highly variable genes
#k_geom: k nearest neighbors
#k_geom larger- niche would also be larger
xenium <- subset(
  xenium1,
  nFeature_Xenium >= 100
)
banksy_x <- RunBanksy(xenium,
                      lambda = 0.55, verbose = TRUE,
                    assay = "Xenium", slot = "counts",  features = "variable",
                    k_geom = 50
)

banksy_x <- RunPCA(banksy_x, assay = "BANKSY", reduction.name = "pca.banksy", npcs = 30, features = rownames(banksy_x))
banksy_x <- RunUMAP(banksy_x, dims = 1:30, reduction = "pca.banksy")
# elbow point is when the variance explained sudden drops, if following PC values show relatively 
# low values, then PC is probably saturated and could be considered abolished.
ElbowPlot(banksy_x, ndims = 20, reduction = 'pca')

#find clusters in banksy
banksy_x <- FindNeighbors(banksy_x, dims = 1:20, assay = "BANKSY", reduction = "pca.banksy")
banksy_x <- FindClusters(banksy_x, resolution = 0.3, algorithm = 4, cluster.name = "banksy_cluster")
#visualization
DimPlot(banksy_x, reduction = "umap", group.by = "banksy_cluster", label = TRUE) +
  labs(title = "UMAP of xenium banksy clusters")
ImageDimPlot(banksy_x, group.by = "banksy_cluster", cols = kelly) +
  labs(title = "Spatial clusters identified by Banksy on Xenium image")


xcell_cluster<-table(
  banksy_x$banksy_cluster,
  xenium$first_type
)
View(as.data.frame.matrix(xcell_cluster))

#heatmap for cell type composition of each cluster
mat <- as.matrix(xcell_cluster)
log_mat <- log2(mat + 1)
pheatmap(
  log_mat,
  cluster_rows = TRUE,     
  cluster_cols = FALSE,    
  scale = "none"
)




#find marker genes
#min.pct could be smaller for interface markers/niche-specific genes/gradient-like expression
#only.pos could be FALSE, to see immune-suppressed marker genes
banksy_markers <- FindAllMarkers(banksy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = "BANKSY", pvalue.cutoff = 0.05)
banksy_top_markers <- banksy_markers %>%
  group_by(cluster) %>%
  top_n(1, avg_log2FC) %>%
  ungroup()
FeaturePlot(banksy, features = banksy_top_markers$gene, ncol = 5) +
  plot_annotation(title = "Top marker genes for each Banksy cluster")
ImageFeaturePlot(banksy, features = banksy_top_markers$gene) +
  plot_annotation(title = "Top marker genes for each Banksy cluster on image data")


#crop and zoom in to see more details
ImageDimPlot(xenium, group.by = "first_type",cols = celltype_colours, 
             dark.background = FALSE,size = 1.5)+ coord_cartesian(
  xlim = c(2800, 5200),
  ylim = c(2000, 5000)
)
ImageDimPlot(xenium, group.by = "first_type",cols = celltype_colours, 
             dark.background = FALSE,size = 1.5)+ coord_cartesian(
               xlim = c(4800, 7200),
               ylim = c(100, 2000)
             )
celltype_colours <- c("B cells Memory" = "grey90", "B cells Naive" = "grey90", "Plasmablasts" = "grey90",
                                 "T cells CD4+" = "grey90", "T cells CD8+" = "magenta", "Cycling T-cells" = "grey90",
                                 "NK cells" = "grey90", "NKT cells" = "grey90",
                                 "Macrophage" = "magenta", "Monocyte" = "grey90", "Cycling_Myeloid" = "grey90", "DCs" = "grey90",
                      "Cancer Basal SC" = "lightblue", "Cancer Cycling", "lightblue", "Cancer Her2 SC" = "lightblue",
                      "Cancer LumA SC" = "lightblue", "Cancer LumB SC" = "lightblue",
                      "Luminal Progenitors" = "grey90", "Mature Luminal" = "grey90", "Myoepithelial" = "grey90",
                      "CAFs MSC iCAF-like" = "gold", "CAFs myCAF-like" = "black",
                                 "Endothelial ACKR1" = "black", "Endothelial CXCL12" = "black", "Endothelial Lymphatic LYVE1" = "black",
                                 "Endothelial RGS5" = "black",
                                 "Cycling PVL" = "grey90", "PVL Differentiated" = "grey90", "PVL Immature" = "grey90")               
            
