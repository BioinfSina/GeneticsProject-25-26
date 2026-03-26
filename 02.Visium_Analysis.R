#load necessary libraries
library(devtools)
library(gprofiler2)
library(remotes)
library(dplyr)
library(ggplot2)
library(Rfast2)
library(Seurat)
library(patchwork)
library(stringr)
library(tidyverse)
library(semla)
library(BiocNeighbors)
library(ComplexHeatmap)
library(SeuratWrappers)
library(Banksy)
library(gprofiler2)
library(org.Hs.eg.db)


# Load the spatial transcriptomics data
visium <- Load10X_Spatial("/mnt/beegfs7/home2/reid/rx241/practical/Ultima_Visium", bin.size = "polygons")

# Explore the Seurat object
str(visium)
View(visium)
ImageDimPlot(visium)

# Visualize the number of transcripts per spot
SpatialFeaturePlot(visium, features = "nCount_Spatial")


##Preprocessing
# Calculate QC metrics
visium[["percent.mt"]] <- PercentageFeatureSet(visium, pattern = "^MT-")
visium[["nCount_Spatial"]] <- Matrix::colSums(visium@assays$Spatial) 


# Visualize QC metrics (low-quality spots and high mitochondrial content)
VlnPlot(visium, features = c("nCount_Spatial", "percent.mt"), ncol = 2)
# Remove low-quality cells and check the difference in metrics
#dependent on sample preparation method/sequencing depth, counts could be different
# could use percentage as QC, rather than absolute counts?
low <- quantile(visium$nCount_Spatial, 0.05)
high <- quantile(visium$nCount_Spatial, 0.99)
#percent.mt could be higher in some types of breast cancer cells-such as epithelial cells where OXPHOS acitivity is upregulated
# also in some hypoxia area, percent.mt could be higher due to metabolic control, therefore thershold for percent.mt should be modified
mt_high <- quantile(na.omit(visium$percent.mt), 0.99)
visium <- subset(
  visium,
  subset = nCount_Spatial > low &
    nCount_Spatial < high &
    percent.mt < mt_high
)

VlnPlot(
  visium,
  features = c("nCount_Spatial", "percent.mt"),
  ncol = 2
)


visium <- NormalizeData(visium)
#apply the normalized factors back to data
visium <- ScaleData(visium)


###Dimensionality Reduction
# Perform PCA on the preprocessed data
visium <- FindVariableFeatures(visium)
visium <- RunPCA(visium, assay = "Spatial.Polygons", npcs = 30, verbose = FALSE)
#Check elbow plot to determine the number of significant PCs
ElbowPlot(visium, ndims = 30, reduction = 'pca')

# Perform UMAP on the PCA results
visium <- RunUMAP(visium, reduction = "pca", dims = 1:30)
#Visualize the UMAP results
umap <- DimPlot(visium, reduction = "umap", label = TRUE,label.size=3) + ggtitle("UMAP of Spatial Transcriptomics Data") +
  theme(
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
umap






###Clustering
#Perform default clustering using Leiden algorithm on the PCA results
visium <- FindNeighbors(visium, reduction = "pca", dims = 1:30)
#as signals are mixed per spot for visium, resolution could also be slightly smaller to avoid subclustering
visium <- FindClusters(visium, resolution = 0.55, algorithm = 4, cluster.name="clusters_0.5") 

#UMAP plot of clusters
visium <- RunUMAP(
  visium,
  dims = 1:30,
  reduction = "pca"
)
umap_clust_0.5 <- DimPlot(visium, reduction = "umap", group.by = "clusters_0.5", cols = kelly,label = TRUE)
#Spatial plot of clusters (we will reduce the saturation of the background image to see better)
spatial_clust_0.5 <- SpatialDimPlot(visium, group.by = "clusters_0.5", label = TRUE, image.alpha = 0.5)

umap_clust_0.5 +spatial_clust_0.5



#visualize cell type composition of each cluster
tab <- table(visium$clusters_0.5, visium$first_type)

cell_cluster_pct <- prop.table(tab, margin = 1)
dfv <- as.data.frame(cell_cluster_pct)
dfv$fill_group <- ifelse(dfv$Freq >= 0.1, as.character(dfv$Var2), "Other")
pals_use <- setNames(kelly, celltypes)
pals_use <- c(pals_use, Other = "grey80")

vv_gg <-ggplot(dfv, aes(x = Var1, y = Freq, fill = fill_group)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = pals_use) +
  labs(x = "Visium cluster", y = "Cell type proportion") +
  theme_classic()







###Spatially Variable Features

# Identify spatially variable features using the "moransi" method
#nfeatures could be more-like 300-500, 
visium <- FindSpatiallyVariableFeatures(visium, selection.method = "moransi", nfeatures = 300)
# View the top spatially variable features
head(SpatiallyVariableFeatures(visium), 10)

# Visualize the top three spatially variable features on a spatial plot
top3_spatial_features <- head(SpatiallyVariableFeatures(visium), 3)
SpatialFeaturePlot(visium, features = top3_spatial_features, ncol = 3) 




###Spatial clustering
# Run BANKSY on the Seurat object
visium <- FindVariableFeatures(visium)

#force to only valid features
hvgs <- VariableFeatures(visium)
visium <- subset(visium, features = hvgs) 

#lambda: weight of the spatial information matrix on the gene expression matrix
#lambda: 0-1 larger more spatial
#lambda smaller would be more about cell types, larger would be more of spatial niches
#features=variable: only use highly variable genes
#k_geom: k nearest neighbors
#k_geom larger- niche would also be larger, for visium this could be smaller
#if we modify lambda to be larger, k_geom should better be also larger.
banksy <- RunBanksy(visium, 
                    assay = "Spatial.Polygons", 
                    lambda = 0.5, 
                    features = hvgs, 
                    k_geom = 50)

#dimensionality reduction and visualization
banksy <- RunPCA(banksy, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(banksy), npcs = 30)
banksy <- FindNeighbors(banksy, reduction = "pca.banksy", dims = 1:30)
banksy <- FindClusters(banksy, resolution = 0.5, algorithm = 4, cluster.name = "banksy_0.5")
banksy <- RunUMAP(banksy, reduction = "pca.banksy", dims = 1:30, reduction.name = "umap.banksy")
DimPlot(banksy, reduction = "umap.banksy", group.by = "banksy_0.5", label = TRUE) + ggtitle("BANKSY Clusters on UMAP")


#visualization of banksy clusters
kelly <- setNames(
  kelly[1:length(levels(banksy$banksy_0.5))],
  levels(banksy$banksy_0.5)
)
polychrome <- pals::polychrome(36)
SpatialDimPlot(
  banksy,
  group.by = "banksy_0.5",
  label = TRUE,
  cols = kelly,
  image.alpha = 0.5
)

#set color to clusters
banksy$banksy_0.5 <- factor(banksy$banksy_0.5)
levels(banksy$banksy_0.5)
cluster_levels <- levels(banksy$banksy_0.5)

banksy_cols <- setNames(
  kelly[1:length(cluster_levels)],
  cluster_levels
)




#heatmap
cell_cluster<-table(
  banksy$banksy_0.5,
  visium$first_type
)
View(as.data.frame.matrix(cell_cluster))


mat <- prop.table(cell_cluster, 1)

log_count <- log2(cell_cluster + 1)

pheatmap(
  log_count,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Cell Counts of BANKSY Clusters(log2 transformed)"
)

prop_cluster <- prop.table(cell_cluster, 1)
df <- as.data.frame(prop_cluster)
df$fill_group <- ifelse(df$Freq >= 0.05, as.character(df$Var2), "Other")
pals_use <- setNames(kelly, celltypes)
pals_use <- c(pals_use, Other = "grey80")

bv_gg <-ggplot(df, aes(x = Var1, y = Freq, fill = fill_group)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = pals_use) +
  labs(x = "BANKSY cluster", y = "Cell type proportion") +
  theme_classic()
#to compare cell type composition in normal clustering and BANKSY clustering
vv_gg / bv_gg


#calculate main cancer cell type to determine sample BC subtype
cancer_cells <- visium$first_type[
  grepl("Cancer", visium$first_type)
]
counts <- table(cancer_cells)
counts
chisq.test(counts)






cluster_cols <- setNames(kelly[1:length(visium$first_type)], visium$first_type)
#cluster cell identities for each banksy subset(iterated for all subsets)
meta <- banksy_13456@meta.data
meta <- meta[!is.na(meta$first_type), ]
tab <- table(meta$banksy_0.5, meta$first_type)
cell_cluster_pct <- prop.table(tab, margin = 1)

dfv <- as.data.frame(cell_cluster_pct)
dfv <- dfv[dfv$Var1 %in% c("3", "1","4","6", "5") & dfv$Freq > 0, ]
dfv$Var1 <- factor(as.character(dfv$Var1), levels = c("3",  "1","4","6","5"))
dfv$fill_group <- ifelse(dfv$Freq >= 0.06, as.character(dfv$Var2), "Other")

dfv <- dfv[!is.na(dfv$fill_group), ]

ggl <- ggplot(dfv, aes(x = Var1, y = Freq, fill = fill_group)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = pals_use, na.translate = FALSE) +
  scale_x_discrete(drop = TRUE) + 
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "BANKSY cluster", y = "Cell type proportion", fill = "Cell Type") +
  theme_classic()

#subset to see more details in local neighborhood
banksy_5710 <- subset(banksy, idents = c("5","7","10"))
DimPlot(banksy_5710, reduction = "umap.banksy", group.by = "banksy_0.5", cols = banksy_cols,label = TRUE) + ggtitle("BANKSY Cluster-5,7&10 on UMAP")
SpatialDimPlot(banksy_5710, group.by = "banksy_0.5", label = TRUE, cols = banksy_cols,image.alpha = 0.3) + ggtitle("BANKSY Cluster-5,7&10 on Spatial Plot")
bsp5710



###Finding Marker Genes
# mixed signals--logfc could be diluted therefore logfc.threshold could not be too harsh.
# spatially-restricted expression pattern could be low in min.pct
markers5710 <- FindAllMarkers(banksy_5710, only.pos = TRUE, logfc.threshold = 0.1, pvalue.cutoff = 0.05)
# View top markers for each cluster
# avg_log2FC could be high due to gene expressed in very restricted area
# could add pct.1>0.1 and p values to filter noise genes
top10 <- markers5710 %>% 
  filter(pct.1 > 0.1) %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 50) %>% 
  ungroup()
top10$gene <- gsub("\\.m[0-9]+$", "", top10$gene)

markers_10 <- markers5710 %>%
  filter(cluster == 5)

#GO enrichment analysis
#gprofiler2 package, gost function
genes5 <- markers_10$gene
genes5 <- gsub("\\.m0$", "", genes5)

gost_res <- gost(
  query = genes5,
  organism = "hsapiens",
  sources = "GO:BP"
)
top_go5 <- gost_res$result %>%
  filter(source == "GO:BP") %>%
  arrange(p_value) %>%
  head(9)
c5<-ggplot(top_go5,
       aes(x = -log10(p_value),
           y = reorder(term_name, p_value))) +
  geom_point(size = 3) +
  labs(
    x = "-log10(p-value)",
    y = "GO Biological Process",
    title = "GO enrichment of cluster 5 markers"
  ) +
  theme_classic()

c5 / c7/ c10


#filter marker genes for Xenium tests(top 15 with large fold change)
View(as.data.frame(markers_7))
markers_10 <- markers_10 %>%
  filter((pct.1 - pct.2) > 0.2) %>%       
  arrange(desc(avg_log2FC)) %>%          
  head(15)
mark10 <- markers_10[, 7]
print(mark10)

#Visualise the 4 top markers for cluster 1
fp <- SpatialFeaturePlot(visium, features = gost_res$result, ncol = 5)
# Visualise the spatial location of cluster 1
sp <- SpatialDimPlot(visium, cells.highlight = WhichCells(visium, idents = "5")) + ggtitle("Cluster 1 Spatial Location")
sp
fp


markers_all <- FindAllMarkers(visium_0.5, logfc.threshold = 0.1, pvalue.cutoff = 0.05)
# View top markers for each cluster
# avg_log2FC could be high due to gene expressed in very restricted area
# could add pct.1>0.1 and p values to filter noise genes
top10 <- markers_all %>% filter(p_val_adj < 0.05, pct.1 > 0.1) %>% group_by(cluster_0.3) %>% top_n(n = 10, wt = avg_log2FC)






#Cellchat analysis
# Expression matrix
input <- visium[["Spatial.Polygons"]]$counts

# Set identities (cell type labels)
Idents(visium) <- "first_type"
labels <- Idents(visium)  
meta <- data.frame(labels = labels, row.names = names(labels)) 

spatial_locs <- GetTissueCoordinates(visium, scale = NULL)[ , c("x", "y")] 
# For segmented Visium: use cell centroid coordinates from metadata
cellChat <- createCellChat(object = visium, group.by = "first_type", assay = "Spatial.Polygons",coordinates = spatial_locs)

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
cellChat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellChat)
table(Idents(visium))


# Identify overexpressed genes and interactions
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#The following function computes the communication probability at the signaling pathway level
#If it takes a long time to compute, so you can skip it and load precomputed results instead
cellchat <- computeCommunProb(cellchat, distance.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
#only consider interactions that are supported by at least 10 cells
cellchat <- filterCommunication(cellchat, min.cells = 10) 
cellchat <- aggregateNet(cellchat)




#visualization
#1.circle plot
#find all significant signalling pathways
pathways.show.all <- cellchat@netP$pathways
pathways.selected <- c("COLLAGEN", "LAMININ")
laminin <- c("LAMININ")

# Circle plos of the overall communication network
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T  , label.edge= F, title.name = "Number of interactions", top = 0.3)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T  , label.edge= F, title.name = "Interaction weights/strength", top = 0.3)

#circle plots of specific cell types outgoing/ingoing communications
mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))

celltype <- "Macrophage"   
celltype <- "Cycling T-cells"
celltype <- "T cells CD8+"
celltype <- "T cells CD4+"
celltype <- "CAFs MSC iCAF-like"
mat_out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat),
                  dimnames = dimnames(mat))
mat_out[celltype, ] <- mat[celltype, ]

netVisual_circle(mat_out,  vertex.weight = groupSize,
                 weight.scale = TRUE,    
                 label.edge = FALSE,
                 vertex.label.cex = 0.7,
                 edge.width.max = 10,title.name = paste(celltype, "→ others"))
mat_in <- matrix(0, nrow = nrow(mat), ncol = ncol(mat),
                 dimnames = dimnames(mat))
mat_in[, celltype] <- mat[, celltype]

netVisual_circle(mat_in, ,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,    
                 label.edge = FALSE,
                 vertex.label.cex = 0.7,
                 edge.width.max = 10,title.name = paste("others →", celltype))

netVisual_circle(
  mat2,
  vertex.weight = groupSize,
  weight.scale = TRUE,    
  label.edge = FALSE,
  vertex.label.cex = 1.2,
  edge.width.max = 10,
  title.name = paste(celltype, "↔ others"),
  vertex.size = 25
)



#iterated for several cell types of interest
target_cells <- c("Macrophage","Cycling T-cells", "T cell CD4+", "T cell CD8+")
par(mfrow = c(1, 2), mar = c(1, 1, 1, 1)) 
for (i in 1:length(target_cells)) {
  current_cell <- target_cells[i]
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  
  if (current_cell %in% rownames(mat)) {
    mat2[current_cell, ] <- mat[current_cell, ]
    
    netVisual_circle(
      mat2,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      label.edge = FALSE,
      edge.weight.max = max(mat),
      title.name = paste(current_cell, "→ others"),
      vertex.label.cex = 0.8,     
      vertex.size = 20
    )
  } else {
    print(paste("Warning: ", current_cell, " not found in mat rownames."))
  }
}

）
par(mfrow = c(1, 1))


#2.bubble plot
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(
  cellchat, 
  sources.use = colnames(cellchat@net$weight),
  targets.use = "T cells CD8+",
  remove.isolate = TRUE
)
netVisual_bubble(
  cellchat, 
  sources.use = c("CAFs MSC iCAF-like"), 
  targets.use = c(29:30),
  remove.isolate = TRUE,
  thresh = 0.01 
)


#3.spatial plot of cellchat results
netVisual_aggregate(cellchat, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
netVisual_aggregate(cellchat, signaling = laminin, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)



