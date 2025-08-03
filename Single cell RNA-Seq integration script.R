# =============================
# 1. Load required packages
# =============================

if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")

library(Seurat)
library(dplyr)
library(patchwork)
library(plotly)

# =============================
# 2. Read the downloaded datasets
# =============================
# Download and extract example dataset from 10x Genomics if needed
untar("C:/Users/User/Downloads/GSM7077865_D1_filtered_feature_bc_matrix.tar.gz", exdir = "pbmc_data_2")
untar("C:/Users/User/Downloads/GSM7077866_G1_filtered_feature_bc_matrix.tar.gz", exdir = "pbmc_data_3")

sample1_dir <- "pbmc_data_2"  # Resting PBMCs (GSM7077865)
sample2_dir <- "pbmc_data_3"  # LPS-stimulated PBMCs (GSM7077866)

sample1_data <- Read10X(data.dir = sample1_dir)
sample2_data <- Read10X(data.dir = sample2_dir)

sample1 <- CreateSeuratObject(counts = sample1_data, project = "Resting", min.cells = 3, min.features = 200)
sample2 <- CreateSeuratObject(counts = sample2_data, project = "LPS", min.cells = 3, min.features = 200)

sample1$condition <- "Resting"
sample2$condition <- "LPS"

pbmc.list <- list(sample1, sample2)

# =============================
# 3. Preprocessing each sample
# =============================

pbmc.list <- lapply(pbmc.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# =============================
# 4. Integration
# =============================

anchors <- FindIntegrationAnchors(object.list = pbmc.list, dims = 1:20)
pbmc.integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

DefaultAssay(pbmc.integrated) <- "integrated"

# =============================
# 5. Downstream analysis
# =============================

pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunPCA(pbmc.integrated, npcs = 30, verbose = FALSE)
pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:10)
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:10)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.5)

# =============================
# 6. Cluster Annotation (Automated + Manual Fallback)
# =============================

top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(top_markers)

DefaultAssay(pbmc.integrated) <- "RNA"
pbmc.integrated <- NormalizeData(pbmc.integrated)
pbmc.integrated <- ScaleData(pbmc.integrated, features = rownames(pbmc.integrated))

marker_db <- list(
  "B Cells"             = c("MS4A1"),
  "CD8 T Cells"         = c("CD8A"),
  "NK Cells"            = c("NCR1"),
  "CD4 T Cells"         = c("CD5", "CD28"),
  "Memory CD4 T Cells"  = c("GATA3", "MAF", "CD40LG"),
  "Naive CD4 T Cells"   = c("LEF1"),
  "Monocytes"           = c("MSR1", "VSIG4"),
  "Dendritic Cells"     = c("FCER1A", "CLEC10A", "CST6"),
  "Platelets"           = c("MYL9", "FMOD"),
  "Regulatory T Cells"  = c("GALM", "RGS1"),
  "Mast Cells"          = c("HDC", "AKAP12", "CLC", "CPA3"),
  "Neutrophils"         = c("LTF", "BPI", "CAMP")
)

Idents(pbmc.integrated) <- "seurat_clusters"
markers <- tryCatch({
  FindAllMarkers(pbmc.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
}, error = function(e) NULL)

if (!is.null(markers) && "cluster" %in% colnames(markers)) {
  top10 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)
  cluster_annotation <- sapply(split(top10, top10$cluster), function(df) {
    scores <- sapply(names(marker_db), function(celltype) {
      sum(df$gene %in% marker_db[[celltype]])
    })
    best <- names(which.max(scores))
    if (all(scores == 0)) return("Unclassified")
    return(best)
  })
  pbmc.integrated <- RenameIdents(pbmc.integrated, cluster_annotation)
  pbmc.integrated$celltype <- Idents(pbmc.integrated)
} else {
  warning("Automatic annotation failed. Consider assigning cluster labels manually.")
}


# =============================
# 7. Visualization
# =============================

DimPlot(pbmc.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
ggsave("UMAP_Celltype_Annotated.png", width = 7, height = 5) #Option 1

DimPlot(pbmc.integrated, reduction = "umap", group.by = "celltype", label = FALSE, repel = TRUE)
ggsave("UMAP_Celltype_Annotated.png", width = 9, height = 5) #Option 2

DimPlot(pbmc.integrated, reduction = "umap", group.by = "condition", label = FALSE, repel = TRUE)
ggsave("UMAP_Celltype_by_Condition.png", width = 9, height = 5)

DimPlot(pbmc.integrated, reduction = "umap", split.by = "condition", group.by = "celltype", label = FALSE, repel = TRUE)
ggsave("UMAP_Celltype_by_Condition_clusters.png", width = 9, height = 5)

DefaultAssay(pbmc.integrated) <- "RNA"

inflammatory_genes <- c("CD14", "IL1B", "TNF", "CXCL10", "IL6", "IFIT1", "ISG15")

for (gene in inflammatory_genes) {
  p <- FeaturePlot(pbmc.integrated, features = gene, split.by = "condition",label = TRUE)
  ggsave(filename = paste0("FeaturePlot_", gene, ".png"), plot = p, width = 7, height = 5)
}

for (gene in inflammatory_genes) {
  p <- VlnPlot(pbmc.integrated, features = gene, split.by = "condition", group.by = "celltype")
  ggsave(filename = paste0("VlnPlot_", gene, ".png"), plot = p, width = 12, height = 5)
}

# =============================
# 8. Differential expression analysis
# =============================

pbmc.integrated <- JoinLayers(pbmc.integrated)  # Required in Seurat v5+
Idents(pbmc.integrated) <- "condition"

#Example with cluster 2
Idents(pbmc.integrated) <- "seurat_clusters"
cells_cluster2 <- WhichCells(pbmc.integrated, idents = 2)
Idents(pbmc.integrated) <- "condition"
markers_cluster2 <- FindMarkers(pbmc.integrated, ident.1 = "LPS", ident.2 = "Resting", cells = cells_cluster2)
head(markers_cluster2)

# =============================
# 9. Save final integrated object
# =============================

saveRDS(pbmc.integrated, file = "pbmc_resting_vs_lps_integrated.rds")

# =============================
# 10. Optional: Interactive UMAP using plotly
# =============================

plotly_umap <- DimPlot(pbmc.integrated, reduction = "umap", group.by = "celltype")
ggplotly(plotly_umap)


