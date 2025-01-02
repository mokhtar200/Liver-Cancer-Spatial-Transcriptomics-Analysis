# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(SpatialExperiment)

# Step 1: Data Preparation
# Xenium Human Multi-Tissue and Cancer Panel 1: Human Liver Data
data_path <- "https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
image_path <- "https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"

# Load spatial data
spatial_data <- Load10X_Spatial(data.dir = data_path, image.dir = image_path)

# View dataset summary
print(spatial_data)

# Step 2: Quality Control
# Filter cells with low gene count or extreme expression
spatial_data <- subset(spatial_data, subset = nFeature_Spatial > 200 & nCount_Spatial < 25000)

# Step 3: Normalization
# Normalize the dataset
spatial_data <- SCTransform(spatial_data, assay = "Spatial", verbose = FALSE)

# Step 4: Spatially Variable Gene Detection
# Identify spatially variable genes
spatial_features <- FindSpatiallyVariableFeatures(
  spatial_data, assay = "SCT", selection.method = "markvariogram"
)
top_spatial_genes <- head(spatial_features, 5)

# Visualize top spatially variable genes
SpatialFeaturePlot(spatial_data, features = top_spatial_genes)

# Step 5: Clustering and Annotation
# Perform clustering
spatial_data <- FindNeighbors(spatial_data, reduction = "pca")
spatial_data <- FindClusters(spatial_data, resolution = 0.5)

# Visualize clusters
SpatialDimPlot(spatial_data, group.by = "seurat_clusters")

# Annotate clusters
spatial_data <- RenameIdents(spatial_data, `0` = "Tumor Core", `1` = "Immune Infiltrate", `2` = "Stroma")

# Step 6: Gene Expression Analysis
# Plot expression of specific genes of interest
SpatialFeaturePlot(spatial_data, features = c("HNF4A", "AFP"))

# Step 7: Save Outputs
# Save processed dataset
saveRDS(spatial_data, file = "xenium_human_liver_data.rds")

# Save visualizations
ggsave("top_spatial_genes.png", plot = last_plot())
ggsave("spatial_clusters.png", plot = SpatialDimPlot(spatial_data, group.by = "seurat_clusters"))

# Export metadata
write.csv(as.data.frame(spatial_data@meta.data), "metadata.csv")
