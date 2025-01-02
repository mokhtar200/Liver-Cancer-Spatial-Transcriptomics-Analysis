# Liver Cancer Spatial Transcriptomics Analysis

## Project Overview
This project performs an in-depth analysis of liver cancer spatial transcriptomics data using the **Xenium Human Multi-Tissue and Cancer Panel 1: Human Liver Data** dataset from [10x Genomics](https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard).

The analysis identifies spatially variable genes, clusters tissue regions based on gene expression patterns, and visualizes the results to understand the molecular landscape of liver cancer.

## Goals
1. Identify spatially variable genes and their biological relevance.
2. Cluster tissue regions based on spatial transcriptomics profiles.
3. Visualize spatial maps of gene expression and tissue regions.

## Dataset
- **Dataset Name**: Xenium Human Multi-Tissue and Cancer Panel 1: Human Liver Data
- **Source**: [10x Genomics](https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard)

## Analysis Steps
1. **Data Loading**:
   - Load the spatial transcriptomics dataset using the `Seurat` package.
2. **Quality Control**:
   - Filter cells and genes based on predefined thresholds.
3. **Normalization**:
   - Normalize gene expression data using SCTransform.
4. **Spatial Gene Analysis**:
   - Identify spatially variable genes using variogram-based methods.
5. **Clustering and Annotation**:
   - Perform clustering to identify regions of interest and annotate clusters with biological insights.
6. **Visualization**:
   - Generate spatial maps for specific genes and clusters.

## Outputs
- **Top Spatially Variable Genes**: `top_spatial_genes.png`
- **Cluster Visualization**: `spatial_clusters.png`
- **Processed Metadata**: `metadata.csv`

## Requirements
- R version >= 4.0.0
- Libraries: `Seurat`, `ggplot2`, `dplyr`, `SpatialExperiment`

## How to Run
1. Download the dataset:
   [Xenium Human Liver Data](https://www.10xgenomics.com/datasets/human-liver-data-xenium-human-multi-tissue-and-cancer-panel-1-standard)
2. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/liver-cancer-spatial-transcriptomics.git
