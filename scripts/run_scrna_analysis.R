#!/usr/bin/env Rscript

################################################################################
# scRNA-seq Analysis Workflow
#
# Project:
#   Public scRNA-seq workflow demonstration using two 10x Genomics-style samples
#   from GSE262288.
#
# Data source:
#   Luo L., Yang P., Mastoraki S., et al.
#   Single-cell RNA sequencing identifies molecular biomarkers predicting late
#   progression to CDK4/6 inhibition in patients with HR+/HER2- metastatic
#   breast cancer.
#   Molecular Cancer, 2025.
#   PMID: 39955556 | DOI: 10.1186/s12943-025-02226-9
#
# Samples used:
#   GSM8162620 — Patient 3 ascites scRNA-seq sample 1
#   GSM8162621 — Patient 3 ascites scRNA-seq sample 2
#
# Purpose:
#   This script demonstrates a reproducible Seurat-based scRNA-seq workflow:
#     1. Load 10x matrices
#     2. Create Seurat objects
#     3. Perform per-sample QC
#     4. Merge samples
#     5. Normalize, scale, and run PCA
#     6. Compare uncorrected PCA UMAP, Harmony UMAP, and CCA-integrated UMAP
#     7. Cluster cells using Harmony embeddings
#     8. Identify cluster marker genes
#     9. Annotate cell types using SingleR
#    10. Save figures, tables, Seurat objects, and session info
################################################################################


# ==============================================================================
# 1. Parameters
# ==============================================================================

set.seed(1234)

sample_dirs <- c(
  sam1 = "sample/sample_1",
  sam2 = "sample/sample_2"
)

outdir <- "results"

min_cells <- 3
min_features_create <- 200

qc_min_features <- 200
qc_max_features <- 6500
qc_max_percent_mt <- 15

n_variable_features <- 2000
n_pcs <- 30
elbow_ndims <- 50

cluster_resolution <- 0.1

selected_marker_genes <- c(
  "S100P", "SLC11A1", "C1R", "KLRB1",
  "TRAC", "NUSAP1", "HLA-DQB2", "TCL1A"
)


# ==============================================================================
# 2. Load libraries
# ==============================================================================

required_packages <- c(
  "Seurat",
  "dplyr",
  "patchwork",
  "harmony",
  "ggplot2",
  "celldex",
  "SingleR"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are missing: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(harmony)
  library(ggplot2)
  library(celldex)
  library(SingleR)
})


# ==============================================================================
# 3. Create output directories
# ==============================================================================

figure_dir <- file.path(outdir, "figures")
qc_dir <- file.path(figure_dir, "qc")
dimred_dir <- file.path(figure_dir, "dimensionality_reduction")
integration_dir <- file.path(figure_dir, "integration")
marker_dir <- file.path(figure_dir, "markers")
annotation_dir <- file.path(figure_dir, "annotation")

table_dir <- file.path(outdir, "tables")
marker_table_dir <- file.path(table_dir, "markers")
annotation_table_dir <- file.path(table_dir, "annotation")

object_dir <- file.path(outdir, "objects")
log_dir <- file.path(outdir, "logs")

dirs_to_create <- c(
  figure_dir,
  qc_dir,
  dimred_dir,
  integration_dir,
  marker_dir,
  annotation_dir,
  table_dir,
  marker_table_dir,
  annotation_table_dir,
  object_dir,
  log_dir
)

invisible(lapply(dirs_to_create, dir.create, recursive = TRUE, showWarnings = FALSE))


# ==============================================================================
# 4. Helper functions
# ==============================================================================

read_10x_gene_expression <- function(data_dir) {
  if (!dir.exists(data_dir)) {
    stop("Input directory does not exist: ", data_dir)
  }

  counts <- Read10X(data.dir = data_dir)

  # Some 10x outputs are returned as a list, especially multiome/multimodal data.
  # For this workflow, we use the Gene Expression matrix.
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
      warning(
        "Read10X returned a list, but no 'Gene Expression' element was found. ",
        "Using the first matrix."
      )
    }
  }

  return(counts)
}


save_ggplot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}


save_plot_expression <- function(filename, expr, width = 10, height = 8, res = 300) {
  png(
    filename = filename,
    width = width,
    height = height,
    units = "in",
    res = res
  )
  on.exit(dev.off(), add = TRUE)
  force(expr)
}


get_qc_summary <- function(seurat_obj, sample_name, stage) {
  metadata <- seurat_obj@meta.data

  data.frame(
    sample = sample_name,
    stage = stage,
    n_cells = ncol(seurat_obj),
    median_nFeature_RNA = median(metadata$nFeature_RNA),
    median_nCount_RNA = median(metadata$nCount_RNA),
    median_percent_mt = median(metadata$percent.mt),
    stringsAsFactors = FALSE
  )
}


join_layers_if_available <- function(seurat_obj, assay = "RNA") {
  # Seurat v5 uses layers. JoinLayers is needed before some downstream operations
  # after merging multi-sample objects.
  if (exists("JoinLayers", mode = "function")) {
    seurat_obj <- JoinLayers(seurat_obj, assay = assay)
  }

  return(seurat_obj)
}


get_assay_matrix <- function(seurat_obj, assay = "RNA", layer_or_slot = "data") {
  # Compatible with Seurat v5 and older Seurat versions.
  if (exists("LayerData", mode = "function")) {
    mat <- LayerData(seurat_obj, assay = assay, layer = layer_or_slot)
  } else {
    mat <- GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
  }

  return(mat)
}


# ==============================================================================
# 5. Load samples and create Seurat objects
# ==============================================================================

message("Loading 10x matrices and creating Seurat objects...")

sample_objects <- list()
qc_summaries <- list()

for (sample_name in names(sample_dirs)) {
  message("Processing sample: ", sample_name)

  counts <- read_10x_gene_expression(sample_dirs[[sample_name]])

  obj <- CreateSeuratObject(
    counts = counts,
    min.cells = min_cells,
    min.features = min_features_create,
    project = sample_name
  )

  obj <- RenameCells(obj, add.cell.id = sample_name)

  obj$sample_id <- sample_name
  obj$orig.ident <- sample_name

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  qc_summaries[[paste0(sample_name, "_before_qc")]] <- get_qc_summary(
    obj,
    sample_name,
    "before_qc"
  )

  p_before_qc <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  ) +
    ggtitle(paste0(sample_name, " before QC"))

  save_ggplot(
    p_before_qc,
    file.path(qc_dir, paste0(sample_name, "_before_QC_vlnplot.png")),
    width = 12,
    height = 5
  )

  obj <- subset(
    obj,
    subset = nFeature_RNA > qc_min_features &
      nFeature_RNA < qc_max_features &
      percent.mt < qc_max_percent_mt
  )

  qc_summaries[[paste0(sample_name, "_after_qc")]] <- get_qc_summary(
    obj,
    sample_name,
    "after_qc"
  )

  p_after_qc <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  ) +
    ggtitle(paste0(sample_name, " after QC"))

  save_ggplot(
    p_after_qc,
    file.path(qc_dir, paste0(sample_name, "_after_QC_vlnplot.png")),
    width = 12,
    height = 5
  )

  sample_objects[[sample_name]] <- obj
}

qc_summary_table <- bind_rows(qc_summaries)

write.csv(
  qc_summary_table,
  file = file.path(table_dir, "qc_summary_by_sample.csv"),
  row.names = FALSE
)


# ==============================================================================
# 6. Merge samples
# ==============================================================================

message("Merging samples...")

if (length(sample_objects) == 1) {
  merged_sam <- sample_objects[[1]]
} else {
  merged_sam <- merge(
    x = sample_objects[[1]],
    y = sample_objects[-1],
    project = "Merged_sam"
  )
}

DefaultAssay(merged_sam) <- "RNA"

saveRDS(
  merged_sam,
  file = file.path(object_dir, "merged_raw_qc_filtered_seurat_object.rds")
)


# ==============================================================================
# 7. Preprocessing: normalization, variable features, scaling, PCA
# ==============================================================================

message("Running normalization, variable feature selection, scaling, and PCA...")

merged_sam <- NormalizeData(merged_sam)
merged_sam <- FindVariableFeatures(
  merged_sam,
  selection.method = "vst",
  nfeatures = n_variable_features
)
merged_sam <- ScaleData(merged_sam)
merged_sam <- RunPCA(
  merged_sam,
  npcs = max(n_pcs, elbow_ndims),
  verbose = FALSE
)

p_pca_12 <- DimPlot(
  merged_sam,
  reduction = "pca",
  group.by = "orig.ident"
) +
  ggtitle("PCA: PC1 vs PC2 by sample")

p_pca_15 <- DimPlot(
  merged_sam,
  reduction = "pca",
  dims = c(1, 5),
  group.by = "orig.ident"
) +
  ggtitle("PCA: PC1 vs PC5 by sample")

p_pca_loadings <- VizDimLoadings(
  merged_sam,
  dims = 1:2,
  reduction = "pca"
)

p_elbow <- ElbowPlot(
  merged_sam,
  reduction = "pca",
  ndims = elbow_ndims
)

save_ggplot(
  p_pca_12,
  file.path(dimred_dir, "merged_pca_PC1_PC2_by_sample.png"),
  width = 7,
  height = 5
)

save_ggplot(
  p_pca_15,
  file.path(dimred_dir, "merged_pca_PC1_PC5_by_sample.png"),
  width = 7,
  height = 5
)

save_ggplot(
  p_pca_loadings,
  file.path(dimred_dir, "merged_pca_loadings_PC1_PC2.png"),
  width = 10,
  height = 6
)

save_ggplot(
  p_elbow,
  file.path(dimred_dir, "merged_pca_elbowplot.png"),
  width = 7,
  height = 5
)


# ==============================================================================
# 8. UMAP before integration
# ==============================================================================

message("Running uncorrected PCA UMAP...")

merged_sam <- RunUMAP(
  merged_sam,
  reduction = "pca",
  reduction.name = "umap.pca",
  dims = 1:n_pcs,
  plot_convergence = FALSE
)

p_umap_pca_sample <- DimPlot(
  merged_sam,
  reduction = "umap.pca",
  group.by = "orig.ident"
) +
  ggtitle("Before integration: PCA UMAP by sample")

save_ggplot(
  p_umap_pca_sample,
  file.path(integration_dir, "umap_before_integration_by_sample.png"),
  width = 7,
  height = 5
)


# ==============================================================================
# 9. Harmony integration
# ==============================================================================

message("Running Harmony integration...")

merged_sam <- RunHarmony(
  object = merged_sam,
  group.by.vars = "orig.ident",
  dims.use = 1:n_pcs,
  plot_convergence = FALSE
)

merged_sam <- RunUMAP(
  merged_sam,
  reduction = "harmony",
  dims = 1:n_pcs,
  reduction.name = "umap.harmony",
  reduction.key = "UMAP_HARMONY_"
)

p_umap_harmony_sample <- DimPlot(
  merged_sam,
  reduction = "umap.harmony",
  group.by = "orig.ident"
) +
  ggtitle("After Harmony integration: Harmony UMAP by sample")

save_ggplot(
  p_umap_harmony_sample,
  file.path(integration_dir, "umap_after_harmony_by_sample.png"),
  width = 7,
  height = 5
)


# ==============================================================================
# 10. CCA integration for comparison
# ==============================================================================

message("Running CCA integration for comparison...")

cca_sample_objects <- sample_objects

cca_sample_objects <- lapply(
  cca_sample_objects,
  function(obj) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(
      obj,
      selection.method = "vst",
      nfeatures = n_variable_features,
      verbose = FALSE
    )
    return(obj)
  }
)

anchors <- FindIntegrationAnchors(
  object.list = cca_sample_objects,
  dims = 1:n_pcs
)

cca_integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:n_pcs
)

DefaultAssay(cca_integrated) <- "integrated"

cca_integrated <- ScaleData(cca_integrated, verbose = FALSE)
cca_integrated <- RunPCA(cca_integrated, npcs = n_pcs, verbose = FALSE)

cca_integrated <- RunUMAP(
  cca_integrated,
  reduction = "pca",
  dims = 1:n_pcs,
  reduction.name = "umap.cca",
  reduction.key = "UMAP_CCA_"
)

cca_integrated <- FindNeighbors(
  cca_integrated,
  reduction = "pca",
  dims = 1:n_pcs
)

cca_integrated <- FindClusters(
  cca_integrated,
  resolution = cluster_resolution
)

p_umap_cca_sample <- DimPlot(
  cca_integrated,
  reduction = "umap.cca",
  group.by = "orig.ident"
) +
  ggtitle("After CCA integration: CCA UMAP by sample")

p_umap_cca_clusters <- DimPlot(
  cca_integrated,
  reduction = "umap.cca",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("CCA-integrated UMAP by cluster")

save_ggplot(
  p_umap_cca_sample,
  file.path(integration_dir, "umap_after_cca_by_sample.png"),
  width = 7,
  height = 5
)

save_ggplot(
  p_umap_cca_clusters,
  file.path(integration_dir, "umap_after_cca_by_cluster.png"),
  width = 7,
  height = 5
)

p_integration_comparison <- p_umap_pca_sample + p_umap_harmony_sample + p_umap_cca_sample +
  plot_layout(ncol = 3)

save_ggplot(
  p_integration_comparison,
  file.path(integration_dir, "integration_comparison_pca_harmony_cca.png"),
  width = 18,
  height = 5
)

saveRDS(
  cca_integrated,
  file = file.path(object_dir, "cca_integrated_seurat_object.rds")
)


# ==============================================================================
# 11. Clustering using Harmony embeddings
# ==============================================================================

message("Clustering cells using Harmony embeddings...")

merged_sam <- FindNeighbors(
  merged_sam,
  reduction = "harmony",
  dims = 1:n_pcs
)

merged_sam <- FindClusters(
  merged_sam,
  resolution = cluster_resolution
)

Idents(merged_sam) <- "seurat_clusters"

p_umap_harmony_clusters <- DimPlot(
  merged_sam,
  reduction = "umap.harmony",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +
  ggtitle("Harmony-integrated UMAP by cluster")

p_umap_harmony_split_sample <- DimPlot(
  merged_sam,
  reduction = "umap.harmony",
  group.by = "seurat_clusters",
  split.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) +
  ggtitle("Harmony clusters split by sample")

save_ggplot(
  p_umap_harmony_clusters,
  file.path(integration_dir, "umap_harmony_by_cluster.png"),
  width = 7,
  height = 5
)

save_ggplot(
  p_umap_harmony_split_sample,
  file.path(integration_dir, "umap_harmony_clusters_split_by_sample.png"),
  width = 12,
  height = 5
)


# ==============================================================================
# 12. Marker gene detection
# ==============================================================================

message("Finding cluster marker genes...")

DefaultAssay(merged_sam) <- "RNA"
merged_sam <- join_layers_if_available(merged_sam, assay = "RNA")

Idents(merged_sam) <- "seurat_clusters"

all_markers <- FindAllMarkers(
  merged_sam,
  only.pos = FALSE,
  logfc.threshold = 0.25
)

logfc_col <- if ("avg_log2FC" %in% colnames(all_markers)) {
  "avg_log2FC"
} else {
  "avg_logFC"
}

write.csv(
  all_markers,
  file = file.path(marker_table_dir, "all_cluster_markers.csv"),
  row.names = FALSE
)

upregulated_markers <- all_markers %>%
  filter(.data[[logfc_col]] > 0.25) %>%
  arrange(cluster, desc(.data[[logfc_col]]))

downregulated_markers <- all_markers %>%
  filter(.data[[logfc_col]] < -0.25) %>%
  arrange(cluster, .data[[logfc_col]])

write.csv(
  upregulated_markers,
  file = file.path(marker_table_dir, "upregulated_markers_by_cluster.csv"),
  row.names = FALSE
)

write.csv(
  downregulated_markers,
  file = file.path(marker_table_dir, "downregulated_markers_by_cluster.csv"),
  row.names = FALSE
)

top10_upregulated_markers <- upregulated_markers %>%
  group_by(cluster) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top10_upregulated_markers,
  file = file.path(marker_table_dir, "top10_upregulated_markers_by_cluster.csv"),
  row.names = FALSE
)

present_marker_genes <- selected_marker_genes[
  selected_marker_genes %in% rownames(merged_sam)
]

if (length(present_marker_genes) > 0) {
  p_feature_markers <- FeaturePlot(
    merged_sam,
    features = present_marker_genes,
    reduction = "umap.harmony",
    ncol = 4
  )

  save_ggplot(
    p_feature_markers,
    file.path(marker_dir, "featureplot_selected_marker_genes.png"),
    width = 14,
    height = 8
  )
} else {
  warning("None of the selected marker genes were found in the dataset.")
}

top_heatmap_genes <- unique(top10_upregulated_markers$gene)

if (length(top_heatmap_genes) > 0) {
  p_heatmap_top10 <- DoHeatmap(
    merged_sam,
    features = top_heatmap_genes
  ) +
    NoLegend() +
    ggtitle("Top 10 upregulated marker genes per cluster")

  save_ggplot(
    p_heatmap_top10,
    file.path(marker_dir, "heatmap_top10_upregulated_markers_by_cluster.png"),
    width = 12,
    height = 10
  )
}


# ==============================================================================
# 13. Cell type annotation with SingleR
# ==============================================================================

message("Running SingleR annotation...")

DefaultAssay(merged_sam) <- "RNA"

norm_counts <- get_assay_matrix(
  merged_sam,
  assay = "RNA",
  layer_or_slot = "data"
)

ref <- celldex::HumanPrimaryCellAtlasData()

singleR_predictions <- SingleR(
  test = norm_counts,
  ref = ref,
  labels = ref$label.main,
  de.method = "wilcox"
)

singleR_labels <- singleR_predictions$pruned.labels
singleR_labels[is.na(singleR_labels)] <- "Unassigned"
names(singleR_labels) <- rownames(singleR_predictions)

merged_sam <- AddMetaData(
  merged_sam,
  metadata = singleR_labels,
  col.name = "SingleR_HCA"
)

singleR_annotation_table <- data.frame(
  cell_id = rownames(singleR_predictions),
  SingleR_label = singleR_predictions$labels,
  SingleR_pruned_label = singleR_labels,
  stringsAsFactors = FALSE
)

write.csv(
  singleR_annotation_table,
  file = file.path(annotation_table_dir, "singler_cell_annotations.csv"),
  row.names = FALSE
)

singleR_counts <- as.data.frame(table(singleR_labels))
colnames(singleR_counts) <- c("SingleR_pruned_label", "cell_count")

write.csv(
  singleR_counts,
  file = file.path(annotation_table_dir, "singler_celltype_counts.csv"),
  row.names = FALSE
)

p_singler_umap <- DimPlot(
  merged_sam,
  reduction = "umap.harmony",
  group.by = "SingleR_HCA",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) +
  NoLegend() +
  ggtitle("SingleR cell type annotation on Harmony UMAP")

save_ggplot(
  p_singler_umap,
  file.path(annotation_dir, "umap_harmony_singler_annotation.png"),
  width = 9,
  height = 6
)

save_plot_expression(
  file.path(annotation_dir, "singler_score_heatmap.png"),
  {
    print(plotScoreHeatmap(singleR_predictions))
  },
  width = 10,
  height = 8
)

save_plot_expression(
  file.path(annotation_dir, "singler_delta_distribution.png"),
  {
    print(plotDeltaDistribution(
      singleR_predictions,
      ncol = 4,
      dots.on.top = FALSE
    ))
  },
  width = 10,
  height = 6
)


# ==============================================================================
# 14. Save final objects and session info
# ==============================================================================

message("Saving final Seurat object and session info...")

saveRDS(
  merged_sam,
  file = file.path(object_dir, "merged_harmony_singler_seurat_object.rds")
)

sink(file.path(log_dir, "sessionInfo.txt"))
sessionInfo()
sink()

message("Analysis complete.")
message("Figures saved to: ", figure_dir)
message("Tables saved to: ", table_dir)
message("Objects saved to: ", object_dir)
