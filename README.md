
# 🧬 Single-cell RNA-seq Analysis Pipeline  
**Author:** Harshitha Muddamsetty  
**Tools:** Seurat | Harmony | SingleR | CCA | R  

## 📘 Overview  
This repository contains a comprehensive single-cell RNA-seq (scRNA-seq) analysis pipeline using R and Seurat. The workflow is designed to process 10x Genomics datasets through quality control, batch effect correction (Harmony or CCA), clustering, marker gene identification, and cell type annotation using SingleR.  

The project compares different integration strategies and includes two robust approaches: one using **Harmony** for fast batch correction, and another using **Seurat’s CCA-based IntegrateLayers**.

---

## 🛠️ Key Features

- QC per sample and post-merge filtering  
- Harmony- and CCA-based batch correction workflows  
- Dimensionality reduction (PCA, UMAP, t-SNE)  
- Unsupervised clustering and visualization  
- Marker gene detection  
- Automated cell type annotation using **SingleR**  

---

## 📑 Workflow Summary

### ✨ Main Pipeline (`scRNAseq_pipeline.Rmd`)
- Per-sample QC and filtering  
- Normalization, variable feature selection, scaling  
- PCA → **Harmony** integration on `orig.ident`  
- Clustering, UMAP embedding  
- Marker gene detection (`avg_log2FC > 1`)  
- Cell type annotation with **SingleR** and `HumanPrimaryCellAtlasData()`  
- Output: annotated Seurat object, plots, top marker CSV  

### 🔁 Alternate Batch Effect Correction (`batch_correction_CCA.R`)
- Merges 3 biological replicates  
- Performs CCA-based integration using `IntegrateLayers()`  
- Compares clustering/UMAP plots before vs. after batch correction  
- Highlights integration flexibility across pipelines  

---

## 🧪 Input Requirements  
- 10x Genomics expression matrices (`Read10X`)  
- Pre-organized folder structure by sample or replicate  

---

## 📊 Outputs  
- Filtered and annotated Seurat object  
- Marker gene table: `top_markers.csv`  
- UMAP/tSNE visualizations of clusters and batch groups  
- Feature plots for marker genes  
- SingleR cell type heatmap and annotations  

---

## 📦 Dependencies  
- `Seurat`  
- `harmony`  
- `SingleR`  
- `celldex`  
- `scran`  
- `patchwork`  
- `ggplot2`  
- `dplyr`

---

## 💡 Notes on Integration Strategy  
Multiple trials were conducted to determine the cleanest and most interpretable pipeline:
- QC **before merging** gives cleaner violin plots and better downstream metrics  
- Harmony provides faster batch correction for exploratory workflows  
- CCA (via `IntegrateLayers`) offers higher precision for closely related samples or publication-ready datasets  

---

## 🧠 Citation  
Please cite Seurat, Harmony, and SingleR packages appropriately if this pipeline is used in publications.
