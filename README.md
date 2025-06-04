# Bulk RNA-seq Analysis Shiny App

This repository contains an interactive [Shiny](https://shiny.posit.co/) application for exploring bulk RNA‑sequencing data.  The app provides several analysis modules, including data filtering, quality control, PCA/correlation analyses, differential expression comparisons, clustering, and single‑gene visualization.

## Features

- **Data Filtering** – Upload an expression matrix and group information, then filter genes based on expression level and sample counts.
- **Quality Control** – Visualize distributions of expression values and count statistics for selected samples or groups.
- **PCA & Correlation** – Perform principal component analysis and compute correlation matrices to assess sample relationships.
- **Pairwise Comparison** – Run differential expression tests between experimental groups using `DESeq2`.
- **Clustering** – Identify variable genes and display interactive heatmaps with hierarchical clustering or k‑means clusters.
- **Individual Gene** – Plot expression of specific genes across samples or groups.

Example expression and group tables are available under `example_data/` for demonstration.

## Getting Started

1. Install R (version 4.0 or higher is recommended).
2. Install required packages:

```r
install.packages(c(
  "shiny", "DT", "readxl", "readr", "dplyr", "ggplot2", "plotly",
  "pheatmap", "heatmaply", "DESeq2", "shinyWidgets", "logger"
))
```

3. Launch the application from the repository directory:

```r
source("app.R")
```

The Shiny interface will open in your browser. Upload an expression table (CSV or XLSX with columns `Symbol`, `Gene_Symbol`, then one column per sample) and a group table (CSV/XLSX with columns `Sample`, `Group`, `Color`). Follow the tabs to run analyses and download results.

Log files are written to the `logs/` directory each time the app starts.

## Repository Structure

- `app.R` – entry point that launches the Shiny app
- `ui.R` / `server.R` – top‑level Shiny definitions
- `R/` – module code for the individual analysis tabs
- `example_data/` – sample expression and group files

## License

This project is distributed under the [Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0) license](https://creativecommons.org/licenses/by-nc/4.0/). You may copy and modify the code for non-commercial purposes, provided that you give appropriate credit.

