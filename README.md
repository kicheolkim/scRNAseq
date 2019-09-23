# single-cell RNA-seq analysis using Seurat

Analysis code for single-cell RNA-Seq dataset using 10xGenomics platform, collaboration with Andrés Cruz-Herranz
Manuscript in preparation, result was presented in ECTRIMS 2018

## Analysis workflow
### 1. CellRanger for mapping and gene counting
- Mapping with reference genome, quantify gene expression (counting reads), and QC
- This is the script for run on QB3 cluster (UCSF cluster server)

### 2. Seurat (R package) for normalization and clustering
- Normalization of gene counts in each cell
- Discover cell clusters and identify cell types of each cluster.
- Comparison of change in cell number between different time-point.
- Differential gene expression analysis for cell sub-population.

-----
Created by Kicheol Kim, PhD (Feb, 2018)

Baranzini Lab. (https://github.com/baranzini-lab), Department of Neurology, UCSF
