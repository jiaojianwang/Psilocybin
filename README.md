# Psilocybin analysis code

This repository contains three analysis scripts: `Psilocybin_snRNAseq_analysis.R` for single-nucleus RNA-seq, `network_metrics.py` for functional network measures, and `permutation_lmm.m` for permutation-based linear mixed models. Update the data paths for your environment and install the packages used by each script before running them.

Direct dependencies

### R analysis

The following packages are imported directly by `Psilocybin_snRNAseq_analysis.R`.

| Package | Supplied version |
| --- | --- |
| AnnotationDbi | 1.66.0 |
| AUCell | 1.26.0 |
| Biobase | 2.64.0 |
| biomaRt | 2.60.1 |
| CellChat | 2.2.0 |
| clusterProfiler | 4.12.6 |
| CytoTRACE | 0.3.3 |
| decontX | 1.2.0 |
| DoubletFinder | 2.0.4 |
| dplyr | 1.1.4 |
| future | 1.34.0 |
| gghalves | 0.1.4 |
| ggplot2 | 3.5.1 |
| ggpubr | 0.6.0 |
| hdWGCNA | 0.4.00 |
| Matrix | 1.7-0 |
| monocle | 2.32.0 |
| org.Hs.eg.db | 3.19.1 |
| patchwork | 1.3.0 |
| Seurat | 5.1.0 |
| stringr | 1.5.1 |
| tidyr | 1.3.1 |

### Python analysis

| Package | Installed version | Used for |
| --- | --- | --- |
| NumPy | 1.26.4 | Array operations |
| nibabel | 5.1.0 | Neuroimaging file input |
| bctpy (`import bct`) | 0.6.1 | Brain Connectivity Toolbox network measures |

### MATLAB analysis

| Software | Installed version | Used for |
| --- | --- | --- |
| Statistics and Machine Learning Toolbox | 12.2 (R2021b) | Linear mixed-effects fitting with `fitlme` |
| Image Processing Toolbox | 11.4 (R2021b) | NIfTI input/output with `niftiread`, `niftiinfo`, and `niftiwrite` when NIfTI mode is used |
