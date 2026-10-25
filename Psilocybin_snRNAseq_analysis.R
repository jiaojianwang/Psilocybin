# Data loading and sample information ----------------------------------------------
library(Seurat)
library(decontX)
library(DoubletFinder)
library(Matrix)
library(ggplot2)

set.seed(123)
setwd("/home/yzy")

folder_path <- "./proj/medicine"
sample_dirs <- list.dirs(folder_path, full.names = TRUE, recursive = FALSE)
sample_names <- basename(sample_dirs)
data_dirs <- file.path(sample_dirs, "filtered_feature_bc_matrix")

stopifnot(length(sample_dirs) > 0, all(dir.exists(data_dirs)))

sample_objects <- lapply(seq_along(sample_names), function(i) {
  counts <- Read10X(data.dir = data_dirs[i])
  if (is.list(counts)) counts <- counts[["Gene Expression"]]
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_names[i],
    min.cells = 3,
    min.features = 200
  )
  obj <- RenameCells(obj, add.cell.id = sample_names[i])
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj
})
names(sample_objects) <- sample_names

# Ambient RNA decontamination ------------------------------------------------------
sample_objects <- lapply(sample_objects, function(obj) {
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  decontx_result <- decontX(counts)
  
  metadata <- obj@meta.data
  metadata$raw_nCount_RNA <- metadata$nCount_RNA
  metadata$raw_nFeature_RNA <- metadata$nFeature_RNA
  metadata$raw_percent.mt <- metadata$percent.mt
  metadata$Contamination <- decontx_result$contamination
  metadata <- metadata[, setdiff(
    colnames(metadata),
    c("nCount_RNA", "nFeature_RNA", "percent.mt")
  ), drop = FALSE]
  
  clean_obj <- CreateSeuratObject(
    counts = decontx_result$decontXcounts,
    project = Project(obj),
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
  clean_obj$percent.mt <- PercentageFeatureSet(clean_obj, pattern = "^mt-")
  clean_obj
})

# Quality control -----------------------------------------------------------------
sample_objects <- lapply(sample_objects, function(obj) {
  subset(
    obj,
    subset = raw_nCount_RNA <= 25000 &
      raw_nFeature_RNA <= 6000 &
      raw_percent.mt <= 5 &
      Contamination <= 0.4
  )
})

sample_sizes <- vapply(sample_objects, ncol, integer(1))
if (any(sample_sizes == 0L)) {
  stop(
    "No nuclei retained for sample(s): ",
    paste(names(sample_objects)[sample_sizes == 0L], collapse = ", ")
  )
}

# Sample-level doublet detection ---------------------------------------------------
doublet_rate_per_1000 <- 0.008
doublet_summary <- vector("list", length(sample_objects))

for (i in seq_along(sample_objects)) {
  sample_id <- names(sample_objects)[i]
  obj <- sample_objects[[i]]
  
  doublet_rate <- ncol(obj) / 1000 * doublet_rate_per_1000
  if (!is.finite(doublet_rate) || doublet_rate >= 1) {
    stop("Invalid expected doublet rate for sample: ", sample_id)
  }
  
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 3000
  )
  obj <- ScaleData(obj, features = VariableFeatures(obj))
  obj <- RunPCA(obj, features = VariableFeatures(obj))
  
  pca_variance <- Stdev(obj, reduction = "pca")^2
  cumulative_variance <- cumsum(pca_variance) / sum(pca_variance)
  dims_use <- seq_len(which(cumulative_variance >= 0.80)[1])
  
  obj <- FindNeighbors(obj, dims = dims_use)
  obj <- FindClusters(
    obj,
    algorithm = 4,
    resolution = 0.2,
    random.seed = 123
  )
  
  sweep_results <- paramSweep(obj, PCs = dims_use, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_results, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)
  bcmvn <- bcmvn[is.finite(bcmvn$BCmetric), , drop = FALSE]
  
  if (nrow(bcmvn) == 0L) {
    stop("No valid pK identified for sample: ", sample_id)
  }
  
  best_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic_prop <- modelHomotypic(as.character(obj$seurat_clusters))
  n_exp <- round(doublet_rate * ncol(obj))
  n_exp_adj <- round(n_exp * (1 - homotypic_prop))
  
  if (n_exp_adj > 0L) {
    old_columns <- colnames(obj@meta.data)
    obj <- doubletFinder(
      obj,
      PCs = dims_use,
      pN = 0.25,
      pK = best_pK,
      nExp = n_exp_adj,
      reuse.pANN = FALSE,
      sct = FALSE
    )
    
    new_columns <- setdiff(colnames(obj@meta.data), old_columns)
    class_column <- grep("^DF.classifications_", new_columns, value = TRUE)
    pann_column <- grep("^pANN_", new_columns, value = TRUE)
    stopifnot(length(class_column) == 1L, length(pann_column) == 1L)
    
    obj$Doublet <- obj@meta.data[[class_column]]
    obj$pANN <- obj@meta.data[[pann_column]]
    obj@meta.data <- obj@meta.data[, setdiff(
      colnames(obj@meta.data),
      c(class_column, pann_column)
    ), drop = FALSE]
  } else {
    obj$Doublet <- "Singlet"
    obj$pANN <- NA_real_
  }
  
  doublet_summary[[i]] <- data.frame(
    sample = sample_id,
    nuclei_before = ncol(obj),
    expected_doublets = n_exp,
    expected_heterotypic_doublets = n_exp_adj,
    detected_doublets = sum(obj$Doublet == "Doublet"),
    nuclei_retained = sum(obj$Doublet == "Singlet"),
    pK = best_pK
  )
  
  sample_objects[[i]] <- subset(obj, subset = Doublet == "Singlet")
  message("Finished: ", sample_id)
}

doublet_summary <- do.call(rbind, doublet_summary)
print(doublet_summary)

# Singlet retention and sample merging ---------------------------------------------
sample_objects <- lapply(sample_objects, function(obj) {
  stopifnot(!anyNA(obj$Doublet), all(obj$Doublet == "Singlet"))
  
  metadata <- obj@meta.data
  metadata <- metadata[, !grepl(
    "^(RNA_snn_res\\.|seurat_clusters$)", colnames(metadata)
  ), drop = FALSE]
  
  CreateSeuratObject(
    counts = GetAssayData(obj, assay = "RNA", layer = "counts"),
    project = Project(obj),
    meta.data = metadata,
    min.cells = 0,
    min.features = 0
  )
})

alldata <- if (length(sample_objects) == 1L) {
  sample_objects[[1]]
} else {
  merge(sample_objects[[1]], y = sample_objects[-1], merge.data = FALSE)
}

alldata <- JoinLayers(alldata, assay = "RNA")
alldata$orig.ident <- factor(alldata$orig.ident, levels = sample_names)

counts <- GetAssayData(alldata, assay = "RNA", layer = "counts")
retained_genes <- rownames(counts)[Matrix::rowSums(counts > 0) >= 3]
alldata <- subset(alldata, features = retained_genes)
counts <- GetAssayData(alldata, assay = "RNA", layer = "counts")

stopifnot(!anyNA(alldata$Doublet), all(alldata$Doublet == "Singlet"))
print(table(alldata$orig.ident, alldata$Doublet))


qc_summary <- data.frame(
  nuclei = ncol(alldata),
  genes = nrow(alldata),
  median_raw_genes = median(alldata$raw_nFeature_RNA),
  median_raw_UMIs = median(alldata$raw_nCount_RNA)
)
print(qc_summary)

# Cell-type marker visualization ---------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)

markers7 <- list(
  "L2-3" = c("CUX2", "ADAMTS12"),
  "L3-5" = c("KCNH5", "RORB"),
  "L5 ET" = c("VAT1L", "FEZF2"),
  "L5/6 IT" = c("HS3ST2", "THEMIS"),
  "L5/6 IT Car3" = c("NTNG2", "OPRK1"),
  "L5-6 NP" = c("NXPH2", "HTR2C"),
  "L6 CT" = c("SYT6", "SEMA3E"),
  "L6B" = c("ZFHX3", "NXPH4", "ARSJ"),
  "IN-LAMP5" = c("LAMP5", "ADARB2"),
  "IN-VIP" = "VIP",
  "IN-PV" = "PVALB",
  "IN-PV_ChC" = "UNC5B",
  "IN-SST" = "SST",
  "Astrocyte" = c("GFAP", "AQP4"),
  "Oligo" = c("MOBP", "MOG"),
  "OPC" = c("OLIG2", "PDGFRA"),
  "Microglia" = c("PTPRC", "FLI1", "APBB1IB", "INPP5D"),
  "Endothelial" = c("DCN", "COL1A2")
)

stopifnot("cell_type3" %in% colnames(scRNA_celltype@meta.data))
gene_annot <- stack(markers7)
colnames(gene_annot) <- c("gene", "group")
gene_annot <- gene_annot[!duplicated(gene_annot$gene), ]

missing_genes <- setdiff(gene_annot$gene, rownames(scRNA_celltype))
if (length(missing_genes) > 0) {
  message("Genes not found: ", paste(missing_genes, collapse = ", "))
}
gene_annot <- gene_annot[
  gene_annot$gene %in% rownames(scRNA_celltype), , drop = FALSE
]
stopifnot(nrow(gene_annot) > 0)

gene_order <- rev(gene_annot$gene)
group_order <- unique(as.character(gene_annot$group))
gene_annot$gene <- factor(gene_annot$gene, levels = gene_order)
gene_annot$group <- factor(gene_annot$group, levels = group_order)
stopifnot(all(group_order %in% names(celltype_col)))
group_colors <- celltype_col[group_order]

p_dot <- DotPlot(
  scRNA_celltype,
  features = gene_order,
  group.by = "cell_type3",
  assay = "RNA",
  dot.scale = 7
) +
  scale_x_discrete(limits = gene_order, expand = expansion(add = 0.6)) +
  scale_color_gradientn(
    colors = rev(c(
      "#810000", "#CE1212", "#F05454", "#FE8A7D", "#FFD06F",
      "#FFE6B7", "#72BCD5", "#528FAD", "#376795", "#1E466E"
    )),
    name = "Average expression"
  ) +
  labs(x = NULL, y = NULL, size = "Percent expressed") +
  coord_flip() +
  theme_classic(base_family = "serif") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

p_annotation <- ggplot(gene_annot, aes(x = 1, y = gene, fill = group)) +
  geom_tile(width = 0.7, height = 1) +
  scale_y_discrete(limits = gene_order, expand = expansion(add = 0.6)) +
  scale_x_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  scale_fill_manual(values = group_colors) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(5.5, 0, 5.5, 5.5)
  )

final_plot <- (p_annotation + p_dot) +
  plot_layout(widths = c(0.15, 5))
print(final_plot)

# Differential expression and functional enrichment --------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(123)

scRNA_celltype$group <- factor(
  scRNA_celltype$group,
  levels = c("WT", "Psilocybin")
)
DefaultAssay(scRNA_celltype) <- "RNA"
scRNA_celltype <- JoinLayers(scRNA_celltype, assay = "RNA")

Celltypes <- unique(as.character(scRNA_celltype$Celltype))
marker_list <- list()
print(table(scRNA_celltype$Celltype, scRNA_celltype$group))

for (ct in Celltypes) {
  cells_use <- colnames(scRNA_celltype)[scRNA_celltype$Celltype == ct]
  sub_seu <- subset(scRNA_celltype, cells = cells_use)
  group_counts <- table(sub_seu$group)
  
  if (any(group_counts < 3)) {
    message("Skipped: ", ct, "; fewer than three cells in one group.")
    next
  }
  
  markers <- FindMarkers(
    object = sub_seu,
    ident.1 = "Psilocybin",
    ident.2 = "WT",
    group.by = "group",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0
  )
  
  markers$Celltype <- rep(ct, nrow(markers))
  markers$gene <- rownames(markers)
  rownames(markers) <- NULL
  marker_list[[ct]] <- markers
  message("Finished differential expression: ", ct)
}

all_markers <- bind_rows(marker_list) %>%
  mutate(
    is_sig = !is.na(p_val_adj) &
      p_val_adj < 0.05 &
      abs(avg_log2FC) > 0.5
  )

significant_degs <- all_markers %>%
  group_by(Celltype) %>%
  summarise(
    up = sum(is_sig & avg_log2FC > 0.5, na.rm = TRUE),
    down = sum(is_sig & avg_log2FC < -0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(total = up + down) %>%
  arrange(desc(total), Celltype)

Celltype_order <- significant_degs$Celltype

bar_data <- significant_degs %>%
  pivot_longer(
    cols = c(up, down),
    names_to = "direction",
    values_to = "count"
  ) %>%
  mutate(
    direction = factor(direction, levels = c("up", "down"), labels = c("Up", "Down")),
    Celltype = factor(Celltype, levels = Celltype_order)
  )

all_markers <- all_markers %>%
  mutate(Celltype = factor(Celltype, levels = Celltype_order))

p_count <- ggplot(bar_data, aes(Celltype, count, fill = direction)) +
  geom_col(width = 0.55) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c("Up" = "#EF7677", "Down" = "#87B2D4")) +
  labs(x = NULL, y = "Number of DEGs", fill = NULL) +
  theme_classic(base_family = "serif") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "right"
  )

p_fc <- ggplot(all_markers, aes(Celltype, avg_log2FC)) +
  geom_hline(
    yintercept = c(-0.5, 0.5),
    linetype = "dashed", color = "grey65", linewidth = 0.3
  ) +
  geom_point(
    data = filter(all_markers, !is_sig),
    position = position_jitter(width = 0.3, height = 0, seed = 123),
    color = "grey80", size = 0.25, alpha = 0.7
  ) +
  geom_point(
    data = filter(all_markers, is_sig),
    aes(color = Celltype),
    position = position_jitter(width = 0.3, height = 0, seed = 123),
    size = 0.3, alpha = 0.8, show.legend = FALSE
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = Celltype_col) +
  labs(x = NULL, y = expression(log[2]~fold~change)) +
  theme_classic(base_family = "serif") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "none"
  )

combined_plot <- (p_count / p_fc) + plot_layout(heights = c(1, 3))
print(combined_plot)

org_db <- org.Hs.eg.db
kegg_organism <- "hsa"
all_results <- list()

gene_mapping <- bitr(
  geneID = unique(all_markers$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org_db
)

ratio_to_numeric <- function(x) {
  vapply(strsplit(x, "/", fixed = TRUE), function(z) {
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

for (ct in Celltype_order) {
  ct_markers <- all_markers %>%
    filter(as.character(Celltype) == ct)
  genes_up <- ct_markers %>%
    filter(is_sig, avg_log2FC > 0.5) %>%
    pull(gene)
  
  if (length(genes_up) == 0L) next
  
  background_ids <- unique(
    gene_mapping$ENTREZID[gene_mapping$SYMBOL %in% ct_markers$gene]
  )
  foreground_ids <- unique(
    gene_mapping$ENTREZID[gene_mapping$SYMBOL %in% genes_up]
  )
  
  if (length(foreground_ids) == 0L) next
  
  go_result <- enrichGO(
    gene = foreground_ids,
    universe = background_ids,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "ALL",
    pool = FALSE,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  
  kegg_result <- enrichKEGG(
    gene = foreground_ids,
    universe = background_ids,
    organism = kegg_organism,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  go_df <- NULL
  kegg_df <- NULL
  
  if (!is.null(go_result)) {
    go_df <- as.data.frame(go_result) %>%
      mutate(
        GeneRatio2 = ratio_to_numeric(GeneRatio),
        Change = "Positive",
        Celltype = ct,
        group_enrich = "GO"
      )
  }
  
  if (!is.null(kegg_result)) {
    kegg_result <- setReadable(
      kegg_result, OrgDb = org_db, keyType = "ENTREZID"
    )
    kegg_df <- as.data.frame(kegg_result) %>%
      mutate(
        ONTOLOGY = "KEGG",
        GeneRatio2 = ratio_to_numeric(GeneRatio),
        Change = "Positive",
        Celltype = ct,
        group_enrich = "KEGG"
      )
  }
  
  all_results[[ct]] <- bind_rows(go_df, kegg_df)
  message("Finished enrichment: ", ct)
}

final_enrich <- bind_rows(all_results)
final_enrich2 <- final_enrich %>%
  filter(!is.na(p.adjust), p.adjust < 0.05)

Celltype_plot_order <- c(
  "Microglia", "Oligo", "Astrocyte", "OPC", "Endothelial",
  "IN-LAMP5", "IN-PV", "IN-PV_ChC", "IN-SST", "IN-VIP",
  "L2-3", "L3-5", "L5 ET", "L5-6 NP", "L6 CT", "L6B",
  "L5/6 IT", "L5/6 IT Car3"
)
Celltype_plot_order <- c(
  Celltype_plot_order,
  setdiff(Celltype_order, Celltype_plot_order)
)

up_regulated <- final_enrich2 %>%
  filter(group_enrich == "GO") %>%
  mutate(Celltype = factor(Celltype, levels = Celltype_plot_order)) %>%
  group_by(Celltype, ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 2, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Celltype, p.adjust) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

p_go <- ggplot(up_regulated, aes(Celltype, Description)) +
  geom_point(aes(size = Count, color = p.adjust), alpha = 0.85) +
  scale_color_gradientn(
    colors = c(
      "#810000", "#CE1212", "#F05454", "#FE8A7D", "#FFD06F",
      "#FFE6B7", "#72BCD5", "#528FAD", "#376795", "#1E466E"
    )
  ) +
  scale_size_continuous(range = c(3, 8)) +
  labs(
    x = NULL, y = NULL,
    color = "Adjusted P", size = "Gene count",
    title = "GO enrichment of upregulated genes"
  ) +
  coord_flip() +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

print(p_go)
# Neuronal co-expression network analysis ------------------------------------------
library(Seurat)
library(hdWGCNA)
library(dplyr)
library(patchwork)

set.seed(123)

neuronal_Celltypes <- c(
  "IN-LAMP5", "IN-PV", "IN-PV_ChC", "IN-SST", "IN-VIP",
  "L2-3", "L3-5", "L5 ET", "L5-6 NP", "L6 CT", "L6B",
  "L5/6 IT", "L5/6 IT Car3"
)
network_Celltype <- "L3-5"

EX_IN <- subset(
  scRNA_celltype,
  subset = Celltype %in% neuronal_Celltypes
)
DefaultAssay(EX_IN) <- "RNA"
EX_IN <- JoinLayers(EX_IN, assay = "RNA")
EX_IN <- NormalizeData(
  EX_IN,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
EX_IN <- FindVariableFeatures(
  EX_IN,
  selection.method = "vst",
  nfeatures = 3000
)
EX_IN <- ScaleData(EX_IN, features = VariableFeatures(EX_IN))
EX_IN <- RunPCA(EX_IN, features = VariableFeatures(EX_IN))

seurat_obj <- SetupForWGCNA(
  EX_IN,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "neuron"
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Celltype", "orig.ident"),
  reduction = "pca",
  k = 25,
  max_shared = 10,
  ident.group = "Celltype",
  min_cells = 100
)
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = network_Celltype,
  group.by = "Celltype",
  assay = "RNA",
  slot = "data"
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed"
)
plot_list <- PlotSoftPowers(seurat_obj)
print(wrap_plots(plot_list, ncol = 2))

power_table <- GetPowerTable(seurat_obj)
print(power_table)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 6,
  networkType = "signed",
  tom_outdir = file.path(tempdir(), "neuron_TOM"),
  tom_name = "neuron"
)

PlotDendrogram(
  seurat_obj,
  main = "L3-5 neuronal hdWGCNA dendrogram"
)

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  assay = "RNA",
  group.by.vars = NULL
)

seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "Celltype",
  group_name = network_Celltype,
  assay = "RNA",
  harmonized = FALSE
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Module"
)

MEs <- GetMEs(seurat_obj, harmonized = FALSE)

p_kme <- PlotKMEs(
  seurat_obj,
  text_size = 1.5,
  ncol = 4,
  n_hubs = 10
)
print(p_kme)

modules <- GetModules(seurat_obj) %>%
  filter(module != "grey")

hub_df <- GetHubGenes(seurat_obj, n_hubs = 100)

head(MEs)
head(modules)
head(hub_df)
# Module functional enrichment -----------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)

enrich_GO_KEGG <- function(
    genes,
    group_name = "Set",
    OrgDb = org.Hs.eg.db,
    pvalue_cut = 1,
    qvalue_cut = 1,
    minGSSize = 5,
    maxGSSize = 500,
    organism = "hsa") {
  
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  gene_mapping <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = OrgDb
  )
  entrez_ids <- unique(gene_mapping$ENTREZID)
  
  go_result <- enrichGO(
    gene = entrez_ids,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    readable = TRUE
  )
  
  kegg_result <- enrichKEGG(
    gene = entrez_ids,
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  kegg_result <- setReadable(
    kegg_result,
    OrgDb = OrgDb,
    keyType = "ENTREZID"
  )
  
  go_df <- as.data.frame(go_result) %>%
    mutate(group = group_name, group_enrich = "GO")
  
  kegg_df <- as.data.frame(kegg_result) %>%
    mutate(
      ONTOLOGY = "KEGG",
      group = group_name,
      group_enrich = "KEGG"
    )
  
  bind_rows(go_df, kegg_df) %>%
    mutate(
      GeneRatio = as.character(GeneRatio),
      BgRatio = as.character(BgRatio),
      GeneRatio2 = vapply(
        strsplit(GeneRatio, "/", fixed = TRUE),
        function(x) as.numeric(x[1]) / as.numeric(x[2]),
        numeric(1)
      )
    ) %>%
    filter(
      pvalue <= pvalue_cut,
      is.na(qvalue) | qvalue <= qvalue_cut
    ) %>%
    arrange(group, group_enrich, p.adjust, desc(Count))
}

res <- enrich_GO_KEGG(
  genes = hub_df10$gene_name,
  group_name = "Module10",
  OrgDb = org.Hs.eg.db,
  pvalue_cut = 1,
  qvalue_cut = 1,
  minGSSize = 3,
  organism = "hsa"
)

GO <- res %>%
  filter(group_enrich == "GO", !is.na(p.adjust), p.adjust < 0.05)

use_pathway <- GO %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, desc(Count), ID, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("MF", "CC", "BP"))) %>%
  arrange(ONTOLOGY, desc(p.adjust), Count) %>%
  mutate(
    index = row_number(),
    enrichment_score = -log10(pmax(p.adjust, .Machine$double.xmin)),
    pathway_label = str_wrap(Description, width = 65),
    gene_label = str_wrap(str_replace_all(geneID, "/", ", "), width = 90)
  )

pal <- c(
  "BP" = "#7BC4E2",
  "CC" = "#ACD372",
  "MF" = "#FBB05B",
  "KEGG" = "#ED6CA4"
)

xaxis_max <- max(use_pathway$enrichment_score) * 1.05
annotation_width <- xaxis_max * 0.08

rect_data <- use_pathway %>%
  group_by(ONTOLOGY) %>%
  summarise(
    ymin = min(index) - 0.45,
    ymax = max(index) + 0.45,
    .groups = "drop"
  ) %>%
  mutate(
    xmin = -3 * annotation_width,
    xmax = -2 * annotation_width
  )

enrichment_plot <- ggplot(use_pathway) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = enrichment_score,
      ymin = index - 0.12,
      ymax = index + 0.38,
      fill = ONTOLOGY
    ),
    alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.02 * xaxis_max, y = index + 0.15, label = pathway_label),
    family = "serif",
    fontface = "bold",
    hjust = 0,
    size = 3.8
  ) +
  geom_text(
    aes(
      x = 0.02 * xaxis_max,
      y = index - 0.3,
      label = gene_label,
      color = ONTOLOGY
    ),
    family = "serif",
    fontface = "italic",
    hjust = 0,
    size = 2.8,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -annotation_width, y = index, size = Count, fill = ONTOLOGY),
    shape = 21,
    color = "grey30",
    stroke = 0.4,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = -annotation_width, y = index, label = Count),
    family = "serif",
    size = 3
  ) +
  geom_rect(
    data = rect_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = rect_data,
    aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = ONTOLOGY
    ),
    family = "serif",
    fontface = "bold",
    size = 3.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_color_manual(values = pal, guide = "none") +
  scale_size_continuous(range = c(7, 11)) +
  scale_x_continuous(
    breaks = pretty(c(0, xaxis_max), n = 5),
    expand = expansion(mult = c(0.02, 0.15))
  ) +
  scale_y_continuous(
    limits = c(0.4, nrow(use_pathway) + 0.6),
    expand = c(0, 0)
  ) +
  labs(
    x = expression(-log[10]~"(adjusted P)"),
    y = NULL,
    title = "Module10 GO enrichment"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "serif") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 11),
    axis.title.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(10, 30, 10, 10)
  )

print(enrichment_plot)

# HTR2A expression visualization ---------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggpubr)

plot_gene <- function(object, gene) {
  plot_df <- FetchData(
    object,
    vars = c(gene, "Celltype"),
    layer = "data"
  )
  colnames(plot_df)[colnames(plot_df) == gene] <- "Expression"
  
  mean_expression <- plot_df %>%
    group_by(Celltype) %>%
    summarise(
      mean_expression = mean(expm1(Expression), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_expression))
  
  plot_df <- plot_df %>%
    left_join(mean_expression, by = "Celltype") %>%
    mutate(Celltype = factor(Celltype, levels = mean_expression$Celltype))
  
  ggplot(plot_df, aes(Celltype, Expression, fill = mean_expression)) +
    geom_violin(scale = "width", alpha = 0.9, color = NA) +
    geom_boxplot(
      width = 0.25,
      color = "black",
      fill = NA,
      linewidth = 0.4,
      outlier.shape = NA
    ) +
    scale_fill_gradientn(
      colors = rev(c(
        "#810000", "#CE1212", "#F05454", "#FE8A7D", "#FFD06F",
        "#FFE6B7", "#72BCD5", "#528FAD", "#376795", "#1E466E"
      )),
      name = "Mean expression\n(linear scale)"
    ) +
    labs(title = gene, x = NULL, y = "Log-normalized expression") +
    theme_bw(base_size = 12, base_family = "serif") +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(
        angle = 45, hjust = 1, color = "black", size = 12
      ),
      axis.text.y = element_text(color = "black", size = 12),
      axis.ticks = element_blank(),
      legend.position = "right"
    )
}

DefaultAssay(scRNA_celltype) <- "RNA"
scRNA_celltype <- JoinLayers(scRNA_celltype, assay = "RNA")

p_all <- plot_gene(scRNA_celltype, gene = "HTR2A")
print(p_all)

post_dose <- subset(scRNA_celltype, subset = plot_group == "Post-dose")
p_post <- plot_gene(post_dose, gene = "HTR2A")
print(p_post)

meta_data <- FetchData(
  scRNA_celltype,
  vars = c("HTR2A", "Celltype", "plot_group"),
  layer = "data"
) %>%
  filter(plot_group %in% c("Pre-dose", "Post-dose")) %>%
  mutate(plot_group = factor(plot_group, levels = c("Pre-dose", "Post-dose")))

expression_summary <- meta_data %>%
  group_by(plot_group, Celltype) %>%
  summarise(
    mean_expression = mean(HTR2A, na.rm = TRUE),
    median_expression = median(HTR2A, na.rm = TRUE),
    .groups = "drop"
  )

plot_df <- meta_data %>%
  filter(Celltype == "L3-5")

median_df <- expression_summary %>%
  filter(Celltype == "L3-5")

group_colors <- c(
  "Pre-dose" = "#569EA4",
  "Post-dose" = "#FCCDE5"
)
background_df <- data.frame(
  xmin = c(0.5, 1.5),
  xmax = c(1.5, 2.5),
  background = c("#FFE6B7", "#AADB84")
)

expression_range <- range(plot_df$HTR2A, na.rm = TRUE)
label_y <- expression_range[2] + 0.12 * max(diff(expression_range), 0.1)

p_comparison <- ggplot(
  plot_df,
  aes(plot_group, HTR2A, fill = plot_group)
) +
  geom_rect(
    data = background_df,
    aes(
      xmin = xmin, xmax = xmax,
      ymin = -Inf, ymax = Inf,
      fill = background
    ),
    inherit.aes = FALSE,
    alpha = 0.15
  ) +
  geom_vline(
    xintercept = 1.5,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.4
  ) +
  geom_half_violin(
    side = "r",
    scale = "width",
    color = NA,
    alpha = 0.85
  ) +
  geom_boxplot(
    width = 0.25,
    position = position_nudge(x = -0.18),
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  geom_line(
    data = median_df,
    aes(plot_group, median_expression, group = 1),
    inherit.aes = FALSE,
    color = "#810000",
    linewidth = 0.9
  ) +
  geom_point(
    data = median_df,
    aes(plot_group, median_expression),
    inherit.aes = FALSE,
    color = "#810000",
    size = 2
  ) +
  stat_compare_means(
    comparisons = list(c("Pre-dose", "Post-dose")),
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    paired = FALSE,
    label = "p.signif",
    hide.ns = FALSE,
    tip.length = 0,
    bracket.size = 0.4,
    label.y = label_y
  ) +
  scale_fill_manual(
    values = c(
      group_colors,
      "#FFE6B7" = "#FFE6B7",
      "#AADB84" = "#AADB84"
    ),
    guide = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x = NULL,
    y = "Log-normalized expression",
    title = "L3-5 HTR2A expression"
  ) +
  theme_classic(base_size = 12, base_family = "serif") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    axis.text.x = element_text(
      angle = 45, hjust = 1, color = "black", size = 12
    ),
    axis.text.y = element_text(color = "black", size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

print(expression_summary)
print(p_comparison)
# Synaptic gene-set activity scoring ------------------------------------------------
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(AUCell)
library(ggplot2)

set.seed(123)

go_terms <- c(
  "GO:0007268",
  "GO:0050808",
  "GO:0007269",
  "GO:0048169"
)

synaptic_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_terms,
  columns = c("ENTREZID", "SYMBOL", "GENENAME"),
  keytype = "GO"
)

DefaultAssay(scRNA_celltype) <- "RNA"
scRNA_celltype <- JoinLayers(scRNA_celltype, assay = "RNA")
expression_matrix <- GetAssayData(
  scRNA_celltype,
  assay = "RNA",
  layer = "data"
)

synaptic_genes <- unique(na.omit(synaptic_annotation$SYMBOL))
synaptic_genes <- intersect(synaptic_genes, rownames(expression_matrix))

cells_rankings <- AUCell_buildRankings(
  expression_matrix,
  splitByBlocks = TRUE,
  plotStats = FALSE
)

cells_AUC <- AUCell_calcAUC(
  geneSets = list(Synaptic = synaptic_genes),
  rankings = cells_rankings,
  aucMaxRank = ceiling(nrow(expression_matrix) * 0.1)
)

auc_matrix <- getAUC(cells_AUC)
scRNA_celltype$Synaptic <- as.numeric(
  auc_matrix["Synaptic", colnames(scRNA_celltype)]
)

p_synaptic <- FeaturePlot(
  scRNA_celltype,
  features = "Synaptic",
  reduction = "umap",
  pt.size = 0.01,
  order = TRUE,
  label = FALSE,
  raster = FALSE
) +
  scale_color_gradientn(
    colors = rev(c(
      "#810000", "#CE1212", "#F05454",
      "#FFD06F", "#FFE6B7", "grey90"
    ))
  ) +
  labs(title = "Synaptic gene-set activity", color = "AUC") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(family = "serif")
  )

print(p_synaptic)
# BDNF expression heatmap ----------------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)

DefaultAssay(scRNA_celltype) <- "RNA"
scRNA_celltype <- JoinLayers(scRNA_celltype, assay = "RNA")

plot_df <- FetchData(
  scRNA_celltype,
  vars = c("BDNF", "Celltype", "plot_group"),
  layer = "data"
) %>%
  filter(
    grepl("^(L|IN-)", Celltype),
    plot_group %in% c("Pre-dose", "Post-dose")
  ) %>%
  mutate(
    Celltype = as.character(Celltype),
    plot_group = factor(plot_group, levels = c("Pre-dose", "Post-dose"))
  )

df_summary <- plot_df %>%
  group_by(Celltype, plot_group) %>%
  summarise(
    mean_expression = mean(BDNF, na.rm = TRUE),
    .groups = "drop"
  )

Celltype_order <- df_summary %>%
  filter(plot_group == "Post-dose") %>%
  arrange(mean_expression, Celltype) %>%
  pull(Celltype)

Celltype_order <- c(
  Celltype_order,
  setdiff(unique(df_summary$Celltype), Celltype_order)
)

df_summary <- df_summary %>%
  mutate(Celltype = factor(Celltype, levels = Celltype_order))

base_colors <- c(
  "#FFFFFF", "#FFEDC2", "#FFE4B2", "#FFD984", "#FFD06F",
  "#FFB84A", "#F4A582", "#E06651", "#D04A4A", "#C02D3D",
  "#B2182B", "#8E1A2A", "#67001F"
)

p_bdnf <- ggplot(
  df_summary,
  aes(x = plot_group, y = Celltype, fill = mean_expression)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colors = colorRampPalette(base_colors)(256),
    name = "Mean log-normalized\nexpression"
  ) +
  scale_y_discrete(position = "right") +
  coord_equal() +
  labs(x = NULL, y = NULL, title = "BDNF") +
  theme_minimal(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      size = 14, angle = 90, hjust = 1, vjust = 0.5,
      face = "bold", color = "black"
    ),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

print(p_bdnf)
# Neuroinflammation gene-set activity scoring --------------------------------------
library(Seurat)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(AUCell)
library(ggplot2)

set.seed(123)

DefaultAssay(scRNA_celltype) <- "RNA"
scRNA_celltype <- JoinLayers(scRNA_celltype, assay = "RNA")
scRNA_celltype <- NormalizeData(
  scRNA_celltype,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

GeneCards_genes <- read.csv(
  "/home/yzy/proj/medicine/V1/GeneCards-SearchResults.csv",
  header = TRUE
)
GeneCards_symbols <- with(
  GeneCards_genes,
  Gene.Symbol[!is.na(Relevance.score) & Relevance.score > 3]
)

neuroinflammation_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = c("GO:0150078", "GO:0150076"),
  columns = c("ENTREZID", "SYMBOL", "GENENAME"),
  keytype = "GO"
)

expression_matrix <- GetAssayData(
  scRNA_celltype,
  assay = "RNA",
  layer = "data"
)

neuroinflammation_genes <- unique(c(
  neuroinflammation_annotation$SYMBOL,
  GeneCards_symbols
))
neuroinflammation_genes <- neuroinflammation_genes[
  !is.na(neuroinflammation_genes) & nzchar(neuroinflammation_genes)
]
neuroinflammation_genes <- intersect(
  neuroinflammation_genes,
  rownames(expression_matrix)
)

cells_rankings <- AUCell_buildRankings(
  expression_matrix,
  splitByBlocks = TRUE,
  plotStats = FALSE
)
cells_AUC <- AUCell_calcAUC(
  geneSets = list(Neuroinflammation = neuroinflammation_genes),
  rankings = cells_rankings,
  aucMaxRank = ceiling(nrow(expression_matrix) * 0.1)
)

auc_matrix <- getAUC(cells_AUC)
scRNA_celltype$Neuroinflammation <- as.numeric(
  auc_matrix["Neuroinflammation", colnames(scRNA_celltype)]
)

p_neuroinflammation <- FeaturePlot(
  scRNA_celltype,
  features = "Neuroinflammation",
  reduction = "umap",
  pt.size = 0.01,
  order = TRUE,
  label = FALSE,
  raster = FALSE
) +
  scale_color_gradientn(
    colors = rev(c(
      "#810000", "#CE1212", "#F05454",
      "#FFD06F", "#FFE6B7", "grey90"
    ))
  ) +
  labs(title = "Neuroinflammation gene-set activity", color = "AUC") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(family = "serif")
  )

print(p_neuroinflammation)
# Microglial subclustering and gene-set activity scoring ----------------------------
.libPaths(c("/home/yzy/docker_library", .libPaths()))
library(Seurat)
library(AUCell)
library(ggplot2)
library(patchwork)

set.seed(123)

scRNA_celltype <- readRDS("/home/yzy/proj/medicine/V2/scRNA_celltype.rds")
Microglia <- subset(scRNA_celltype, subset = Celltype == "Microglia")
DefaultAssay(Microglia) <- "RNA"
Microglia <- JoinLayers(Microglia, assay = "RNA")
Microglia <- NormalizeData(
  Microglia,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

sample_cells <- split(
  colnames(Microglia),
  as.character(Microglia$orig.ident),
  drop = TRUE
)
variable_genes <- lapply(sample_cells, function(cells_use) {
  obj <- subset(Microglia, cells = cells_use)
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  VariableFeatures(obj)
})
VariableFeatures(Microglia) <- unique(unlist(variable_genes, use.names = FALSE))

Microglia <- ScaleData(
  Microglia,
  features = VariableFeatures(Microglia),
  verbose = FALSE
)
Microglia <- RunPCA(
  Microglia,
  features = VariableFeatures(Microglia),
  npcs = 30,
  seed.use = 123,
  verbose = FALSE
)
Microglia <- FindNeighbors(Microglia, reduction = "pca", dims = 1:20)
Microglia <- FindClusters(Microglia, resolution = 0.5, random.seed = 123)
Microglia <- RunUMAP(
  Microglia,
  reduction = "pca",
  dims = 1:20,
  seed.use = 123
)

p_cluster <- DimPlot(
  Microglia,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  raster = FALSE
)
p_sample <- DimPlot(
  Microglia,
  reduction = "umap",
  group.by = "orig.ident",
  raster = FALSE
)
print(p_cluster)
print(p_sample)

gene_sets <- list(
  M1 = c(
    "CD40", "CXCL9", "MYD88", "TLR2", "CXCL11", "TRAF6", "MAPK8", "TNF",
    "RFXANK", "MAP3K7", "IRF1", "IL1A", "IRF5", "CCR7", "PTGS2", "SOAT1",
    "CD80", "CD36", "ROCK2", "CRYGC", "IRAK4", "NLRP3", "STAT3", "CXCL10",
    "CD47", "FOSB", "KYNU", "NLRC3", "NOS2", "TREM1", "IL23R", "CCL5",
    "FCGR3A", "TLR4", "IRAK1", "IL1B", "MAPK14", "IDO1", "FCGR2B", "BCL10",
    "FCGR2A", "MAPK1", "IL17A", "IL18", "DECR1", "ITGAX", "PTPRC", "CD14",
    "STAT1", "ITGAM", "IL6", "TAOK2", "MAST2", "CD86"
  ),
  M2 = c(
    "MFGE8", "GDNF", "IL10", "RETNLB", "NLRP12", "CTSA", "FASLG", "CD163",
    "MMP14", "TGFB2", "TGFB3", "VEGFA", "STK11", "CCL4", "MSR1", "MMP9",
    "NGF", "IL13", "FN1", "FIGF", "CD276", "IL4", "SOCS1", "TNFSF12",
    "STAT6", "PRKAB1", "BDNF", "IL33", "WNT7B", "CCL18", "CCL24", "GCM2",
    "LYVE1", "IRF4", "HTRA1", "CTSB", "RUNX2", "GH1", "MMP19", "EGF",
    "TNFSF8", "IL4R", "SOCS3", "GRN", "CLEC7A", "VTCN1", "ARG1", "CTSC",
    "NLRC3", "CCL13", "CCL17", "VEGFC", "MRC1", "VEGFB", "TGFB1", "CCL20",
    "CTSD", "KLK1", "TREM2", "AP1B1", "CCL22"
  )
)

expression_matrix <- GetAssayData(
  Microglia,
  assay = "RNA",
  layer = "data"
)

gene_set_summary <- data.frame(
  Gene_set = names(gene_sets),
  Input_genes = lengths(gene_sets),
  Matched_genes = vapply(
    gene_sets,
    function(genes) length(intersect(genes, rownames(expression_matrix))),
    integer(1)
  )
)
print(gene_set_summary)

gene_sets <- lapply(gene_sets, function(genes) {
  intersect(unique(genes), rownames(expression_matrix))
})

cells_rankings <- AUCell_buildRankings(
  expression_matrix,
  splitByBlocks = TRUE,
  plotStats = FALSE
)
cells_AUC <- AUCell_calcAUC(
  geneSets = gene_sets,
  rankings = cells_rankings,
  aucMaxRank = ceiling(nrow(expression_matrix) * 0.1)
)

auc_matrix <- getAUC(cells_AUC)
Microglia$M1_score <- as.numeric(auc_matrix["M1", colnames(Microglia)])
Microglia$M2_score <- as.numeric(auc_matrix["M2", colnames(Microglia)])

p_scores <- FeaturePlot(
  Microglia,
  features = c("M1_score", "M2_score"),
  reduction = "umap",
  pt.size = 0.01,
  ncol = 2,
  order = TRUE,
  label = FALSE,
  raster = FALSE
) &
  scale_color_gradientn(
    colors = rev(c(
      "#810000", "#CE1212", "#F05454",
      "#FFD06F", "#FFE6B7", "grey90"
    ))
  ) &
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(family = "serif")
  )

print(p_scores)
# Human ortholog conversion --------------------------------------------------------
library(Seurat)
library(biomaRt)
library(Matrix)

convert_to_human_homologs <- function(
    seurat_obj,
    from = c("macaque", "mouse"),
    macaque_dataset = "mmulatta_gene_ensembl",
    ensembl_version = NULL,
    batch_size = 200) {
  
  from <- match.arg(from)
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  gene_list <- rownames(counts)
  gene_list <- unique(gene_list[!is.na(gene_list) & nzchar(gene_list)])
  
  dataset <- switch(
    from,
    macaque = macaque_dataset,
    mouse = "mmusculus_gene_ensembl"
  )
  
  mart <- useEnsembl(
    biomart = "genes",
    dataset = dataset,
    version = ensembl_version
  )
  
  gene_batches <- split(
    gene_list,
    ceiling(seq_along(gene_list) / batch_size)
  )
  
  mapping_list <- lapply(gene_batches, function(genes) {
    getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "hsapiens_homolog_ensembl_gene",
        "hsapiens_homolog_associated_gene_name",
        "hsapiens_homolog_orthology_type"
      ),
      filters = "external_gene_name",
      values = genes,
      mart = mart,
      useCache = FALSE
    )
  })
  
  homolog_map <- unique(do.call(rbind, mapping_list))
  colnames(homolog_map) <- c(
    "source_ensembl",
    "source_gene",
    "human_ensembl",
    "human_gene",
    "orthology_type"
  )
  
  homolog_map <- homolog_map[
    !is.na(homolog_map$source_gene) &
      nzchar(homolog_map$source_gene) &
      !is.na(homolog_map$human_gene) &
      nzchar(homolog_map$human_gene) &
      homolog_map$orthology_type == "ortholog_one2one",
    ,
    drop = FALSE
  ]
  
  ambiguous_source <- duplicated(homolog_map$source_gene) |
    duplicated(homolog_map$source_gene, fromLast = TRUE)
  ambiguous_target <- duplicated(homolog_map$human_gene) |
    duplicated(homolog_map$human_gene, fromLast = TRUE)
  
  homolog_map <- homolog_map[
    !ambiguous_source & !ambiguous_target,
    ,
    drop = FALSE
  ]
  homolog_map <- homolog_map[
    match(
      intersect(rownames(counts), homolog_map$source_gene),
      homolog_map$source_gene
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(homolog_map) == 0L) {
    stop("No unambiguous one-to-one human orthologs matched the input genes.")
  }
  
  human_counts <- counts[homolog_map$source_gene, , drop = FALSE]
  rownames(human_counts) <- homolog_map$human_gene
  
  metadata <- seurat_obj@meta.data
  metadata$source_nCount_RNA <- metadata$nCount_RNA
  metadata$source_nFeature_RNA <- metadata$nFeature_RNA
  
  if ("percent.mt" %in% colnames(metadata)) {
    metadata$source_percent.mt <- metadata$percent.mt
  }
  
  metadata <- metadata[, setdiff(
    colnames(metadata),
    c("nCount_RNA", "nFeature_RNA", "percent.mt")
  ), drop = FALSE]
  
  human_obj <- CreateSeuratObject(
    counts = human_counts,
    meta.data = metadata,
    project = paste0(Project(seurat_obj), "_human"),
    min.cells = 0,
    min.features = 0
  )
  
  human_obj@misc$ortholog_conversion <- list(
    source_species = from,
    source_dataset = dataset,
    ensembl_version = ensembl_version,
    target_species = "human",
    orthology_type = "ortholog_one2one",
    input_genes = length(gene_list),
    mapped_genes = nrow(human_obj),
    mapping = homolog_map
  )
  
  message(
    "Conversion complete: ",
    nrow(human_obj), " human genes retained from ",
    length(gene_list), " input genes."
  )
  
  human_obj
}

monkey <- convert_to_human_homologs(
  seurat_obj = scRNA_celltype,
  from = "macaque",
  macaque_dataset = "mmulatta_gene_ensembl"
)

conversion_summary <- monkey@misc$ortholog_conversion
conversion_summary$mapped_genes
conversion_summary$source_dataset
head(conversion_summary$mapping)
dim(monkey)

# Cell–cell communication analysis -------------------------------------------------
library(Seurat)
library(CellChat)
library(future)

set.seed(123)
load("/home/yzy/proj/medicine/V2/monkey.RData")

print(table(monkey$Celltype, monkey$group))

run_cellchat <- function(
    seurat_obj,
    species = c("hsa", "mmu"),
    signal_type = c("All", "Secreted", "ECM", "Contact"),
    cell_group = "Celltype",
    workers = 1,
    comm_method = "triMean",
    min_cells = 10) {
  
  species <- match.arg(species)
  signal_type <- match.arg(signal_type)
  
  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan("multisession", workers = workers)
  
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  data_input <- GetAssayData(
    seurat_obj,
    assay = "RNA",
    layer = "data"
  )
  metadata <- seurat_obj@meta.data[
    colnames(data_input),
    c("group", cell_group),
    drop = FALSE
  ]
  metadata[[cell_group]] <- droplevels(factor(metadata[[cell_group]]))
  
  cellchat <- createCellChat(
    object = data_input,
    meta = metadata,
    group.by = cell_group
  )
  
  database <- switch(
    species,
    hsa = CellChatDB.human,
    mmu = CellChatDB.mouse
  )
  cellchat@DB <- switch(
    signal_type,
    All = database,
    Secreted = subsetDB(database, search = "Secreted Signaling"),
    ECM = subsetDB(database, search = "ECM-Receptor"),
    Contact = subsetDB(database, search = "Cell-Cell Contact")
  )
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(
    cellchat,
    type = comm_method,
    raw.use = TRUE,
    population.size = FALSE,
    seed.use = 123
  )
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat
}

WT <- subset(monkey, subset = group == "WT")
Psilocybin <- subset(monkey, subset = group == "Psilocybin")

cellchat_WT <- run_cellchat(
  seurat_obj = WT,
  species = "hsa",
  signal_type = "All",
  cell_group = "Celltype"
)

cellchat_Psilocybin <- run_cellchat(
  seurat_obj = Psilocybin,
  species = "hsa",
  signal_type = "All",
  cell_group = "Celltype"
)

communication_WT <- subsetCommunication(cellchat_WT)
communication_Psilocybin <- subsetCommunication(cellchat_Psilocybin)

head(communication_WT)
head(communication_Psilocybin)
# Neuron–glia communication differences --------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

reference_group <- "Pre-dose"
comparison_group <- "Post-dose"

neuronal_Celltypes <- c(
  "L3-5", "L5 ET", "L5-6 NP", "L5/6 IT", "L5/6 IT Car3"
)

extract_sender_receiver <- function(
    object_list,
    senders,
    receivers,
    pathways = NULL) {
  
  if (is.null(pathways)) {
    pathways <- sort(unique(unlist(
      lapply(object_list, function(obj) obj@netP$pathways),
      use.names = FALSE
    )))
  }
  
  results <- lapply(names(object_list), function(group_name) {
    obj <- object_list[[group_name]]
    probability <- obj@netP$prob
    
    available_senders <- intersect(senders, dimnames(probability)[[1]])
    available_receivers <- intersect(receivers, dimnames(probability)[[2]])
    
    if (length(available_senders) == 0L ||
        length(available_receivers) == 0L) {
      return(NULL)
    }
    
    pathway_results <- lapply(pathways, function(pathway_name) {
      if (pathway_name %in% dimnames(probability)[[3]]) {
        probability_matrix <- probability[
          available_senders,
          available_receivers,
          pathway_name,
          drop = FALSE
        ]
        dim(probability_matrix) <- c(
          length(available_senders),
          length(available_receivers)
        )
        dimnames(probability_matrix) <- list(
          sender = available_senders,
          receiver = available_receivers
        )
        
        result <- as.data.frame.table(
          probability_matrix,
          responseName = "value",
          stringsAsFactors = FALSE
        )
      } else {
        result <- expand_grid(
          sender = available_senders,
          receiver = available_receivers
        ) %>%
          mutate(value = 0)
      }
      
      result %>%
        mutate(pathway = pathway_name, group = group_name)
    })
    
    bind_rows(pathway_results)
  })
  
  bind_rows(results)
}

calculate_difference <- function(
    object_list,
    senders,
    receivers,
    reference_group,
    comparison_group) {
  
  selected_objects <- object_list[c(reference_group, comparison_group)]
  
  extract_sender_receiver(
    object_list = selected_objects,
    senders = senders,
    receivers = receivers
  ) %>%
    pivot_wider(
      id_cols = c(sender, receiver, pathway),
      names_from = group,
      values_from = value
    ) %>%
    mutate(
      delta = .data[[comparison_group]] - .data[[reference_group]]
    ) %>%
    filter(!is.na(delta))
}

communication_to_astrocyte <- calculate_difference(
  object.list,
  senders = neuronal_Celltypes,
  receivers = "Astrocyte",
  reference_group = reference_group,
  comparison_group = comparison_group
)

communication_from_astrocyte <- calculate_difference(
  object.list,
  senders = "Astrocyte",
  receivers = neuronal_Celltypes,
  reference_group = reference_group,
  comparison_group = comparison_group
)

communication_to_microglia <- calculate_difference(
  object.list,
  senders = neuronal_Celltypes,
  receivers = "Microglia",
  reference_group = reference_group,
  comparison_group = comparison_group
)

communication_from_microglia <- calculate_difference(
  object.list,
  senders = "Microglia",
  receivers = neuronal_Celltypes,
  reference_group = reference_group,
  comparison_group = comparison_group
)

all_pathways <- sort(unique(c(
  communication_to_astrocyte$pathway,
  communication_from_astrocyte$pathway,
  communication_to_microglia$pathway,
  communication_from_microglia$pathway
)))

base_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
pathway_colors <- setNames(
  colorRampPalette(base_colors)(length(all_pathways)),
  all_pathways
)

plot_communication_difference <- function(
    data,
    facet_by,
    title,
    Celltype_order = neuronal_Celltypes) {
  
  plot_data <- data %>%
    mutate(Celltype = factor(.data[[facet_by]], levels = Celltype_order))
  
  pathway_order <- plot_data %>%
    group_by(pathway) %>%
    summarise(mean_delta = mean(delta), .groups = "drop") %>%
    arrange(mean_delta, pathway) %>%
    pull(pathway)
  
  plot_data <- plot_data %>%
    mutate(pathway = factor(pathway, levels = pathway_order))
  
  ggplot(plot_data, aes(x = delta, y = pathway)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#E64B35",
      linewidth = 0.5
    ) +
    geom_segment(
      aes(x = 0, xend = delta, yend = pathway),
      color = "grey65",
      linewidth = 0.6
    ) +
    geom_point(aes(color = pathway), size = 2.5) +
    scale_color_manual(values = pathway_colors, guide = "none") +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    facet_wrap(~Celltype, ncol = 5) +
    labs(
      x = paste0(
        "Communication probability difference (",
        comparison_group, " − ", reference_group, ")"
      ),
      y = NULL,
      title = title
    ) +
    theme_bw(base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(
        color = "black", fill = NA, linewidth = 0.6
      ),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(
        angle = 45, hjust = 1, color = "black", size = 10
      ),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.x = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
}

p_to_astrocyte <- plot_communication_difference(
  communication_to_astrocyte,
  facet_by = "sender",
  title = "Neurons → Astrocyte"
)

p_from_astrocyte <- plot_communication_difference(
  communication_from_astrocyte,
  facet_by = "receiver",
  title = "Astrocyte → Neurons"
)

p_to_microglia <- plot_communication_difference(
  communication_to_microglia,
  facet_by = "sender",
  title = "Neurons → Microglia"
)

p_from_microglia <- plot_communication_difference(
  communication_from_microglia,
  facet_by = "receiver",
  title = "Microglia → Neurons"
)

print(p_to_astrocyte)
print(p_from_astrocyte)
print(p_to_microglia)
print(p_from_microglia)
# Astrocyte trajectory and state-specific analysis ---------------------------------
library(Seurat)
library(monocle)
library(Biobase)
library(CytoTRACE)
library(Matrix)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(123)

scRNA_celltype <- readRDS(
  "/data3/yangzeyuan/project_brain/medicine/v2/scRNA_celltype.rds"
)
seu <- subset(scRNA_celltype, subset = Celltype == "Astrocyte")
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu, assay = "RNA")
seu <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 2000
)

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
retained_genes <- Matrix::rowSums(counts > 0) >= 5
mono_matrix <- counts[retained_genes, , drop = FALSE]

feature_info <- data.frame(
  gene_short_name = rownames(mono_matrix),
  row.names = rownames(mono_matrix)
)
cell_info <- seu@meta.data[colnames(mono_matrix), , drop = FALSE]

cds <- newCellDataSet(
  mono_matrix,
  phenoData = new("AnnotatedDataFrame", data = cell_info),
  featureData = new("AnnotatedDataFrame", data = feature_info),
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

ordering_genes <- intersect(VariableFeatures(seu), rownames(cds))
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(
  cds,
  max_components = 2,
  reduction_method = "DDRTree"
)
cds <- orderCells(cds)

trajectory_metadata <- Biobase::pData(cds)
seu$State <- factor(
  trajectory_metadata[colnames(seu), "State"]
)
seu$Pseudotime <- trajectory_metadata[colnames(seu), "Pseudotime"]

print(table(seu$State, seu$group))

p_group <- plot_cell_trajectory(cds, color_by = "group") +
  scale_color_manual(values = c("WT" = "#4DBBD5", "Psilocybin" = "#E64B35")) +
  theme(
    text = element_text(family = "serif"),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

p_state <- plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap(~group, nrow = 1) +
  theme(
    text = element_text(family = "serif"),
    axis.ticks = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

p_pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
  theme(text = element_text(family = "serif"))

print(p_group)
print(p_state)
print(p_pseudotime)

cytotrace_results <- CytoTRACE(
  as.matrix(mono_matrix),
  enableFast = FALSE
)
seu$CytoTRACE <- as.numeric(
  cytotrace_results$CytoTRACE[colnames(seu)]
)

p_cytotrace <- plot_cell_trajectory(cds, color_by = "State")
Biobase::pData(cds)$CytoTRACE <- seu@meta.data[
  colnames(cds), "CytoTRACE"
]

p_cytotrace <- plot_cell_trajectory(cds, color_by = "CytoTRACE") +
  scale_color_gradientn(
    colors = rev(c(
      "#810000", "#CE1212", "#F05454", "#FE8A7D", "#FFD06F",
      "#FFE6B7", "#72BCD5", "#528FAD", "#376795", "#1E466E"
    ))
  ) +
  theme(text = element_text(family = "serif"))

print(p_cytotrace)

cytotrace_summary <- seu@meta.data %>%
  group_by(State) %>%
  summarise(
    nuclei = n(),
    mean_CytoTRACE = mean(CytoTRACE, na.rm = TRUE),
    median_CytoTRACE = median(CytoTRACE, na.rm = TRUE),
    .groups = "drop"
  )
print(cytotrace_summary)

state_counts <- with(
  seu@meta.data,
  table(plot_group, State)
)
state_counts <- state_counts[
  rowSums(state_counts) > 0,
  colSums(state_counts) > 0,
  drop = FALSE
]
expected_counts <- outer(
  rowSums(state_counts),
  colSums(state_counts)
) / sum(state_counts)

roe_matrix <- state_counts / expected_counts
roe_data <- as.data.frame.table(
  roe_matrix,
  responseName = "ROE",
  stringsAsFactors = FALSE
)
colnames(roe_data)[1:2] <- c("plot_group", "State")

p_roe <- ggplot(roe_data, aes(State, plot_group, fill = ROE)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", ROE)), size = 4) +
  scale_x_discrete(
    labels = function(x) paste0("State-", x),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_gradient2(
    low = "#4DBBD5",
    mid = "white",
    high = "#E64B35",
    midpoint = 1
  ) +
  labs(x = NULL, y = NULL, fill = "O/E", title = "State enrichment") +
  theme_minimal(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14)
  )

print(p_roe)

marker_list <- list()

for (state_id in levels(seu$State)) {
  cells_use <- colnames(seu)[seu$State == state_id]
  sub_seu <- subset(seu, cells = cells_use)
  group_counts <- table(
    factor(sub_seu$group, levels = c("WT", "Psilocybin"))
  )
  
  if (any(group_counts < 3)) {
    message("Skipped State-", state_id, "; insufficient cells.")
    next
  }
  
  markers <- FindMarkers(
    sub_seu,
    ident.1 = "Psilocybin",
    ident.2 = "WT",
    group.by = "group",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0
  )
  
  markers$State <- rep(state_id, nrow(markers))
  markers$Celltype <- rep("Astrocyte", nrow(markers))
  markers$gene <- rownames(markers)
  rownames(markers) <- NULL
  marker_list[[state_id]] <- markers
}

all_markers <- bind_rows(marker_list)

significant_degs <- all_markers %>%
  group_by(Celltype, State) %>%
  summarise(
    up = sum(p_val_adj < 0.05 & avg_log2FC > 0.5, na.rm = TRUE),
    down = sum(p_val_adj < 0.05 & avg_log2FC < -0.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(total = up + down)

all_markers_up <- all_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5)

print(significant_degs)

gene_mapping <- bitr(
  unique(all_markers$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

all_results <- list()

for (state_id in unique(all_markers_up$State)) {
  genes_up <- all_markers_up$gene[all_markers_up$State == state_id]
  genes_tested <- all_markers$gene[all_markers$State == state_id]
  
  foreground_ids <- unique(
    gene_mapping$ENTREZID[gene_mapping$SYMBOL %in% genes_up]
  )
  background_ids <- unique(
    gene_mapping$ENTREZID[gene_mapping$SYMBOL %in% genes_tested]
  )
  
  if (length(foreground_ids) == 0L) next
  
  go_result <- enrichGO(
    gene = foreground_ids,
    universe = background_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500,
    readable = TRUE
  )
  
  kegg_result <- enrichKEGG(
    gene = foreground_ids,
    universe = background_ids,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 5,
    maxGSSize = 500
  )
  
  go_df <- NULL
  kegg_df <- NULL
  
  if (!is.null(go_result)) {
    go_df <- as.data.frame(go_result) %>%
      mutate(group_enrich = "GO")
  }
  
  if (!is.null(kegg_result)) {
    kegg_result <- setReadable(
      kegg_result,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
    kegg_df <- as.data.frame(kegg_result) %>%
      mutate(ONTOLOGY = "KEGG", group_enrich = "KEGG")
  }
  
  state_result <- bind_rows(go_df, kegg_df)
  
  if (nrow(state_result) > 0L) {
    all_results[[state_id]] <- state_result %>%
      mutate(
        Celltype = "Astrocyte",
        State = state_id,
        Change = "Upregulated"
      )
  }
}

final_enrich <- bind_rows(all_results)
final_enrich2 <- final_enrich %>%
  filter(!is.na(p.adjust), p.adjust < 0.05)

head(final_enrich2)