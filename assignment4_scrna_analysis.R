# =============================================================================
# Assignment 4: scRNA-seq Analysis of Influenza A Virus Nasal Infection
# Dataset: Mouse nasal tissue (RM, OM, LNG) at multiple timepoints post-IAV
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 0: Install Required Packages (run once, then comment out)
# -----------------------------------------------------------------------------

# CRAN packages
install.packages(c("Seurat", "devtools", "harmony"))
install.packages("pheatmap")
install.packages('mutoss')
install.packages("metap")

# Bioconductor packages
BiocManager::install(c(
  "SingleR",           # Automatic cell type annotation
  "celldex",           # Reference datasets for SingleR
  "clusterProfiler",   # GSEA / ORA pathway analysis
  "org.Mm.eg.db",      # Mouse gene ID database (needed for clusterProfiler)
  "glmGamPoi"          # Speeds up SCTransform (optional but recommended)
))
BiocManager::install("multtest")

# GitHub packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("neurorestore/Augur")
devtools::install_github("satijalab/seurat-wrappers")

BiocManager::install("speckle") # Easier R-based differential abundance
BiocManager::install("glmGamPoi", force = TRUE)

# -----------------------------------------------------------------------------
# SECTION 1: Load Libraries
# -----------------------------------------------------------------------------

library(Seurat)
library(SeuratWrappers) 
library(harmony)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DESeq2)
library(Augur)
library(speckle)       
library(ggplot2)
library(dplyr)
library(patchwork)     
library(glmGamPoi)
library(pheatmap)
library(metap)
library(multtest)
library(mutoss)
library(cowplot)

options(future.globals.maxSize = 10 * 1024^3)
# -----------------------------------------------------------------------------
# SECTION 2: Load the Data
# -----------------------------------------------------------------------------

# Load the pre-built Seurat object provided 
# Download from: https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds
# Make sure the file is in working directory

setwd("/media/sorbaralab2/Extra_storage/Ailish/Assignment_4")

seurat_obj <- readRDS("seurat_ass4.rds")

# First, inspect the object to understand what's inside
print(seurat_obj)                          # Overview: cell count, assays, etc.
head(seurat_obj@meta.data)                 # What metadata columns exist?
table(seurat_obj$orig.ident)               # How many cells per sample?

# NOTE: Check your metadata column names carefully!
colnames(seurat_obj@meta.data)  

# -----------------------------------------------------------------------------
# SECTION 3: Quality Control (QC)
# -----------------------------------------------------------------------------
# Mitochondrial genes in mouse start with "mt-" (lowercase, unlike human "MT-")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Visualize QC metrics — look for obvious outliers
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3,
        pt.size = 0)  # pt.size=0 removes dots for cleaner plot with many cells

# Scatter plots to see relationships between metrics
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter out low-quality cells
# Adjust thresholds based on what you see in the violin plots above
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                               nFeature_RNA < 4000 &   # removes doublets (2 cells stuck together)
                               percent.mt < 20)

cat("Cells remaining after QC:", ncol(seurat_obj), "\n")

# -----------------------------------------------------------------------------
# SECTION 4: Normalization
# -----------------------------------------------------------------------------

# Regress out mitochondrial percentage as a confounding variable.
DefaultAssay(seurat_obj) <- "RNA"   # must be RNA for these functions

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = "percent.mt")


# -----------------------------------------------------------------------------
# SECTION 5: Dimensionality Reduction (PCA + UMAP)
# -----------------------------------------------------------------------------

seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)

# Elbow plot: helps decide how many PCs to use downstream
ElbowPlot(seurat_obj, ndims = 30)

# UMAP for visualization (uses top PCs)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = .3, verbose = FALSE)

# Quick look at the UMAP before integration
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Before Integration")

# IMPORTANT: Color by sample/condition to check if cells cluster by biology or by batch
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Colored by Sample — Check for Batch Effects")

# -----------------------------------------------------------------------------
# SECTION 7: Cell Type Annotation with SingleR
# -----------------------------------------------------------------------------
# SingleR automatically labels each cell by comparing its expression profile
# to a curated reference dataset of known cell types.
#
# For mouse data, we use the ImmGen reference (immune cells) and
# MouseRNAseqData (broader cell types including non-immune).

# Load reference datasets
immgen_ref <- ImmGenData()           # Immune cell reference (mouse)
mouse_ref  <- MouseRNAseqData()      # Broader mouse cell type reference

# Extract normalized expression matrix from Seurat object
# SingleR needs a matrix, not a Seurat object
expr_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

# Run SingleR annotation
# This compares each cell to reference cell types
singler_immgen <- SingleR(test = expr_matrix,
                          ref  = immgen_ref,
                          labels = immgen_ref$label.main)

singler_mouse  <- SingleR(test = expr_matrix,
                          ref  = mouse_ref,
                          labels = mouse_ref$label.main)

# Add annotations to Seurat object metadata
seurat_obj$singler_immgen <- singler_immgen$labels
seurat_obj$singler_mouse  <- singler_mouse$labels

# Visualize annotations on UMAP
DimPlot(seurat_obj, 
        group.by = "singler_immgen", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types (ImmGen Reference)")

DimPlot(seurat_obj,
        group.by = "singler_mouse", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types (MouseRNAseq Reference)")

# Also look at confidence scores — low scores mean uncertain annotation
plotScoreHeatmap(singler_immgen)
plotScoreHeatmap(singler_mouse)

marker_genes <- c("Cd3e", "Cd4", "Cd8a", #T cells
                  "Cd19", "Ms4a1", #B cells
                  "Csf1r", "Lyz2", #Macrophages
                  "S100a8", "S100a9", #neutrophils
                  "Foxj1", "Scgb1a1", #epithelial
                  "Pecam1", "Col1a1", "Ncr1" #endothelial/fibroblast/NK
                  )

gene_expr <- as.data.frame(
  t(GetAssayData(seurat_obj, assay = "RNA", slot = "data") [marker_genes, ])
)
gene_expr$barcode <- rownames(gene_expr)
write.csv(gene_expr, "marker_gene_expression.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# SECTION 8: Manual Annotation Verification 
# -----------------------------------------------------------------------------

# Set identity to cluster assignments
Idents(seurat_obj) <- "seurat_clusters"

# Find markers for all clusters
# max.cells.per.indent subsamples to 300 cells per cluster for speed
#results are nearly identical to using all cells
all_markers <- FindAllMarkers (seurat_obj,
                               only.pos = TRUE, #only upregulated markers
                               min.pct = 0.25,  #gene in at least 25% of cells
                               logfc.threshold = 0.25, #minimum fold cghange cutoff
                               test.use = "wilcox",    #fastest and well validated
                               max.cells.per.ident = 300)

# Save all markers to a CSV
write.csv(all_markers, "all_cluster_markers.csv", row.names = TRUE)

# View top 5 markers per cluster
top_5 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) 
print(top_5)

#top 10 for a heatmap
top_10 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n =10)

DoHeatmap(seurat_obj,
          features = top_10$gene,
          size = 3,
          angle = 90) +
  ggtitle("Top 10 Marker Genes Per Cluster")

top3 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 3) %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2)

print(top3, n = 100)

#now go thorugh genes and figure out what cells they belong to
#Feature plots for known mouse nasal/immune cell markers
#use them to visually confirm what each cluster is

FeaturePlot(seurat_obj,
            features = c("Ptprc",   # CD45 — all immune cells
                         "Cd3e",    # T cells
                         "Cd19",    # B cells
                         "Csf1r",   # Macrophages/monocytes
                         "Foxj1",   # Ciliated epithelial cells
                         "S100a8"), # Neutrophils
            reduction = "umap",
            ncol = 3)

# After consulting your markers and SingleR output, rename clusters
new.cluster.ids <- c(
  "0"  = "Neurons",                        # Nqo1, Cartpt
  "1"  = "Basal Epithelial Cells",         # Krt15, Serpinb5 — strong basal markers
  "2"  = "Macrophages",                    # C1qa, Ms4a7, Fcrls — very confident
  "3"  = "B Cells",                        # Iglc1, Iglc2, Fcer2a — very confident
  "4"  = "Olfactory Sensory Neurons",      # Nhlh1, Gng8 — OM specific
  "5"  = "Neutrophils",                    # Retnlg, Cxcr2, Mmp8 — very confident
  "6"  = "Endothelial Cells",              # Gpihbp1, Ptprb, Rbp7 — very confident
  "7"  = "Fibroblasts",                    # Col1a1, Lum, Omd — very confident
  "8"  = "Dendritic Cells",                # Cd209a, Ifi205 — medium confidence
  "9"  = "Olfactory Neurons",              # Kcna1, Mdga2 — likely OSN subtype
  "10" = "T Cells",                        # Cd3e, Cd3g, Themis — very confident
  "11" = "Olfactory Sensory Neurons II",      # Abca4, Nrcam — OSN subtype, diff from cluster 4
  "12" = "Monocytes",                      # Adgre4, Treml4, Ear2 — medium confidence
  "13" = "Bowman's Gland Cells",           # Chil6, Bpifb4 — OM/LNG secretory
  "14" = "Mature Olfactory Neurons",       # Adcy3, Umodl1, Cnga4 — mature OSNs
  "15" = "Specialized Epithelial",         # Smbd1, Galnt13, Kl 
  "16" = "Nasal Secretory Epithelial",     # Odam, Bpifb9a, Bpifb5 — LNG specific
  "17" = "Stromal cells",                        # Higd1b 
  "18" = "Neuroendocrine Cells",           # Tac1, Barx2 — neuronal/secretory
  "19" = "Goblet Cells",                   # Muc2, Sec14l3 — very confident
  "20" = "Immature Neutrophils",           # Elane, Mpo, Ms4a3 — neutrophil precursors
  "21" = "Ciliated Epithelial Cells",      # Dnah5, Fam183b — RM ciliated
  "22" = "NK Cells",                       # Ncr1, Klra4, Klra8 — very confident
  "23" = "Nasal Secretory Epithelial II",     # Car6, Scgb2b27 — LNG/OM secretory subtype
  "24" = "Tuft Cells",                     # Hmx3, Il25 — chemosensory tuft cells
  "25" = "Olfactory Ensheathing Cells",    # Mpz, Foxd3 — unique to OM
  "26" = "Ciliated Epithelial Cells II",      # Tmem212, Sntn — highly confident ciliated
  "27" = "Proliferating Cells",            # Pclaf, Hist1h1a — cycling, not a true cell type
  "28" = "Squamous Epithelial Cells",      # Csta1, Clca4a
  "29" = "Immature Olfactory Neurons",     # Neurod1, Neurog1 — neuronal progenitors
  "30" = "Olfactory/Taste Receptor Cells"  # Fezf2, Otop1
)

seurat_obj <- RenameIdents(seurat_obj,new.cluster.ids)
seurat_obj$cell_type <- Idents(seurat_obj)  # Save as permanent metadata column

# Final annotated UMAP
DimPlot(seurat_obj,
        group.by = "cell_type",
        label = TRUE, repel = TRUE,
        pt.size = 0.3,
        alpha = 1) +
  ggtitle("Annotated Cell Types") +
  NoLegend()
# -----------------------------------------------------------------------------
# SECTION 9: Subset to Your Comparison of Interest
# -----------------------------------------------------------------------------

# Check what values exist in your condition columns
table(seurat_obj$organ_custom)      # "RM", "OM", "LNG"
table(seurat_obj$time)   #  D02   D05   D08   D14 Naive 

# Define the correct order 
time_order <- c("Naive", "D02", "D05", "D08", "D14")

# Subset to RM tissue only
rm_obj <- subset(seurat_obj, subset = organ_custom == "RM")

rm <- subset(rm_obj, subset = time %in% c("Naive", "D02","D05", "D08", "D14"))
rm$time  <- factor(rm$time,  levels = time_order)

cat("Cells in comparison subset:", ncol(rm), "\n")
table(rm$time)

# Quick UMAP of this subset colored by time
DimPlot(rm,
        group.by = "time",
        split.by = "time") +
  ggtitle("RM")


# Subset to OM tissue only
om_obj <- subset(seurat_obj, subset = organ_custom == "OM")

om <- subset(om_obj, subset = time %in% c("Naive", "D02","D05", "D08", "D14"))
om$time  <- factor(om$time,  levels = time_order)

cat("Cells in comparison subset:", ncol(om), "\n")
table(om$time)

# Quick UMAP of this subset colored by time
DimPlot(om,
        group.by = "time",
        split.by = "time") +
  ggtitle("OM")


# Subset to LNG tissue only
lng_obj <- subset(seurat_obj, subset = organ_custom == "LNG")

lng <- subset(lng_obj, subset = time %in% c("Naive", "D02","D05", "D08", "D14"))
lng$time <- factor(lng$time, levels = time_order)

cat("Cells in comparison subset:", ncol(lng), "\n")
table(lng$time)

# Quick UMAP of this subset colored by time
DimPlot(lng,
        group.by = "time",
        split.by = "time") +
  ggtitle("LNG")


make_umap <- function(obj, title) {
  
  # Find clusters with > 100 cells in this tissue
  label_these <- names(which(table(obj$cell_type) > 100))
  
  # Annotation reference panel
  annot_panel <- DimPlot(obj,
                         group.by   = "cell_type",
                         label      = FALSE,   # turn off default labels
                         pt.size    = 0.3) +
    ggtitle(paste(title, "— Cell Types")) +
    NoLegend()
  
  # Manually add labels only for large clusters
  annot_panel <- LabelClusters(annot_panel,
                               id        = "cell_type",
                               clusters  = label_these,
                               repel     = TRUE,
                               size      = 3)
  
  # Time split panel
  time_panel <- DimPlot(obj,
                        group.by   = "cell_type",
                        split.by   = "time",
                        label      = FALSE,
                        pt.size    = 0.3) +
    ggtitle(paste(title, "— By Timepoint")) +
    NoLegend()
 
  
  plot_grid(annot_panel, time_panel, 
            nrow        = 2)
}


rmplot  <- make_umap(rm,  "Respiratory Mucosa")
omplot  <- make_umap(om,  "Olfactory Mucosa")
lngplot <- make_umap(lng, "Lateral Nasal Gland")

rmplot
omplot
lngplot

# -----------------------------------------------------------------------------
# SECTION 10: Differential Abundance 
# -----------------------------------------------------------------------------
# Differential abundance asks: does the PROPORTION of each cell type change
# Calculate proportions per mouse then average by timepoint
rmprop_df <- rm@meta.data %>%
  group_by(mouse_id, time, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mouse_id, time) %>%
  mutate(prop = n / sum(n)) %>%
  group_by(time, cell_type) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

# Stacked bar per timepoint
ggplot(rmprop_df, aes(x = time, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "Timepoint", y = "Mean Proportion",
       title = "RM Cell Type Proportions by Timepoint",
       fill  = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

omprop_df <- om@meta.data %>%
  group_by(mouse_id, time, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mouse_id, time) %>%
  mutate(prop = n / sum(n)) %>%
  group_by(time, cell_type) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

# Stacked bar per timepoint
ggplot(omprop_df, aes(x = time, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "Timepoint", y = "Mean Proportion",
       title = "OM Cell Type Proportions by Timepoint",
       fill  = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

lngprop_df <- lng@meta.data %>%
  group_by(mouse_id, time, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(mouse_id, time) %>%
  mutate(prop = n / sum(n)) %>%
  group_by(time, cell_type) %>%
  summarise(mean_prop = mean(prop), .groups = "drop")

# Stacked bar per timepoint
ggplot(lngprop_df, aes(x = time, y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "Timepoint", y = "Mean Proportion",
       title = "LNG Cell Type Proportions by Timepoint",
       fill  = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# # -----------------------------------------------------------------------------
# SECTION 11: Differential Expression — T Cells, Naive vs D08, per tissue
# -----------------------------------------------------------------------------

# Create simple mouse IDs (strips subcompartment from mouse_id)
rm$simple_mouse_id  <- gsub("_RT|_ET|_Sinus", "", rm$mouse_id)
om$simple_mouse_id  <- gsub("_RT|_ET|_Sinus", "", om$mouse_id)
lng$simple_mouse_id <- gsub("_RT|_ET|_Sinus", "", lng$mouse_id)

# Verify — should show m1, m2, m3 across timepoints cleanly
table(rm$simple_mouse_id, rm$time)

run_tcell_de <- function(tissue_obj, tissue_name) {
  
  tcell_obj <- subset(tissue_obj,
                      subset = cell_type == "T Cells" &
                        time %in% c("Naive", "D08"))
  
  cat(tissue_name, "— T cells for DE:", ncol(tcell_obj), "\n")
  print(table(tcell_obj$time, tcell_obj$simple_mouse_id))
  
  DefaultAssay(tcell_obj) <- "RNA"
  pseudo <- AggregateExpression(tcell_obj,
                                assays        = "RNA",
                                return.seurat = TRUE,
                                group.by      = c("simple_mouse_id", "time", "cell_type"))
  
  pseudo$group_label <- paste(pseudo$cell_type, pseudo$time, sep = "_")
  Idents(pseudo) <- "group_label"
  
  cat("Available groups:\n")
  print(table(pseudo$group_label))
  
  de <- FindMarkers(pseudo,
                    ident.1         = "T Cells_D08",
                    ident.2         = "T Cells_Naive",
                    test.use        = "DESeq2",
                    min.cells.group = 3)
  
  de$gene   <- rownames(de)
  de$tissue <- tissue_name
  
  write.csv(de, paste0("DE_Tcells_D08vsNaive_", tissue_name, ".csv"))
  return(de)
}

de_rm  <- run_tcell_de(rm,  "RM")
de_om  <- run_tcell_de(om,  "OM")
de_lng <- run_tcell_de(lng, "LNG")

# Quick look at top hits
cat("\n--- RM ---\n");  print(head(de_rm[order(de_rm$p_val_adj), ], 10))
cat("\n--- OM ---\n");  print(head(de_om[order(de_om$p_val_adj), ], 10))
cat("\n--- LNG ---\n"); print(head(de_lng[order(de_lng$p_val_adj), ], 10))

# -----------------------------------------------------------------------------
# SECTION 12: Volcano Plots
# -----------------------------------------------------------------------------

make_volcano <- function(de, title) {
  de$significant <- de$p_val_adj < 0.05 & abs(de$avg_log2FC) > 0.5
  
  ggplot(de, aes(x = avg_log2FC, y = -log10(p_val_adj),
                 color = significant, label = gene)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("grey70", "firebrick")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    ggrepel::geom_text_repel(data = subset(de, significant),
                             size = 3, max.overlaps = 20) +
    theme_classic() +
    ggtitle(title) +
    xlab("log2 Fold Change") + ylab("-log10(Adjusted P-value)") +
    NoLegend()
}

v_rm  <- make_volcano(de_rm,  "T Cells: D08 vs Naive — RM")
v_om  <- make_volcano(de_om,  "T Cells: D08 vs Naive — OM")
v_lng <- make_volcano(de_lng, "T Cells: D08 vs Naive — LNG")

plot_grid(v_rm, v_om, v_lng, nrow = 1)

# -----------------------------------------------------------------------------
# SECTION 13: Violin Plots — T cell subtype markers + top DE genes
# -----------------------------------------------------------------------------
# Combined condition label for cross-tissue comparison
condition_order <- c("RM_Naive", "RM_D08", "OM_Naive", "OM_D08", "LNG_Naive", "LNG_D08")

# Subset to T cells only
tcell_obj <- subset(seurat_obj, 
                    subset = cell_type == "T Cells" & 
                      time %in% c("Naive", "D08"))
tcell_obj$condition <- paste(tcell_obj$organ_custom, tcell_obj$time, sep = "_")
tcell_obj$condition <- factor(tcell_obj$condition, levels = condition_order)

condition_colors <- c(
  "RM_Naive"  = "#4E9AF1", "RM_D08"    = "#1A5EA8",
  "OM_Naive"  = "#F4A261", "OM_D08"    = "#C1440E",
  "LNG_Naive" = "#57CC99", "LNG_D08"   = "#1A7A4A"
)

# Part A: Known T cell subtype markers
tcell_subtype_markers <- c(
  "Cd4",   # CD4+ helper T cells
  "Cd8a",  # CD8+ cytotoxic T cells
  "Foxp3", # Regulatory T cells (Tregs)
  "Tcf7",  # Stem-like/memory T cells
  "Gzmb",  # Cytotoxic effector T cells
  "Ifng"   # Effector T cells
)

vln_plots_markers <- lapply(tcell_subtype_markers, function(gene) {
  VlnPlot(tcell_obj,
          features  = gene,
          group.by  = "condition",
          cols      = condition_colors,
          pt.size   = 0) +
    theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
          axis.title.x   = element_blank(),
          legend.position = "none") +
    xlab("") + ylab("Expression")
})

vln_marker_combined <- wrap_plots(vln_plots_markers, ncol = 2) +
  plot_annotation(
    title    = "T Cell Subtype Marker Expression Across Conditions",
    subtitle = "Naive vs D08 | RM, OM, LNG"
  )
print(vln_marker_combined)
ggsave("VlnPlot_Tcell_markers.pdf", vln_marker_combined, width = 12, height = 10)

# Part B: Top DE genes per tissue
top_de_rm  <- de_rm  %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  arrange(p_val_adj) %>% head(6) %>% pull(gene)
top_de_om  <- de_om  %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  arrange(p_val_adj) %>% head(6) %>% pull(gene)
top_de_lng <- de_lng %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  arrange(p_val_adj) %>% head(6) %>% pull(gene)

make_vln_panel <- function(obj, genes, title) {
  genes <- genes[genes %in% rownames(obj)]
  if (length(genes) == 0) { message("No valid genes for: ", title); return(NULL) }
  
  plots <- lapply(genes, function(gene) {
    VlnPlot(obj, features = gene, group.by = "condition",
            cols = condition_colors, pt.size = 0) +
      theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 9),
            axis.title.x   = element_blank(),
            legend.position = "none") +
      xlab("") + ylab("Expression")
  })
  
  wrap_plots(plots, ncol = 2) +
    plot_annotation(title = title,
                    subtitle = "Top upregulated genes at D08 | All conditions shown")
}

if (length(top_de_rm)  > 0) { p_rm  <- make_vln_panel(tcell_obj, top_de_rm,  "Top DE Genes D08 — RM");  print(p_rm);  ggsave("VlnPlot_TopDE_RM.pdf",  p_rm,  width = 12, height = 10) }
if (length(top_de_om)  > 0) { p_om  <- make_vln_panel(tcell_obj, top_de_om,  "Top DE Genes D08 — OM");  print(p_om);  ggsave("VlnPlot_TopDE_OM.pdf",  p_om,  width = 12, height = 10) }
if (length(top_de_lng) > 0) { p_lng <- make_vln_panel(tcell_obj, top_de_lng, "Top DE Genes D08 — LNG"); print(p_lng); ggsave("VlnPlot_TopDE_LNG.pdf", p_lng, width = 12, height = 10) }

# -----------------------------------------------------------------------------
# SECTION 14: GSEA — GO Biological Process + Hallmark, per tissue
# -----------------------------------------------------------------------------

library(msigdbr)

run_gsea <- function(de, tissue_name) {
  
  # Build ranked gene list (avg_log2FC + tiny noise to break ties)
  de_clean <- de %>%
    filter(!is.na(avg_log2FC), !is.na(p_val_adj)) %>%
    mutate(rank_score = avg_log2FC + rnorm(n(), 0, 1e-10))
  
  gene_ranks <- de_clean$rank_score
  names(gene_ranks) <- de_clean$gene
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  # Convert to Entrez IDs
  entrez_map <- bitr(names(gene_ranks),
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Mm.eg.db)
  
  gene_ranks <- gene_ranks[names(gene_ranks) %in% entrez_map$SYMBOL]
  names(gene_ranks) <- entrez_map$ENTREZID[match(names(gene_ranks), entrez_map$SYMBOL)]
  gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]
  
  cat(tissue_name, "— genes going into GSEA:", length(gene_ranks), "\n")
  
  # GO Biological Process GSEA
  gsea_go <- gseGO(geneList     = gene_ranks,
                   OrgDb        = org.Mm.eg.db,
                   ont          = "BP",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
  
  cat(tissue_name, "— GO BP terms found:", nrow(gsea_go), "\n")
  
  # Relax cutoff if nothing found
  if (nrow(gsea_go) == 0) {
    cat(tissue_name, "— relaxing to 0.2\n")
    gsea_go <- gseGO(geneList     = gene_ranks,
                     OrgDb        = org.Mm.eg.db,
                     ont          = "BP",
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = 0.2,
                     verbose      = FALSE)
    cat(tissue_name, "— GO BP terms found (relaxed):", nrow(gsea_go), "\n")
  }
  # Hallmark gene sets
  hallmark_sets <- msigdbr(species = "Mus musculus", collection = "H") %>%
    dplyr::select(gs_name, ncbi_gene) %>%
    dplyr::mutate(ncbi_gene = as.character(ncbi_gene))
  
  gsea_hallmark <- GSEA(geneList     = gene_ranks,
                        TERM2GENE    = hallmark_sets,
                        pvalueCutoff = 0.05,
                        verbose      = FALSE)
  
  cat(tissue_name, "— Hallmark sets found:", nrow(gsea_hallmark), "\n")
  
  write.csv(as.data.frame(gsea_go),
            paste0("GSEA_GO_Tcells_D08vsNaive_", tissue_name, ".csv"), row.names = FALSE)
  write.csv(as.data.frame(gsea_hallmark),
            paste0("GSEA_Hallmark_Tcells_D08vsNaive_", tissue_name, ".csv"), row.names = FALSE)
  
  return(list(go = gsea_go, hallmark = gsea_hallmark))
}

gsea_rm  <- run_gsea(de_rm,  "RM")
gsea_om  <- run_gsea(de_om,  "OM")
gsea_lng <- run_gsea(de_lng, "LNG")

# Plots
plot_gsea <- function(gsea_result, tissue_name) {
  go       <- gsea_result$go
  hallmark <- gsea_result$hallmark
  
  if (!is.null(go) && nrow(go) > 0) {
    p1 <- dotplot(go, showCategory = 15, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste("GO BP GSEA: T Cells D08 vs Naive —", tissue_name))
    print(p1)
    ggsave(paste0("GSEA_GO_dotplot_", tissue_name, ".pdf"), p1, width = 14, height = 8)
  }
  
  if (!is.null(hallmark) && nrow(hallmark) > 0) {
    p2 <- dotplot(hallmark, showCategory = 15, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(paste("Hallmark GSEA: T Cells D08 vs Naive —", tissue_name))
    print(p2)
    ggsave(paste0("GSEA_Hallmark_dotplot_", tissue_name, ".pdf"), p2, width = 14, height = 8)
  }
}

plot_gsea(gsea_rm,  "RM")
plot_gsea(gsea_om,  "OM")
plot_gsea(gsea_lng, "LNG")

# -----------------------------------------------------------------------------
# SECTION 15: Save processed object
# -----------------------------------------------------------------------------

saveRDS(seurat_obj, "seurat_ass4_processed.rds")