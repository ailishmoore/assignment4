# Assignment 4: scRNA-seq Analysis of Influenza A Virus Nasal Infection
Mouse nasal tissue scRNA-seq — [PMC11324402](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324402/)

---

## Table of Contents

1. [Repository Structure](#repository-structure)
2. [Introduction](#introduction)
3. [Methods](#methods)
4. [Results](#results)
5. [Discussion](#discussion)
6. [References](#references)

---

## Repository Structure

```
Assignment_4/
├── README.md                                  # This file
├── assignment4_scrna_analysis.R               # Full annotated R analysis script
│
├── figures/
│   ├── UMAP_annotated_celltypes.png          # Fig 1  — Annotated UMAP of all 31 clusters
│   ├── UMAP_RM_by_timepoint.png              # Fig 2  — RM UMAP by cell type and timepoint
│   ├── UMAP_OM_by_timepoint.png              # Fig 3  — OM UMAP by cell type and timepoint
│   ├── UMAP_LNG_by_timepoint.png             # Fig 4  — LNG UMAP by cell type and timepoint
│   ├── CellProportion_RM.png                 # Fig 5  — Cell type proportions by timepoint, RM
│   ├── CellProportion_OM.png                 # Fig 6  — Cell type proportions by timepoint, OM
│   ├── CellProportion_LNG.png                # Fig 7  — Cell type proportions by timepoint, LNG
│   ├── T_cells_D08_vs_Naive.png              # Fig 8  — Volcano plots: T cells D08 vs Naive, all tissues
│   ├── VlnPlot_TopDE_RM.png                  # Fig 9  — Violin plots: top upregulated DE genes, RM
│   ├── VlnPlot_TopDE_LNG.png                 # Fig 10 — Violin plots: top upregulated DE genes, LNG
│   ├── VlnPlot_Tcell_markers.png             # Fig 11 — Violin plots:  T cell subtype markers across conditions  
│   ├── GSEA_GO_dotplot_RM.png                # Fig 12 — GO BP GSEA dotplot, RM
│   ├── GSEA_Hallmark_dotplot_RM.png          # Fig 13 — Hallmark GSEA dotplot, RM
│   ├── GSEA_GO_dotplot_OM.png                # Fig 14 — GO BP GSEA dotplot, OM
│   ├── GSEA_Hallmark_dotplot_OM.png          # Fig 15 — Hallmark GSEA dotplot, OM 
│   ├── GSEA_GO_dotplot_LNG.png               # Fig 16 — GO BP GSEA dotplot, LNG
│   └── GSEA_Hallmark_dotplot_LNG.png         # Fig 17 — Hallmark GSEA dotplot, LNG
│
├── output_tables/
│   ├── all_cluster_markers.csv                 # Marker genes for all clusters
│   ├── DE_Tcells_D08vsNaive_RM.csv             # DE results: T cells, RM
│   ├── DE_Tcells_D08vsNaive_OM.csv             # DE results: T cells, OM
│   ├── DE_Tcells_D08vsNaive_LNG.csv            # DE results: T cells, LNG
│   ├── GSEA_GO_Tcells_D08vsNaive_RM.csv        # GO GSEA results — RM
│   ├── GSEA_GO_Tcells_D08vsNaive_OM.csv        # GO GSEA results — OM
│   ├── GSEA_GO_Tcells_D08vsNaive_LNG.csv       # GO GSEA results — LNG
│   ├── GSEA_Hallmark_Tcells_D08vsNaive_RM.csv  # Hallmark GSEA results - RM
│   ├── GSEA_Hallmark_Tcells_D08vsNaive_OM.csv  # Hallmark GSEA results - OM
│   └── GSEA_Hallmark_Tcells_D08vsNaive_LNG.csv # Hallmark GSEA results - LNG
│
└── seurat_ass4_processed.rds                  # Final processed Seurat object (not tracked by Git — >100 MB)
```

> **Note:** Raw data (`seurat_ass4.rds`, `seurat_ass4_processed.rds`) are not tracked by Git due to file size. The raw object can be downloaded from [this link](https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds).

---

## Introduction

Influenza A Virus (IAV) infection of the upper respiratory tract triggers a dynamic immune response that varies both temporally and spatially across nasal tissue compartments. The nasal cavity is composed of functionally distinct regions: the **respiratory mucosa (RM)**, which lines the nasal turbinates and septum and serves as the primary site of viral entry; the **olfactory mucosa (OM)**, a specialized neuroepithelial tissue involved in chemosensation that is particularly vulnerable to infection-associated damage; and the **lateral nasal gland (LNG)**, a secretory gland contributing to mucosal defense (Zazhytska et al., 2024). These compartments differ substantially in cell type composition, barrier function, and immune cell residency, making them likely to mount distinct responses to infection.

It has been demonstrated that IAV infection of the nasal mucosa drives myeloid cell recruitment by day 5 post-infection (dpi), followed by T cell infiltration peaking around 8 dpi and viral clearance by 14 dpi (Zazhytska et al., 2024). However, the tissue-compartment-specific transcriptional dynamics of these immune populations — particularly T cells — remain incompletely characterized at single-cell resolution. This analysis uses scRNA-seq to profile cell type composition and gene expression across all three nasal compartments at five timepoints (naive, D02, D05, D08, D14), focusing on T cell responses during active infection (naive vs. D08).

Seurat (v5) was used as the primary framework for data processing, dimensionality reduction, clustering, and visualization. Seurat is the most widely adopted toolkit for single-cell analysis in R, with extensive documentation and community support (Hao et al., 2021). Its integrated workflow, from QC through UMAP, enables reproducible, modular analysis. Compared to alternatives such as Scanpy (Python-based; Wolf et al., 2018), Seurat offers tighter integration with downstream R-based tools including DESeq2 and clusterProfiler, making it the natural choice for an R-centric pipeline.

Normalization was performed using the standard log-normalization approach (`NormalizeData`, `FindVariableFeatures`, `ScaleData`), with mitochondrial percentage regressed out during scaling to reduce technical confounding. SCTransform uses a regularized negative binomial model to stabilize variance across the dynamic range of expression, and may outperform log-normalization for datasets with large cell count variation (Hafemeister & Satija, 2019). Log-normalization was selected here for interpretability and speed given the dataset size.

Dimensionality reduction relied on PCA followed by Uniform Manifold Approximation and Projection (UMAP). UMAP preserves both local and some global structure in high-dimensional data and is the current standard for scRNA-seq visualization, having largely replaced t-SNE due to better scalability and more meaningful inter-cluster distances (McInnes et al., 2018; Becht et al., 2019). The number of PCs used (30) was selected based on the elbow plot.

Cell type annotation combined automatic labeling via SingleR with manual verification using canonical marker genes. SingleR assigns cell identities by correlating each cell's expression profile to curated reference datasets — here, ImmGen (immune cells) and MouseRNAseqData (broader cell types) from the `celldex` package (Aran et al., 2019). This dual-reference approach is beneficial because ImmGen provides fine-grained immune cell resolution while MouseRNAseqData captures non-immune populations such as epithelial and neuronal cells. SingleR is faster and less subjective than purely manual annotation, though it is limited to cell types represented in its reference and can misclassify rare or tissue-specific populations. Manual verification using `FeaturePlot` and `DoHeatmap` was therefore used to confirm and adjust assignments for tissue-specific clusters such as olfactory sensory neurons and Bowman's gland cells not well-represented in standard references. Additionally, PanglaoDB was used to further validate cell type annotations based on curated marker gene sets (Franzén et al., 2019; Clarke et al., 2021).

Differential expression analysis was performed using a pseudobulk approach via `AggregateExpression` followed by DESeq2 (Love et al., 2014). Pseudobulk DE aggregates counts per sample before testing, which properly models biological replication and avoids the inflated false positive rates that result from treating individual cells as independent observations — a known limitation of single-cell-level tests such as Wilcoxon rank-sum (Squair et al., 2021). DESeq2 uses a negative binomial generalized linear model and has well-established performance for bulk and pseudobulk RNA-seq. The comparison was restricted to T cells (cluster 10) between naive and D08 timepoints, per tissue.

Gene Set Enrichment Analysis (GSEA) was carried out using clusterProfiler (Wu et al., 2021) with two gene set collections: GO Biological Process terms (via `gseGO`) and MSigDB Hallmark gene sets (via `GSEA` with `msigdbr`). GSEA ranks all tested genes by log2 fold change and tests whether gene sets are enriched at the top or bottom of this ranked list, making it more sensitive than over-representation analysis (ORA) for detecting coordinated pathway-level changes when no single gene clears a strict significance threshold. Hallmark gene sets represent well-defined biological states with reduced redundancy compared to GO, making them complementary for biological interpretation. clusterProfiler was chosen over alternatives like fgsea or GAGE for its rich visualization functions and tight integration with Bioconductor annotation databases.

Cell type proportion analysis was performed by calculating mean per-mouse proportions across timepoints and visualizing changes as stacked bar charts. This approach provides an intuitive summary of compositional dynamics across the infection time course. For formal statistical testing of differential abundance, tools such as `speckle` (propeller) or Milo (which uses k-nearest-neighbor graphs) could be applied in follow-up analysis; the descriptive proportion plots here provide sufficient context for interpreting DE and GSEA results.

---

## Methods

### Data and Quality Control

A pre-processed Seurat object containing scRNA-seq data from mouse nasal tissue (RM, OM, LNG) across five timepoints (naive, D02, D05, D08, D14; n = 3 mice per group) was provided. Data originate from Zazhytska et al. (2024; PMC11324402). Quality control filtering retained cells with 200–4,000 detected features and < 20% mitochondrial reads, removing likely empty droplets, doublets, and dying cells. Mitochondrial content was calculated using `PercentageFeatureSet` targeting genes matching the `^mt-` pattern.

### Normalization and Feature Selection

Normalized gene expression was computed using log-normalization (scale factor 10,000) via `NormalizeData`. The top 2,000 highly variable features were identified with `FindVariableFeatures`. Data were scaled using `ScaleData` with regression of mitochondrial percentage (`vars.to.regress = "percent.mt"`).

### Dimensionality Reduction and Clustering

PCA was run on the top 2,000 variable features (30 PCs). UMAP was computed on the top 30 PCs (`RunUMAP`, `dims = 1:30`). Shared nearest-neighbor graphs were built with `FindNeighbors` (dims = 1:30) and clusters identified using the Louvain algorithm at resolution 0.3 (`FindClusters`), yielding 31 initial clusters (0–30).

### Cell Type Annotation

Automatic annotation was performed with SingleR v2 using ImmGen and MouseRNAseqData references (celldex v1). Manual annotation was subsequently applied using `FindAllMarkers` (Wilcoxon test, `min.pct = 0.25`, `logfc.threshold = 0.25`, `only.pos = TRUE`) and cross-referenced with canonical marker genes visualized via `FeaturePlot` and `DoHeatmap`. PanglaoDB was used to further validate cell type annotations. Final cell type labels were assigned using `RenameIdents`, identifying 31 populations including immune, epithelial, stromal, and neuronal cell types.

### Tissue Subsetting and Composition Analysis

Cells were subset by tissue (RM, OM, LNG) and timepoints confirmed. Cell type proportions were calculated per mouse, then averaged across mice within each timepoint, and visualized as stacked bar charts using ggplot2.

### Differential Expression

T cells (cluster 10) were subset from each tissue and compared between naive and D08 timepoints. Pseudobulk counts were generated with `AggregateExpression` (grouped by mouse ID, timepoint, and cell type). DESeq2 was used for DE testing via `FindMarkers` (`test.use = "DESeq2"`, `min.cells.group = 3`). Results were visualized as volcano plots with significance thresholds of adjusted p < 0.05 and |log2FC| > 0.5.

### Gene Set Enrichment Analysis

Ranked gene lists (ordered by log2FC) were generated from each tissue's T cell DE results and mapped to Entrez IDs via `bitr` (org.Mm.eg.db). GSEA was run using `gseGO` (GO Biological Process, `pvalueCutoff = 0.05`, minGSSize = 10, maxGSSize = 500) and `GSEA` with MSigDB Hallmark gene sets (Mus musculus) from `msigdbr`. Results were visualized as dotplots split by enrichment direction.

### Software

All analyses were performed in R (v4.3). Key packages: Seurat v5, SingleR v2, celldex v1, clusterProfiler v4, DESeq2 v1.42, org.Mm.eg.db v3.18, msigdbr v7.5, ggplot2 v3.5, patchwork v1.2, cowplot v1.1.

Full code is provided in `assignment4_scrna_analysis.R`.

---

## Results

### Quality Control and Cell Recovery

After applying QC thresholds (200–4,000 detected features per cell and less than 20% mitochondrial reads), a high-quality set of cells was retained for downstream analysis. Violin plots of nFeature_RNA, nCount_RNA, and percent.mt confirmed removal of low-complexity and stressed cells without excessive data loss.

### UMAP and Cell Type Identification

Clustering at resolution 0.3 identified 31 distinct clusters, which were annotated using a combination of SingleR automated labeling and manual marker-gene verification (Figure 1). The full dataset captures remarkable cellular diversity spanning immune, epithelial, neuronal, and stromal compartments. Immune populations include T cells, B cells, macrophages, dendritic cells, monocytes, neutrophils, immature neutrophils, and NK cells. Epithelial populations encompass ciliated epithelial cells (two subtypes), goblet cells, basal epithelial cells, squamous epithelial cells, nasal secretory epithelial cells (two subtypes), specialized epithelial cells, and tuft cells. Multiple neuronal populations were identified — olfactory sensory neurons (two subtypes), mature olfactory neurons, immature olfactory neurons, olfactory ensheathing cells, olfactory/taste receptor cells, and a broader neuron cluster — reflecting the substantial olfactory tissue contribution from the OM. Stromal populations included fibroblasts, endothelial cells, Bowman's gland cells, neuroendocrine cells, and a stromal cell cluster (Cluster 17, whose identity warrants further investigation). A proliferating cell cluster was also identified, likely capturing cycling cells across multiple lineages rather than a distinct cell type.
![UMAP](figures/UMAP_annotated_celltypes.png)
**Figure 1.** UMAP of all cells colored by annotated cell type. Clustering at resolution 0.3 identified 31 populations spanning immune, epithelial, neuronal, and stromal lineages. Labels are placed at cluster centroids; Cluster 17 remains provisionally unannotated.

Tissue-specific populations were clearly apparent when the object was subset by compartment (Figures 2–4). The OM was enriched for olfactory sensory neurons and olfactory ensheathing cells, while the LNG was dominated at baseline by secretory epithelial populations including Bowman's gland cells and nasal secretory epithelial cells. The RM contained the most diverse immune infiltrate. These tissue-specific UMAP layouts are consistent with the known biology of each compartment.
![UMAP](figures/UMAP_RM_by_timepoint.png)
**Figure 2.** UMAP of respiratory mucosa (RM) cells. Top panel: cells colored by annotated cell type. Bottom panel: same UMAP split by timepoint (naive, D02, D05, D08, D14).
![UMAP](figures/UMAP_OM_by_timepoint.png)
**Figure 3.** UMAP of olfactory mucosa (OM) cells. Top panel: cells colored by annotated cell type, highlighting the enrichment of olfactory sensory neuron populations. Bottom panel: UMAP split by timepoint (naive, D02, D05, D08, D14).
![UMAP](figures/UMAP_LNG_by_timepoint.png)
**Figure 4.** UMAP of lateral nasal gland (LNG) cells. Top panel: cells colored by annotated cell type, showing the secretory epithelial dominance of this compartment. Bottom panel: UMAP split by timepoint (naive, D02, D05, D08, D14).


### Cell Type Composition Across Timepoints

Stacked bar charts of mean per-mouse cell type proportions at each timepoint revealed dynamic compositional shifts following IAV infection across all three tissues (Figures 5–7). In the RM, neutrophil and immature neutrophil proportions expanded substantially at D05, consistent with the peak myeloid recruitment phase described in the source paper, and then decreased by D08–D14. Most strikingly, and consistent with the T cell infiltration phase annotated in the experimental design (D08), T cell proportions increased visibly at D08 in the RM relative to naive and earlier timepoints. A similar pattern was observed in the OM and LNG, where T cells also rose at D08, though from lower baseline proportions. The LNG showed the highest epithelial cell dominance at naive, with immune proportions remaining relatively modest even at peak infiltration. Given this clear increase in T cell abundance at the D08 timepoint — the timepoint designated as "T cell infiltration" in the source paper — T cells were selected as the focus for differential expression analysis comparing naive and D08.
![Cell Proportions](figures/CellProportion_RM.png)
**Figure 5.** Stacked bar chart of mean cell type proportions per timepoint in the RM. Each bar represents the average proportion across n=3 mice. The expansion of neutrophils at D05 and T cells at D08 are consistent with the myeloid recruitment and T cell infiltration phases of infection respectively.
![Cell Proportions](figures/CellProportion_OM.png)
**Figure 6.** Stacked bar chart of mean cell type proportions per timepoint in the OM. T cell proportions rise at D08 from a lower baseline than in the RM, consistent with secondary immune infiltration of this neuroepithelial compartment.
![Cell Proportions](figures/CellProportion_LNG.png)
**Figure 7.** Stacked bar chart of mean cell type proportions per timepoint in the LNG. Secretory epithelial populations dominate at baseline; immune proportions remain modest relative to RM and OM even at peak infiltration.

### Differential Expression in T Cells: Naive vs. D08

Pseudobulk DESeq2 analysis comparing D08 to naive T cells revealed strongly tissue-specific patterns of transcriptional response (Figure 8; full results in output_tables/). The RM showed the most extensive differential expression, with a large number of genes significantly upregulated at D08 (adjusted p < 0.05, |log2FC| > 0.5), including Ly6c2, Nkg7, Lgals1, Cd3d, Ms4a6b, and the interferon-stimulated gene Gm42418. The OM yielded no significantly upregulated genes at these thresholds, with the few significant hits being downregulated ribosomal genes. The LNG produced only a single significant upregulated gene, Gm42418. These tissue-specific differences likely reflect both the higher direct viral burden in the RM and the lower T cell numbers recovered from OM and LNG, which reduce statistical power under the pseudobulk framework.
![Volcano Plot](figures/T_cells_D08_vs_Naive.png)
**Figure 8.** Volcano plots of pseudobulk DESeq2 differential expression results for T cells at D08 vs naive, shown for all three tissues (RM, OM, LNG). Significant genes (adjusted p < 0.05, |log2FC| > 0.5) are highlighted in red. The RM displays the largest number of significantly upregulated genes; no significant upregulation was detected in the OM and only Gm42418 reached significance in the LNG.

Violin plots of the top RM upregulated genes showed that while Gm42418, Cd3d, Nkg7, Lgals1, and Ms4a6b are broadly expressed across tissues and conditions, Ly6c2 stands out as more selectively elevated in RM_D08 and to a lesser degree OM_Naive, consistent with recruitment or expansion of Ly6c2-expressing effector T cells specifically in the respiratory mucosa at day 8 (Figure 9). For the LNG, the strong upregulation of Gm42418 at D08 relative to the notably low naive baseline represents the clearest Naive-to-D08 transcriptional shift in that compartment (Figure 10). OM did not have any top upregulated T cell genes. 
![Violin Plot](figures/VlnPlot_TopDE_RM.png)
**Figure 9.** Violin plots of the top significantly upregulated genes in T cells at D08 vs naive in the RM (Gm42418, Ly6c2, Cd3d, Nkg7, Lgals1, Ms4a6b), shown across all six tissue-timepoint conditions. Expression is highest at RM D08, confirming tissue-specific induction.
![Violin Plot](figures/VlnPlot_TopDE_LNG.png)
**Figure 10.** Violin plots of the top significantly upregulated gene in T cells at D08 vs naive in the LNG (Gm42418), shown across all six tissue-timepoint conditions.

Violin plots of T cell subtype markers confirmed Cd8a expression within the cluster, with the largest violin bodies in RM_D08 and OM_D08, consistent with enrichment of cytotoxic T cells at these sites following infection (Figure 11). Tcf7, a marker of stem-like and memory T cells, was broadly expressed across RM_D08, OM_Naive, OM_D08, and LNG_Naive, suggesting a substantial memory or progenitor component throughout tissues. Gzmb was markedly and specifically elevated in RM_D08, reinforcing the site-specific induction of cytotoxic effector activity in the respiratory mucosa. Ifng remained low across all conditions, and both Foxp3 and Cd4 were minimally detected, indicating that regulatory and CD4+ helper T cell populations are not prominently represented in this cluster at this timepoint.
![Violin Plot](figures/VlnPlot_Tcell_markers.png)
**Figure 11.** Violin plots of T cell subtype marker expression (Cd4, Cd8a, Foxp3, Tcf7, Gzmb, Ifng) across six tissue-timepoint conditions (RM, OM, and LNG at naive and D08). Gzmb and Ifng are elevated at D08 in the RM, with more modest increases in the OM and LNG.

### Gene Set Enrichment Analysis

GSEA using GO Biological Process terms and MSigDB Hallmark gene sets revealed strongly tissue-dependent pathway enrichment profiles (Figures 12–17). In the RM, the most significantly enriched activated terms included defense response to virus, response to virus, negative regulation of viral genome replication, chromosome segregation, and nuclear division, consistent with a broad antiviral transcriptional program in infiltrating T cells (Figure 12). Hallmark GSEA in the RM identified enrichment of E2F targets, Interferon Alpha Response, DNA Repair, MYC targets, G2M Checkpoint, Interferon Gamma Response, and mTORC1 signaling as activated gene sets, alongside suppression of metabolic and stromal programs (Figure 13). These results together point to actively proliferating, interferon-stimulated effector T cells in the RM at D08.
![GSEA GO](figures/GSEA_GO_dotplot_RM.png)
**Figure 12.** GO Biological Process GSEA dotplot for T cells D08 vs naive in the RM. Dot size reflects gene count; color reflects adjusted p-value. Activated terms (left panel) include defense response to virus, response to interferon-β, cell killing, chromosome segregation, and cytoplasmic translation. Suppressed terms (right panel) include synaptic translation and neuronal signaling pathways.
![GSEA Hallmark](figures/GSEA_Hallmark_dotplot_RM.png)
**Figure 13.** MSigDB Hallmark GSEA dotplot for T cells D08 vs naive in the RM. Activated sets (left) include E2F Targets, Interferon Alpha Response, DNA Repair, MYC Targets, G2M Checkpoint, Interferon Gamma Response, and mTORC1 Signaling. Suppressed sets (right) include metabolic and stromal programs (Oxidative Phosphorylation, TNFα Signaling, Adipogenesis, Epithelial Mesenchymal Transition).

In the OM, GO BP GSEA identified only suppressed terms, synaptic translation pathways (translation at presynapse, postsynapse, and synapse) and sensory perception of chemical stimulus, with no significantly activated immune terms (Figure 14). This suppression of neuronal gene sets in T cells within the OM either reflects contamination of the T cell gene expression profile with olfactory neuron-associated transcripts due to the tissue context, or a genuine downregulation of these pathways in the local T cell population. Hallmark GSEA in the OM did recover activated immune terms including Interferon Gamma Response, Interferon Alpha Response, E2F Targets, and G2M Checkpoint, alongside suppression of P53 pathway, KRAS signaling, and TNFα signaling via NFκB (Figure 15).
![GSEA GO](figures/GSEA_GO_dotplot_OM.png)
**Figure 14.** GO Biological Process GSEA dotplot for T cells D08 vs naive in the OM. Only suppressed terms were detected, including synaptic translation (pre-, post-, and synapse), sensory perception of chemical stimulus, and regulation of fat cell differentiation, likely reflecting the neuronal tissue context of the OM.
![GSEA Hallmark](figures/GSEA_Hallmark_dotplot_OM.png)
**Figure 15.** MSigDB Hallmark GSEA dotplot for T cells D08 vs naive in the OM. Activated sets include Interferon Gamma Response, Interferon Alpha Response, E2F Targets, and G2M Checkpoint. Suppressed sets include P53 Pathway, KRAS Signaling Up, Estrogen Response Early, and TNFα Signaling via NFκB.

In the LNG, GO BP GSEA revealed exclusively activated terms centered on cell division and mitotic processes — including cell division, spindle organization, mitotic sister chromatid segregation, and negative regulation of mitotic nuclear division. This suggests the dominant transcriptional signature of D08 T cells in the LNG is one of active proliferation (Figure 16). Hallmark GSEA in the LNG was sparse, identifying only G2M Checkpoint (activated, p.adjust ~0.04) and Oxidative Phosphorylation (weakly activated), consistent with a proliferating but metabolically engaged T cell population (Figure 17).
![GSEA GO](figures/GSEA_GO_dotplot_LNG.png)
**Figure 16.** GO Biological Process GSEA dotplot for T cells D08 vs naive in the LNG. All enriched terms are in the activated direction and are centered on mitotic and cell division processes, including cell division, spindle organization, and negative regulation of mitotic nuclear division, consistent with a proliferating T cell population.
![GSEA Hallmark](figures/GSEA_Hallmark_dotplot_LNG.png)
**Figure 17.** MSigDB Hallmark GSEA dotplot for T cells D08 vs naive in the LNG. Only two sets reached significance: G2M Checkpoint (activated, p.adjust ~0.04) and Oxidative Phosphorylation (weakly activated), consistent with a proliferating, metabolically active T cell state.

---

## Discussion

This analysis demonstrates that IAV infection of the mouse nasal mucosa drives tissue-compartment-specific changes in both immune cell composition and T cell transcriptional state, with the D08 timepoint capturing the peak of T cell infiltration across all three tissue sites.

The proportional increase in T cells at D08 — observed in the RM, OM, and LNG — aligns directly with the experimental design of the source paper, which characterizes this timepoint as the "T cell infiltration" phase of the nasal immune response (Zazhytska et al., 2024). This temporal pattern is well-established for IAV infection more broadly: following viral replication and peak myeloid recruitment, antigen-specific CD8+ and CD4+ T cells are recruited to sites of infection to mediate viral clearance (Hogan et al., 2001). The presence of both Cd4 and Cd8a expressing cells in the D08 T cell cluster confirms that both lineages are represented, though the high expression of cytotoxic markers Gzmb and Ifng at this timepoint suggests cytotoxic effector T cells are the dominant activated population.

The RM showed the most transcriptionally active T cell response at D08, with a large number of DE genes and the richest GSEA enrichment profile. The activation of antiviral GO terms (defense response to virus, negative regulation of viral genome replication) and Hallmark interferon response gene sets in RM T cells is consistent with IFN-γ-driven effector function. IFN-γ produced by activated CD8+ T cells plays a well-characterized role in amplifying innate antiviral responses and facilitating viral clearance (Schoenborn & Wilson, 2007). The concurrent enrichment of E2F targets, MYC targets, and G2M Checkpoint gene sets indicates that these T cells are also actively proliferating — a key feature of clonal expansion during an adaptive immune response (Kaech & Wherry, 2007).

The olfactory mucosa presents a more complex picture. The absence of activated immune GO terms in the OM GSEA, contrasted with the appearance of suppressed neuronal translation pathways, warrants careful interpretation. The OM is heavily enriched for olfactory sensory neurons (OSNs), which are transcriptionally distinctive non-immune cells. The suppression of synaptic translation and sensory perception gene sets in OM T cells at D08 may reflect a shift away from a neuronally-influenced baseline transcriptional state, rather than active suppression of these pathways within T cells per se. Nonetheless, Hallmark GSEA did recover interferon response and proliferation signatures in OM T cells, confirming that antiviral T cell activation occurs here as well. Prior work has established that IAV infection can cause OSN loss and olfactory dysfunction, and it has been proposed that immune-mediated bystander damage may contribute alongside direct viral cytopathology (Zazhytska et al., 2022). The infiltration of cytotoxic T cells into the OM at D08, as observed here, is consistent with this hypothesis.

The LNG T cell GSEA profile, dominated by proliferation-associated GO terms (cell division, spindle organization, mitotic regulation) and the G2M Checkpoint Hallmark, suggests that T cells in this compartment at D08 are primarily in an active proliferative state rather than fully differentiated effectors. This may reflect that the LNG, as a glandular tissue, is a site of secondary T cell expansion following initial activation, or that T cell recruitment to this compartment is slightly delayed relative to the RM. The relatively sparse number of DE genes in the LNG compared to the RM is consistent with a T cell population that is expanding but not yet fully transcriptionally activated.

Taken together, these results paint a coherent picture of a spatially graded T cell response to nasal IAV infection: the RM mounts the strongest and most transcriptionally mature effector response, the OM shows evidence of T cell activation with signatures of neuronal tissue context, and the LNG shows a predominantly proliferative signature consistent with active clonal expansion. Limitations of this analysis include the descriptive rather than statistically formal treatment of cell type composition changes, the use of a uniform clustering resolution across all tissues, and the small number of biological replicates (n=3 per group) which limits the power of the pseudobulk DE approach — particularly for the OM and LNG where T cell numbers are lower. Future work could apply Milo (Dann et al., 2022) for rigorous differential abundance testing, or trajectory analysis tools such as Monocle3 to map T cell differentiation states across the infection time course.

---

## References

- Aran D et al. (2019). Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*, 20, 163–172. https://doi.org/10.1038/s41590-018-0276-y
- Becht E et al. (2019). Dimensionality reduction for visualizing single-cell data using UMAP. *Nature Biotechnology*, 37, 38–44. https://doi.org/10.1038/nbt.4314
- Clarke ZA et al. (2021). Tutorial: guidelines for annotating single-cell transcriptomic maps using automated and manual methods. *Nature Protocols* 16, 2749–2764. https://doi.org/10.1038/s41596-021-00534-0
- Dann E et al. (2022). Differential abundance testing on single-cell data using k-nearest neighbor graphs. *Nature Biotechnology*, 40, 245–253. https://doi.org/10.1038/s41587-021-01033-z
- Franzén O, Gan L-M & Björkegren JLM. (2019). PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. *Database* baz046. https://doi.org/10.1093/database/baz046
- Hafemeister C & Satija R (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. *Genome Biology*, 20, 296. https://doi.org/10.1186/s13059-019-1874-1
- Hao Y et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*, 184(13), 3573–3587. https://doi.org/10.1016/j.cell.2021.04.048
- Hogan RJ et al. (2001). Activated antigen-specific CD8+ T cells persist in the lungs following recovery from respiratory virus infections. *Journal of Immunology*, 166(3), 1813–1822. https://doi.org/10.4049/jimmunol.166.3.1813
- Kaech SM & Wherry EJ (2007). Heterogeneity and cell-fate decisions in effector and memory CD8+ T cell differentiation during viral infection. *Immunity*, 27(3), 393–405. https://doi.org/10.1016/j.immuni.2007.08.007
- Love MI, Huber W & Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550. https://doi.org/10.1186/s13059-014-0550-8
- McInnes L, Healy J & Melville J (2018). UMAP: Uniform Manifold Approximation and Projection for dimension reduction. *arXiv:1802.03426*. https://doi.org/10.48550/arXiv.1802.03426
- Schoenborn JR & Wilson CB (2007). Regulation of interferon-γ during innate and adaptive immune responses. *Advances in Immunology*, 96, 41–101. https://doi.org/10.1016/S0065-2776(07)96002-2
- Squair JW et al. (2021). Confronting false discoveries in single-cell differential expression. *Nature Communications*, 12, 5692. https://doi.org/10.1038/s41467-021-25960-2
- Wolf FA, Angerer P & Theis FJ (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15. https://doi.org/10.1186/s13059-017-1382-0
- Wu T et al. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *Innovation*, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141
- Zazhytska M et al. (2024). Primary nasal influenza infection rewires tissue-scale memory response dynamics. *Immunity*. 57(8),1955–1974.e8. https://doi.org/10.1016/j.immuni.2024.06.005

