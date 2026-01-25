library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratWrappers)
devtools::load_all("D:/BCRdataset/PRJNA890633/monocle_2.30.0/monocle")
library(patchwork)
library(umap)
library(anndata)
library(dplyr)
library(tidyr)
library(magrittr)
library(dsb)
library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
library(patchwork)

create_volcano_plot <- function(markers_df, cluster_id) {
  cluster_data <- markers_df[markers_df$cluster == cluster_id, ]
  
  # Add significance categories
  cluster_data$significance <- ifelse(cluster_data$p_val_adj < 0.05 & abs(cluster_data$avg_log2FC) > 0.25, 
                                      "Significant", "Not Significant")
  
  # Select top genes by both p-value and fold change
  sig_genes <- cluster_data[cluster_data$significance == "Significant", ]
  
  # Get top upregulated and downregulated genes
  top_up <- head(sig_genes[sig_genes$avg_log2FC > 0, ][order(sig_genes[sig_genes$avg_log2FC > 0, ]$p_val_adj), ], 5)
  top_down <- head(sig_genes[sig_genes$avg_log2FC < 0, ][order(sig_genes[sig_genes$avg_log2FC < 0, ]$p_val_adj), ], 5)
  top_genes <- rbind(top_up, top_down)
  
  ggplot(cluster_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = significance), alpha = 0.6) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    geom_text_repel(data = top_genes,
                    aes(label = gene), 
                    size = 3, 
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    force = 2,           # Increase repulsion force
                    min.segment.length = 0) +
    theme_minimal() +
    labs(title = paste("Volcano Plot - Cluster", cluster_id),
         x = "Average Log2 Fold Change",
         y = "-log10(Adjusted P-value)") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)
}



datapath <- "D:/BCRdataset/CNS_GSE203552/samples/p1_biopsy"
embpath <- "D:/unifiedBCR/benchmark/Embeddings/CNS_GSE203552/p1_biopsy/mean_inner01.pth"
outdir <- "D:/BCRdataset/CNS_GSE203552"

bcrfile <- read.table(file=file.path(datapath, "processed_bcr.csv"), header = TRUE, row.names="barcode", sep=",")
annotation <- read.table(file=file.path(datapath, "annotation.csv"), sep=",", header=TRUE, row.names="barcode")
bcremb <- read.table(file=file.path(embpath, "embedding.csv"), header = TRUE, sep=",")
gexemb <- read.table(file=file.path(embpath, "gexembedding.csv"), header = TRUE, sep=",")

immune.combined <- readRDS(file=file.path(outdir, "rds_save/p1_biopsy.rds"))
immune.combined <- AddMetaData(immune.combined, annotation)
indices <- which(immune.combined@meta.data$celltypes == "")
immune.combined@meta.data[indices, "celltypes"] <- "NA"


immune.combined@meta.data$cdr3 <- bcrfile$cdr3
# Concatenate the embeddings along axis 1
embeddings <- cbind(bcremb, gexemb)
row_norms <- sqrt(rowSums(embeddings^2))
normalized_matrix <- embeddings / row_norms
rownames(normalized_matrix) <- rownames(bcrfile)
inte_assay <- CreateAssayObject(counts = t(normalized_matrix))
immune.combined[["inteemb"]] <- inte_assay
# keep leiden is not NA
keep_cells <- !is.na(immune.combined@meta.data$leiden)
immune.combined <- immune.combined[, keep_cells]

DefaultAssay(immune.combined) <- "RNA"
immune.combined@meta.data$leiden <- as.character(immune.combined@meta.data$leiden )
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE, features=rownames(immune.combined))
immune.combined <- FindNeighbors(immune.combined, dims = 1:30, annoy.metric = "cosine")
immune.combined <- FindClusters(immune.combined, resolution = 0.5, annoy.metric = "cosine")
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
DimPlot(immune.combined, reduction = "umap", group.by="celltypes", label=TRUE)
DimPlot(immune.combined, reduction = "umap", group.by="leiden", label=TRUE)
DimPlot(immune.combined, reduction = "umap", group.by="celltypeleiden", label=TRUE)
immune.combined@meta.data$`celltypeleiden` <- paste(immune.combined@meta.data$celltypes, immune.combined@meta.data$leiden, sep = ":")

# Trajectory Analysis monocle2
devtools::load_all("D:/BCRdataset/PRJNA890633/monocle_2.30.0/monocle")

DefaultAssay(immune.combined) <- "integrated"
immune.combined 
express_genes <- VariableFeatures(immune.combined)
data <- as(as.matrix(immune.combined@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = immune.combined@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150
HSMM <- setOrderingFilter(HSMM, express_genes)
HSMM <- reduceDimension(HSMM, norm_method="none", max_components = 2, method = "tSNE",auto_param_selection = F) # 可调参数，选择降维方法和主成分
HSMM <- orderCells(HSMM)
p1 <- plot_cell_trajectory(HSMM, color_by = "State")#+ scale_color_manual(values = c("gray","red", "blue"))
p1
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/MBCS_state.pdf", p1, width = 8, height = 8, dpi = 100)
p1 <- plot_cell_trajectory(HSMM, color_by = "celltypes") + 
  scale_color_manual(values = c("#ee4035","#f37736","#0392cf", "#fdf498", "#D3D3D300","#7bc043","purple"))
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/MBCS_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)

HSMM$leiden <- as.character(HSMM$leiden)
p1 <- plot_cell_trajectory(HSMM, color_by = "leiden")+ scale_color_manual(values = c( "#9E0142",  # Dark red
                                                                                "#D53E4F",  # Red
                                                                                "#F46D43",  # Orange-red
                                                                                "#FDAE61",  # Orange
                                                                                "#FEE08B",  # Yellow
                                                                                "#E6F598",  # Light yellow-green
                                                                                "#ABDDA4",  # Light green
                                                                                "#66C2A5",  # Green
                                                                                "#3288BD",  # Blue
                                                                                "#D3D3D300" ))
p1

ggsave("D:/unifiedBCR/figures/CNS_GSE203552/p1_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)

HSMM$highlight <- ifelse(HSMM$leiden %in% c("0", "2", "4","8", "7"), HSMM$leiden, "Others")
plot_cell_trajectory(HSMM, color_by = "highlight") + scale_color_manual(values = c("red","blue","purple", "yellow", "green", "#D3D3D300"))

HSMM$highlight <- ifelse(HSMM$leiden %in% c("0", "1", "5"), HSMM$leiden, "Others")
p1 <- plot_cell_trajectory(HSMM, color_by = "highlight") + scale_color_manual(values = c("#9E0142","#D53E4F","#E6F598","#D3D3D300"))
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/MBCS510_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)
HSMM$highlight <- ifelse(HSMM$leiden %in% c("2", "3"), HSMM$leiden, "Others")
p1 <- plot_cell_trajectory(HSMM, color_by = "highlight") + scale_color_manual(values = c("#F46D43","#FDAE61","#D3D3D300"))
p1
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/MBCS32_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)
HSMM$highlight <- ifelse(HSMM$leiden %in% c("4", "6"), HSMM$leiden, "Others")
p1 <- plot_cell_trajectory(HSMM, color_by = "highlight") + scale_color_manual(values = c("#FEE08B","#ABDDA4","#D3D3D300"))
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/MBCS64_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)
HSMM$highlight <- ifelse(HSMM$leiden %in% c("7", "8"), HSMM$leiden, "Others")
p1 <- plot_cell_trajectory(HSMM, color_by = "highlight") + scale_color_manual(values = c("#66C2A5","#3288BD","#D3D3D300"))
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/MBCS78_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)

HSMM$celltype<-immune.combined@meta.data$celltypes
HSMM$celltype_leiden <- paste(HSMM$celltype, HSMM$leiden, sep = "_")
group1 <- c("mBc1_5", "mBc1_6", "mBc1_3") 
group2 <- c("mBc1_1", "mBc1_4", "mBc1_2") 
HSMM$highlight <- case_when(
  HSMM$celltype_leiden %in% group2 ~ "Group2",  
  HSMM$celltype_leiden %in% group1 ~ "Group1",
  TRUE ~ "Others"
)
HSMM$highlight <- factor(HSMM$highlight, 
                         levels = c("Others", "Group1", "Group2"))
HSMM$point_size <- case_when(
  HSMM$highlight == "Group2" ~ 1.5,  # Group2 最大
  HSMM$highlight == "Group1" ~ 1.0,  # Group1 中等
  TRUE ~ 0.5                         # Others 最小
)
color_values <- c(
  "Others" = "#D3D3D300",    
  "Group1" = "#9E0142",      
  "Group2" = "#3288BD"       
)
p1 <- plot_cell_trajectory(HSMM, color_by = "highlight", cell_size = HSMM$point_size) + 
  scale_color_manual(values = color_values)
p1
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/mBc1_distribution_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)
HSMM$highlight <- ifelse(HSMM$celltype_leiden %in% c("mBc1_2"), HSMM$celltype_leiden, "Others")
p1 <- plot_cell_trajectory(HSMM, color_by = "highlight") + scale_color_manual(values = c("#F46D43","#D3D3D300"))
p1
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/mBc1_2_distribution_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)
group1 <- c("mBc1_5", "mBc1_6", "mBc1_3") 
HSMM$highlight <- case_when(
  HSMM$celltype_leiden %in% group1 ~ "Group1",
  TRUE ~ "Others"
)
p1 <- plot_cell_trajectory(HSMM, color_by = "highlight") + scale_color_manual(values = c("#66C2A5","#D3D3D300"))
p1
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/mBc1_563_distribution_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)


HSMM$mBc3_subset <- case_when(
  HSMM$celltype_leiden == "mBc3_0" ~ "mBc3_0",
  HSMM$celltype_leiden == "mBc3_1" ~ "mBc3_1", 
  HSMM$celltype_leiden == "mBc3_2" ~ "mBc3_2",
  HSMM$celltype_leiden == "mBc3_4" ~ "mBc3_4",
  TRUE ~ NA_character_  # 其他細胞設為 NA
)

p1 <- plot_cell_trajectory(HSMM, color_by = "mBc3_subset") + 
  facet_wrap(~mBc3_subset, ncol = 2) +
  theme(strip.text = element_text(size = 12)) +
  scale_color_discrete(na.value = "transparent")  
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/mBc3_distribution_trajectory.pdf", p1, width = 8, height = 8, dpi = 100)

# DE Analysis
DefaultAssay(immune.combined) <- "RNA"
subobj <- subset(immune.combined, celltypes == "mBc3")
subobj$target_group <- ifelse(
  subobj$celltypes == "mBc3" & (subobj$leiden == "7" | subobj$leiden == "8"),  
  "Target", 
  "Other")
deg_results <- FindMarkers(
  object = subobj, 
  ident.1 = "Target", 
  ident.2 = "Other", 
  group.by = "target_group")

p1 <- EnhancedVolcano(
  deg_results,
  lab = rownames(deg_results),            # Gene names
  x = 'avg_log2FC',                       # Log2 fold-change column
  y = 'p_val_adj',                        # Adjusted p-value column
  pCutoff = 0.05,                         # Adjusted p-value cutoff
  FCcutoff = 1)
p1
ggsave("D:/unifiedBCR/figures/CNS_GSE203552/mBc3_cluster78_DEG.pdf", p1, width = 8, height = 8, dpi = 100)


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(enrichR)

Idents(subobj) <- subobj@meta.data$target_group
degenes<-DEenrichRPlot(subobj, ident.1 = "Target", ident.2 = "Other", num.pathway = 15, max.genes=200,
                       enrich.database =  "GO_Biological_Process_2023", return.gene.list = TRUE)
#write.table(degenes,'cluster78_vs_othermBc3.tsv', sep = '\t',row.names=FALSE, quote = FALSE)


gene_list <- deg_results$avg_log2FC
names(gene_list) <- rownames(deg_results)
gene_list <- sort(gene_list, decreasing = TRUE)


gene_list_entrez <- bitr(names(gene_list), 
                         fromType = "SYMBOL",
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)


gene_list_df <- data.frame(SYMBOL = names(gene_list), 
                           avg_log2FC = gene_list)
gene_list_merged <- merge(gene_list_df, gene_list_entrez, by = "SYMBOL")
gene_list_final <- gene_list_merged$avg_log2FC
names(gene_list_final) <- gene_list_merged$ENTREZID
gene_list_final <- sort(gene_list_final, decreasing = TRUE)


gsea_go <- gseGO(geneList = gene_list_final,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
gsea_results <- gsea_go@result
gsea_results_sig <- gsea_results[gsea_results$p.adjust < 0.05, ]


positive_nes <- gsea_results_sig[gsea_results_sig$NES > 0, ]
if(nrow(positive_nes) > 0) {

  gsea_positive <- gsea_go
  gsea_positive@result <- head(positive_nes[order(positive_nes$NES, decreasing = TRUE), ], 15)
  
  p2 <- dotplot(gsea_positive, showCategory = 15, font.size = 11) + 
    ggtitle("GSEA - Top 15 Positive NES Pathways")
  print(p2)
  ggsave("D:/unifiedBCR/figures/CNS_GSE203552/cluster_78_positivepathways.pdf", p2, width = 12, height = 8)
}
negative_nes <- gsea_results_sig[gsea_results_sig$NES < 0, ]
if(nrow(negative_nes) > 0) {

  gsea_negative <- gsea_go
  gsea_negative@result <- head(negative_nes[order(negative_nes$NES, decreasing = FALSE), ], 15)
  
  p3 <- dotplot(gsea_negative, showCategory = 15, font.size = 11) + 
    ggtitle("GSEA - Top 15 Negative NES Pathways")
  print(p3)
  ggsave("D:/unifiedBCR/figures/CNS_GSE203552/cluster_78_negativepathways.pdf", p3, width = 12, height = 8)
}
