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

datapath <- "D:/BCRdataset/PRJNA890633/integrated"
embpath <- "D:/unifiedBCR/benchmark/Embeddings/PRJNA890633/mean_inner01.pth"
outdir <- "D:/BCRdataset/PRJNA890633"

bcrfile <- read.table(file=file.path(datapath, "processed_bcr.csv"), header = TRUE, row.names="barcode", sep=",")
annotation <- read.table(file=file.path(outdir, "annotation.csv"), sep=",", header=TRUE, row.names="barcode")
gap <- read.table(file=file.path(outdir, "gap.csv"), header = TRUE, row.names="barcode", sep=",")
bcremb <- read.table(file=file.path(embpath, "embedding.csv"), header = TRUE, sep=",")
gexemb <- read.table(file=file.path(embpath, "gexembedding.csv"), header = TRUE, sep=",")

immune.combined <- readRDS(file=file.path(outdir, "rds_save/immune_combined.rds"))
immune.combined <- AddMetaData(immune.combined, annotation)
immune.combined <- AddMetaData(immune.combined, gap)
immune.combined@meta.data$cdr3 <- bcrfile$cdr3
# Concatenate the embeddings along axis 1
embeddings <- cbind(bcremb, gexemb)
row_norms <- sqrt(rowSums(embeddings^2))
normalized_matrix <- embeddings / row_norms

rownames(normalized_matrix) <- rownames(bcrfile)
inte_assay <- CreateAssayObject(counts = t(normalized_matrix))
immune.combined[["inteemb"]] <- inte_assay
DefaultAssay(immune.combined) <- "inteemb"

DefaultAssay(immune.combined) <- "inteemb"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE, features=rownames(immune.combined))
immune.combined <- FindNeighbors(immune.combined, dims = 1:30, annoy.metric = "cosine")
immune.combined <- FindClusters(immune.combined, resolution = 0.5, annoy.metric = "cosine")
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
DimPlot(immune.combined, reduction = "umap", group.by="batch", label=TRUE)
annotation['seurat_clusters'] <- immune.combined@meta.data['seurat_clusters']
annotation <- data.frame(barcode=rownames(annotation),annotation)
write.table(annotation, file.path(outdir, "annotation.csv"), 
            sep = ',', row.names = F, col.names = T, quote = T) #有commas一定要quote=T
grouped_df <- immune.combined@meta.data %>%
  group_by(celltypes, difference) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = difference, values_from = count, values_fill = 0)
grouped_df['ratio'] <- grouped_df['Large']/(grouped_df['Large'] + grouped_df['Small'])


selectdf <- immune.combined@meta.data[immune.combined@meta.data$seurat_clusters == 15, ]
table(selectdf$celltypes)

markers <-FindAllMarkers(immune.combined, assay="RNA")
spmarkers <- markers[markers$cluster==15,]
immune.combined@meta.data$selected <- 'No'
selected_cells <- rownames(immune.combined@meta.data[immune.combined@meta.data$seurat_clusters %in% c(15,22),])
immune.combined@meta.data[selected_cells,]$selected <- 'Yes'
immune.combined <- SetIdent(immune.combined, value = "selected")
markers <- FindMarkers(immune.combined, ident.1 = "Yes", ident.2 = "No", assay="RNA")

immune.combined <- SetIdent(immune.combined, value = "difference")
markers <- FindMarkers(immune.combined, ident.1 = "Large", assay="RNA")

#Trajectory
subobj <- subset(immune.combined, subset = celltypes == c("Memory"))
subobj <- subset(subobj, subset = batch == c(6))
subobj <- subset(subobj, subset = sample %in% c("CNMC89", "CNMC71"))
DefaultAssay(subobj) <- "integrated"
express_genes <- VariableFeatures(subobj)
data <- as(as.matrix(subobj@assays$integrated@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = subobj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
HSMM <- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = uninormal())# since I have already normalized, thresholded and scalled in Suerat v3.0.0.9150
HSMM <- setOrderingFilter(HSMM, express_genes)
HSMM <- reduceDimension(HSMM, norm_method="none", max_components = 2, method = "DDRTree",
                        ncenter = 300, auto_param_selection = F) 
HSMM <- orderCells(HSMM)
st <- as.character(pData(HSMM)$State)
map <- c("1"="7", "2"="5", "3"="3", "4"="1", "5"="2", "6"="4", "7"="6")
pData(HSMM)$State <- factor(map[st], levels = as.character(1:7))
HSMM <- orderCells(HSMM, root_state = 1)
plot_cell_trajectory(HSMM, color_by = "State")
p1 <- plot_cell_trajectory(HSMM, color_by = "State")
p1
ggsave("D:/unifiedBCR/figures/PRJNA890633/MBCS_state.pdf", p1, width = 8, height = 8, dpi = 100)

plot_cell_trajectory(HSMM, color_by = "difference") + facet_wrap(~State, nrow = 1)
p2 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")
p2
ggsave("D:/unifiedBCR/figures/PRJNA890633/MBCS_pseudotime.pdf", p2, width = 8, height = 8, dpi = 100)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")
p2 <- plot_complex_cell_trajectory(HSMM, x = 1, y = 2,
                                   color_by = "difference")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank()) 
p2
custom_colors <- c("#F39C12", "#9B59B6")
p2 <- plot_cell_trajectory(HSMM, color_by = "difference")+scale_color_manual(values = custom_colors)
p2
ggsave("D:/unifiedBCR/figures/PRJNA890633/MBCS_trajectory.pdf", p2, width = 8, height = 8, dpi = 100)

memoryS_state <- data.frame(barcode=rownames(HSMM@phenoData@data),HSMM@phenoData@data)
stopifnot(identical(colnames(subobj), rownames(memoryS_state)))
write.table(memoryS_state, file.path("D:/BCRdataset/PRJNA890633", "integrated/memoryS_state.csv"), 
            sep = ',', row.names = F, col.names = T, quote = T) 


