# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
# Function to create a Seurat object and update cell barcodes
create_seurat_object <- function(gsm_id, base_dir, metadata) {
  # Path to the data directory
  data_dir <- file.path(base_dir, gsm_id)
  
  # Read 10X data
  tenx <- Read10X(data.dir = data_dir, gene.column = 1, cell.column = 1, unique.features = TRUE)
  counts <- tenx$`Gene Expression`
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)
  seurat_obj$nCount_RNA <- Matrix::colSums(counts)
  seurat_obj$nFeature_RNA <- Matrix::colSums(counts > 0)
  # Ensure metadata has rownames matching updated cell barcodes
  metadata$Cell_barcode <- make.unique(as.character(metadata$Cell_barcode))
  rownames(metadata) <- metadata$Cell_barcode
  # Add 'ID' metadata to the Seurat object
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata["ID"])
  
  return(seurat_obj)
}

# Base directory containing all GSM folders
base_dir <- "/path/to/your/directory"

# Load metadata (replace 'metadata.csv' with your actual file path)
metadata <- read.csv("metadata.csv")

# List of GSM IDs ||| make sure each folder has three files for importing (barcode, matrix, features)
gsm_ids <- c("GSM7349769", "GSM7349770", "GSM7349771", "GSM7349773", "GSM7349775", 
             "GSM7349780", "GSM7349786", "GSM7349787")

# Create a list to hold all Seurat objects
seurat_list <- list()

# Loop through GSM IDs and create Seurat objects
for (gsm_id in gsm_ids) {
  seurat_list[[gsm_id]] <- create_seurat_object(gsm_id, base_dir, metadata)
}

# Merge all Seurat objects into a single object
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)

# View metadata in the merged Seurat object
View(merged_seurat@meta.data)
#head(seurat_obj@meta.data)
# Identify cells with non-NA 'ID' values
valid_cells <- which(!is.na(merged_seurat@meta.data$ID))
valid_cell_names <- rownames(merged_seurat@meta.data)[valid_cells]

# Filter the Seurat object to exclude cells with NA in the 'ID' column
# Subset the Seurat object using valid cell names
merged_seurat2 <- subset(merged_seurat, cells = valid_cell_names)


# Define the IDs to exclude
exclude_ids <- c("029-Monocytes", "033-Monocytes", "037-Monocytes", "048-Monocytes", "061-Monos-Tregs", "100-Monos-Tregs")

# Identify cells that do not have the excluded IDs
valid_cells <- which(!is.na(merged_seurat@meta.data$ID) & !merged_seurat@meta.data$ID %in% exclude_ids)

# Get cell names corresponding to the valid cells
valid_cell_names <- rownames(merged_seurat@meta.data)[valid_cells]

# Subset the Seurat object using valid cell names
merged_seurat3<- subset(merged_seurat, cells = valid_cell_names)

# Verify the unique IDs in the filtered Seurat object
unique(merged_seurat3@meta.data$ID)



# Define the ICU IDs
icu_ids <- c("034-Monocytes", "036-Monocytes", "050-Monocytes", "101-Monos-Tregs", "103-Monos-Tregs")

# Add a new 'condition' column to the metadata of merged_seurat3
merged_seurat3@meta.data$condition <- ifelse(merged_seurat3@meta.data$ID %in% icu_ids, "ICU", "NON-ICU")

# Verify the new column
merged_seurat3@meta.data$condition
merged_seurat3[["percent.mt"]] <- PercentageFeatureSet(merged_seurat3, pattern = "mt-")
lb <- quantile(merged_seurat[["nFeature_RNA"]]$nFeature_RNA, probs = 0.01)
ub <- quantile(merged_seurat[["nFeature_RNA"]]$nFeature_RNA, probs = 0.99)
merged_seurat <- merged_seurat3[, merged_seurat3[["nFeature_RNA"]] > lb & 
                                  merged_seurat3[["nFeature_RNA"]] < ub & merged_seurat3[["percent.mt"]] < 15]

merged_seurat <- NormalizeData(object = merged_seurat, verbose = FALSE)
merged_seurat <- FindVariableFeatures(object = merged_seurat, nfeatures = 3000, 
                                      verbose = FALSE, selection.method = 'vst')
merged_seurat <- ScaleData(merged_seurat, verbose = FALSE)
merged_seurat <- RunPCA(merged_seurat, npcs = 20, verbose = FALSE)
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:20)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.08)
merged_seurat <- RunUMAP(merged_seurat, reduction = "pca", dims = 1:20)
merged_seurat <- subset(merged_seurat, idents = "6", invert = TRUE)
clusters_woa <- DimPlot(
  merged_seurat,
  reduction = "umap",
  label = TRUE
) + theme(
  text = element_text(size = 18, face = "bold"),
  plot.title = element_text(face = "bold")
)

clusters_condition <- DimPlot(merged_seurat, reduction = 'umap', group.by = "condition") + theme(
  text = element_text(size = 18, face = "bold"),
  plot.title = element_text(face = "bold")
)

ggsave("./plots/UMAP_no_anno.pdf",
       plot = clusters_woa,
       width = 8,
       height = 6,
       dpi = 300)
ggsave("./plots/UMAP_CON.pdf",
       plot = clusters_condition,
       width = 8,
       height = 6,
       dpi = 300)

###################################################################
########## Gene ID conversion #########################

data <- merged_seurat

library(biomaRt)
library(httr)
set_config(config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE))

Sys.setenv(BIOMART_CACHE = tempdir())
# Use the Ensembl database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert ENSG IDs to Gene Symbols
#genes <- rownames(data)
genes <- rownames(data[["RNA"]]@counts)
conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)

gene_map <- data.frame(ENSG = genes) %>%
  left_join(conversion, by = c("ENSG" = "ensembl_gene_id"))

gene_map$hgnc_symbol[
  is.na(gene_map$hgnc_symbol) | gene_map$hgnc_symbol == ""
] <- gene_map$ENSG[
  is.na(gene_map$hgnc_symbol) | gene_map$hgnc_symbol == ""
]

gene_map <- gene_map[!duplicated(gene_map$hgnc_symbol), ]

gene_map_vector <- setNames(
  gene_map$hgnc_symbol,
  gene_map$ENSG
)
#############################################################
expr <- data[["RNA"]]@data

new_names <- gene_map_vector[rownames(expr)]
rownames(expr) <- new_names
# Replace rownames in Seurat object
rownames(expr) <- gene_map_vector[rownames(expr)]
######################################################################

######################################################################

#data1 <- JoinLayers(expr)
data1 <- NormalizeData(data)
data1 <- FindVariableFeatures(data1)
data1 <- ScaleData(data1)

# Find markers for each cluster
markers <- FindAllMarkers(data1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# View the top markers
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
View(top_markers)

# Heatmap of top marker genes
library(ggplot2)
htmap <- DoHeatmap(data1, features = top_markers$gene) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(
    text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )
ggsave("./plots/HTMAP_GENE.pdf",
       plot = htmap,
       width = 12,
       height = 10,
       dpi = 300)

ftplot_gene <- FeaturePlot(data1, features = c("S100A8", "S100A9",  "IL32", "CD3E",
                                "FCGR3A", "VMO1", "GNLY", "NKG7", 
                                "CD79A", "MS4A1", "FCER1A", "CD1C")) + theme(
                                  text = element_text(face = "bold"),
                                  plot.title = element_text(face = "bold")
                                )

dtplt<- DotPlot(data1, features = c("S100A8", "S100A9",  "IL32", "CD3E",
                                    "FCGR3A", "VMO1", "GNLY", "NKG7", 
                                    "CD79A", "MS4A1", "FCER1A", "CD1C"), cols = c("yellow", "blue")) +
  RotatedAxis() + theme(
    text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave("./plots/FTPLOT_GENE.pdf",
       plot = ftplot_gene,
       width = 12,
       height = 8,
       dpi = 300)
ggsave("./plots/DOTPLOT_GENE.pdf",
       plot = dtplt,
       width = 12,
       height = 6,
       dpi = 300)
##################################################
##### Doublet Detection ##########################
library(scDblFinder)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(data1)

sce <- scDblFinder(sce)

table(sce$scDblFinder.class)
data1$DoubletStatus <- sce$scDblFinder.class
umap_doublet <- DimPlot(
  data1,
  group.by="DoubletStatus"
) + theme(
  text = element_text(size = 18, face = "bold"),
  plot.title = element_text(face = "bold")
)

ggsave("./plots/UMAP_doublet.pdf",
       plot = umap_doublet,
       width = 8,
       height = 6,
       dpi = 300)

table(Idents(data1), data1$DoubletStatus)
prop.table(table(Idents(data1), data1$DoubletStatus), margin = 1)

data1 <- subset(data1, subset = DoubletStatus == "singlet")
#################################################
###### Azimuth annotation #######################
library(Seurat)
library(Azimuth)
data1 <- RunAzimuth(
  query = data1,
  reference = "pbmcref"
)
prop.table(
  table(data1$seurat_clusters, data1$predicted.celltype.l2),
  margin = 1
)
table(data1$predicted.celltype.l2)
table(data1$seurat_clusters, data1$predicted.celltype.l2)
pdf("./plots/azimuth_anno_l2.pdf", width= 10, height=6)
DimPlot(
  data1,
  reduction = "umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE
)
dev.off()

##################################################
# Assign cluster labels based on marker genes
new_cluster_ids <- c(
  "Classical monocytes (CD14+)",     # Cluster 0
  "T Cells",                  # Cluster 1
  "Non-classical monocytes (CD16+)",# Cluster 2
  "Cytotoxic lymphocytes",     # Cluster 3
  "B cells",       # Cluster 4
  "Dendritic cells"        # Cluster 5
)

data2 <- data1
# Update cluster identities
names(new_cluster_ids) <- levels(data2)
data2 <- RenameIdents(data2, new_cluster_ids)

# Visualize with new labels
annotated <- DimPlot(data2, reduction = "umap", pt.size = 1) + theme(
  text = element_text(size = 18, face = "bold"),
  plot.title = element_text(face = "bold")
)
ggsave("./plots/UMAP_annotated.pdf",
       plot = annotated,
       width = 10,
       height = 6,
       dpi = 300)
# Extract raw matrices
counts_mat <- GetAssayData(data2, assay = "RNA", slot = "counts")
data_mat <- GetAssayData(data2, assay = "RNA", slot = "data")
scale_data_mat <- GetAssayData(data2, assay = "RNA", slot = "scale.data")

# Create a Seurat v4-style Assay object manually
assay_v4 <- CreateAssayObject(counts = counts_mat)
assay_v4@data       <- data_mat
assay_v4@scale.data <- scale_data_mat

data2_clean <- CreateSeuratObject(counts = counts_mat)
data2_clean@assays$RNA <- assay_v4
# Transfer metadata, reductions, etc., if needed
data2_clean@meta.data  <- data2@meta.data
data2_clean@reductions <- data2@reductions
data2_clean@graphs     <- data2@graphs

########################################################################
####### Calculate cell type counts per condition #######################
data2$current_cluster_identity <- Idents(data2)

# Check the result
print("Current cluster identities:")
print(unique(data2$current_cluster_identity))

condition_celltype_counts <- table(data2@meta.data$condition, 
                                   data2@meta.data$current_cluster_identity)
unique_conditions <- unique(data2@meta.data$condition)

###########################################################################
#######  Plot for covid with cell counts     ###################

# Subset data for current condition
icu_data <- subset(data2, subset = condition == "ICU")
# Filter cell type counts for current condition 
icu_celltype_counts <- condition_celltype_counts["ICU", ]

# Convert table to data frame 
icu_celltype_counts_df <- as.data.frame(icu_celltype_counts)

# Modify label text with counts for this condition
icu_data$label_text <- paste0(icu_data$current_cluster_identity, 
                              " (", icu_celltype_counts_df[icu_data$current_cluster_identity, 1], ")")

# Generate UMAP plot with labeled cell counts
icu_plot <- DimPlot(icu_data, reduction = "umap", split.by = "condition", group.by = "label_text") + theme(
  text = element_text(size = 18, face = "bold"),
  plot.title = element_text(face = "bold")
)
ggsave("./plots/UMAP_ICU.pdf",
       plot = icu_plot,
       width = 10,
       height = 6,
       dpi = 300)
#################################################################################################
############  Plot for Normal with cellcounts ############################

# Subset data for current condition
nor_data <- subset(data2, subset = condition == "NON-ICU")
# Filter cell type counts for current condition 
nor_celltype_counts <- condition_celltype_counts["NON-ICU", ]

# Convert table to data frame 
nor_celltype_counts_df <- as.data.frame(nor_celltype_counts)

# Modify label text with counts for this condition
nor_data$label_text <- paste0(nor_data$current_cluster_identity, " (", nor_celltype_counts_df[nor_data$current_cluster_identity, 1], ")")

# Generate UMAP plot with labeled cell counts
nor_plot <- DimPlot(nor_data, reduction = "umap", split.by = "condition", group.by = "label_text") + theme(
  text = element_text(size = 18, face = "bold"),
  plot.title = element_text(face = "bold")
)
ggsave("./plots/UMAP_NON_ICU.pdf",
       plot = nor_plot,
       width = 10.5,
       height = 6,
       dpi = 300)
###################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)

patient_prop <- data2@meta.data %>%
  count(ID, condition, current_cluster_identity, name = "CellCount") %>%
  group_by(ID) %>%
  mutate(
    TotalCells = sum(CellCount),
    Proportion = CellCount / TotalCells
  ) %>%
  ungroup()

write.csv(
  patient_prop,
  "./plots/Patient_celltype_proportions.csv",
  row.names = FALSE
)

p_patient_prop <- ggplot(
  patient_prop,
  aes(
    x = condition,
    y = Proportion,
    group = ID
  )
) +
  geom_line(alpha = 0.4) +
  geom_point(size = 3) +
  facet_wrap(
    ~ current_cluster_identity,
    scales = "free_y"
  ) +
  scale_y_continuous(
    labels = scales::percent_format()
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of cells per patient"
  )

ggsave(
  "./plots/Patient_celltype_proportions.pdf",
  p_patient_prop,
  width = 12,
  height = 8,
  dpi = 300
)

library(glmmTMB)
celltypes <- unique(data2$current_cluster_identity)

patient_counts <- data2@meta.data %>%
  count(ID, condition, current_cluster_identity, name = "CellCount") %>%
  complete(
    ID,
    current_cluster_identity = celltypes,
    fill = list(CellCount = 0)
  ) %>%
  left_join(
    data2@meta.data %>%
      distinct(ID, condition),
    by = "ID",
    suffix = c("", ".new")
  ) %>%
  mutate(
    condition = coalesce(condition, condition.new)
  ) %>%
  select(-condition.new) %>%
  group_by(ID) %>%
  mutate(
    TotalCells = sum(CellCount)
  ) %>%
  ungroup()

patient_counts$condition <- factor(
  patient_counts$condition,
  levels = c("NON-ICU", "ICU")
)

beta_results <- lapply(celltypes, function(ct) {
  
  dat <- patient_counts %>%
    filter(current_cluster_identity == ct) %>%
    mutate(
      Success = CellCount,
      Failure = TotalCells - CellCount
    )
  
  model <- glmmTMB(
    cbind(Success, Failure) ~ condition,
    family = betabinomial(link = "logit"),
    data = dat
  )
  
  cf <- summary(model)$coefficients$cond
  
  data.frame(
    CellType = ct,
    Estimate = cf["conditionICU", "Estimate"],
    SE = cf["conditionICU", "Std. Error"],
    P.Value = cf["conditionICU", "Pr(>|z|)"]
  )
}) %>%
  bind_rows()

beta_results <- beta_results %>%
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    OR = exp(Estimate),
    CI_lower = exp(Estimate - 1.96 * SE),
    CI_upper = exp(Estimate + 1.96 * SE)
  )
beta_results

beta_results_save <- beta_results %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2)),
    P.Value = round(P.Value, 4),
    FDR = round(FDR, 4)
  )

write.csv(
  beta_results_save,
  "./plots/Patient_level_celltype_beta_binomial_results.csv",
  row.names = FALSE
)


###################################################################################
marker_list <- list(
  "Classical monocytes (CD14+)" = c(
    "S100A8", "S100A9", "LYZ", "CTSD"
  ),
  
  "Non-classical monocytes (CD16+)" = c(
    "FCGR3A", "LST1", "IFITM3"
  ),
  
  "T Cells" = c(
    "CD3D", "CD3E", "TRBC1"
  ),
  
  "Cytotoxic lymphocytes" = c(
    "NKG7", "GNLY", "GZMB"
  ),
  
  "B cells" = c(
    "CD79A", "MS4A1", "CD74"
  ),
  
  "Dendritic cells" = c(
    "FCER1A", "CD1C", "CST3"
  )
)

marker_genes <- unique(unlist(marker_list))

marker_genes <- marker_genes[
  marker_genes %in% rownames(data2)
]

marker_genes

library(ggplot2)
library(Seurat)

# Create patient-cell-type grouping
data2$Patient_CellType <- paste(
  data2$current_cluster_identity,
  data2$ID,
  sep = " | "
)

# Make sure marker genes are present
marker_genes <- marker_genes[
  marker_genes %in% rownames(data2)
]

p_patient_markers <- DotPlot(
  data2,
  features = marker_genes,
  group.by = "Patient_CellType",
  cols = c("grey90", "blue")
) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      size = 12
    )
  )

ggsave(
  "./plots/Patient_Markers/Patient_wise_celltype_markers.pdf",
  plot = p_patient_markers,
  width = 18,
  height = 14,
  dpi = 300
)

ggsave(
  "./plots/Patient_Markers/Patient_wise_celltype_markers.jpg",
  plot = p_patient_markers,
  width = 18,
  height = 14,
  dpi = 300
)

####### Extrapolation starts here ################
source("/path/to/miRSCAPE.R")
example <- data2
#example[["RNA"]] <- as(example[["RNA"]], "Assay")
example <- ScaleData(example)
colnames(example@meta.data)[colnames(example@meta.data) == "seurat_clusters"] <- "cluster"
denem = modifySeuratObject(pbmc = example, scaled = TRUE)

pdac_mirna <- read.delim("/path/to/mirna_mirscape.txt", row.names=1)
pdac_mrna <- read.delim("/path/to/mrna_mirscape.txt", row.names=1)

bulkk_mirna = bulkTransform(pdac_mirna, justNorm = TRUE)
bulkk_mrna = bulkTransform(pdac_mrna)

pred = miRSCAPE(bulkmRNA = bulkk_mrna, bulkmiRNA = bulkk_mirna, scmRNA = denem)

df_mir <- pred
# Set the first row as column names
#colnames(df_mir) <- df_mir[1, ]

# Remove the first row (since it's now column names)
#df_mir <- df_mir[-1, ]

# Reorder columns
desired_order <- c("0", "1", "2", "3", "4", "5")
pred_reordered <- df_mir[, desired_order]

# View the reordered data
head(pred_reordered)

##################################################

mir_expr <- pred_reordered
# Transpose mir_matrix so: columns = miRNAs, rows = clusters
mir_matrix_t <- as.data.frame(t(mir_expr))  # rownames = clusters
# Add the cluster ID column
mir_matrix_t$seurat_clusters <- rownames(mir_matrix_t)
#tail(mir_matrix_t)

# Fetch metadata from Seurat
meta <- data2_clean@meta.data
meta <- meta %>% 
  arrange(as.numeric(as.character(seurat_clusters)))

meta$cell <- rownames(meta)
meta <- meta %>%
  arrange(seurat_clusters)



# Ensure cluster ID matches
meta$seurat_clusters <- as.character(meta$seurat_clusters)
mir_matrix_t$seurat_clusters <- as.character(mir_matrix_t$seurat_clusters)

# Merge miRNA predictions into the metadata
meta_merged <- merge(meta, mir_matrix_t, by = "seurat_clusters")
meta_merged <- meta_merged[colnames(data2_clean), , drop = FALSE]
identical(rownames(meta_merged), colnames(data2_clean))
rownames(meta_merged) <- meta_merged$cell
# Set a valid key (prefix) for the RNA assay
data2_clean@assays$RNA@key <- "RNA_"
new.cols <- setdiff(colnames(meta_merged), colnames(data2_clean@meta.data))

length(new.cols)
meta_new <- meta_merged[, new.cols, drop = FALSE]
data2_clean@meta.data <- cbind(data2_clean@meta.data, meta_new)

######### Analysis on sc-miRNA data ########################################
mir_sc <- data2_clean


# Automatically detect miRNA columns from meta.data
miRNA_cols <- grep("^hsa-", colnames(mir_sc@meta.data), value = TRUE)
# Extract and transpose: features (miRNAs) as rows, cells as columns
mir_matrix <- t(mir_sc@meta.data[, miRNA_cols])
# Create an Assay object from the miRNA matrix
mir_assay <- CreateAssayObject(counts = mir_matrix)

# Add it to the Seurat object
mir_sc[["miRNA"]] <- mir_assay

# Set the default assay to miRNA
DefaultAssay(mir_sc) <- "miRNA"
# Set identities (grouping) to mRNA-based clusters
Idents(mir_sc) <- mir_sc$seurat_clusters
################################################################

miRNA_markers <- FindAllMarkers(
  mir_sc,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.0
)

miRNA_markers <- miRNA_markers %>%
  filter(!grepl("\\.x$|\\.y$", gene))

# Top 5 per cluster
library(dplyr)
top_miRNAs <- miRNA_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

top_miRNAs <- top_miRNAs %>% 
  arrange(as.numeric(as.character(cluster)))

top_miRNAs

mir_sc <- ScaleData(mir_sc)

htmap <- DoHeatmap(
  object = mir_sc,
  features = unique(top_miRNAs$gene),  # or top_miRNAs_subset$gene
  group.by = "seurat_clusters"         # Adds clusters
) + scale_fill_gradientn(colors = c("blue", "white", "red"))

ggsave("./plots/HTMAP_miRNA.pdf",
       plot = htmap,
       width = 12,
       height = 10,
       dpi = 300)

ftplot_miRNA <- FeaturePlot(mir_sc, features = c("hsa-miR-125a-5p", "hsa-miR-1226-3p", "hsa-miR-592", "hsa-miR-618",
                                 "hsa-miR-425-3p", "hsa-miR-3173-5p", "hsa-miR-2110", "hsa-miR-20b-5p", 
                                 "hsa-miR-4273", "hsa-miR-517-5p", "hsa-let-7b-5p", "hsa-miR-378a-3p"))

dtplot_miRNA <- DotPlot(mir_sc, features = c("hsa-miR-125a-5p", "hsa-miR-1226-3p", "hsa-miR-592", "hsa-miR-618",
                                             "hsa-miR-425-3p", "hsa-miR-3173-5p", "hsa-miR-2110", "hsa-miR-20b-5p", 
                                             "hsa-miR-4273", "hsa-miR-517-5p", "hsa-let-7b-5p", "hsa-miR-378a-3p"), 
        cols = c("yellow", "blue")) +
  RotatedAxis() +
  scale_y_discrete(limits = as.character(0:5)) + theme(
    text = element_text(size = 18, face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave("./plots/FTPLOT_miRNA.pdf",
       plot = ftplot_miRNA,
       width = 12,
       height = 8,
       dpi = 300)
ggsave("./plots/DTPLOT_miRNA.pdf",
       plot = dtplot_miRNA,
       width = 13,
       height = 5,
       dpi = 300)
##################################################################