library(Seurat)
library(celda)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)

dir.create("./New_plots/Ambient_RNA", showWarnings = FALSE, recursive = TRUE)
dir.create("./New_plots/Patient_Markers", showWarnings = FALSE, recursive = TRUE)


ambient_results <- list()

for (gsm_id in names(seurat_list)) {
  
  cat("Running DecontX for:", gsm_id, "\n")
  
  seu <- seurat_list[[gsm_id]]
  
  # Convert Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(seu)
  
  # Run DecontX
  sce <- decontX(sce)
  
  # Store DecontX results
  ambient_results[[gsm_id]] <- data.frame(
    GSM = gsm_id,
    Cell = colnames(sce),
    DecontX_contamination = sce$decontX_contamination,
    DecontX_cluster = sce$decontX_clusters
  )
}

# Combine all samples
ambient_df <- dplyr::bind_rows(ambient_results)

# Save results
write.csv(
  ambient_df,
  "./New_plots/Ambient_RNA/DecontX_cell_contamination.csv",
  row.names = FALSE
)

ambient_df %>%
  group_by(GSM) %>%
  summarise(
    n_cells = n(),
    mean_contamination = mean(DecontX_contamination, na.rm = TRUE),
    median_contamination = median(DecontX_contamination, na.rm = TRUE),
    max_contamination = max(DecontX_contamination, na.rm = TRUE)
  )


library(dplyr)
metadata <- read.csv("metadata.csv")
ambient_df <- ambient_df %>%
  left_join(
    metadata %>%
      select(Cell_barcode, ID),
    by = c("Cell" = "Cell_barcode")
  )

ambient_df <- na.omit(ambient_df)

patient_contamination <- ambient_df %>%
  group_by(ID) %>%
  summarise(
    n_cells = n(),
    mean_contamination =
      mean(DecontX_contamination, na.rm = TRUE),
    median_contamination =
      median(DecontX_contamination, na.rm = TRUE),
    pct_gt_10 =
      mean(DecontX_contamination > 0.10, na.rm = TRUE) * 100,
    pct_gt_20 =
      mean(DecontX_contamination > 0.20, na.rm = TRUE) * 100,
    pct_gt_30 =
      mean(DecontX_contamination > 0.30, na.rm = TRUE) * 100,
    pct_gt_50 =
      mean(DecontX_contamination > 0.50, na.rm = TRUE) * 100,
    .groups = "drop"
  )

patient_contamination
selected_IDs <- unique(data2$ID)
ambient_df_selected <- ambient_df %>%
  filter(ID %in% selected_IDs)

patient_contamination_selected <- ambient_df_selected %>%
  group_by(ID) %>%
  summarise(
    n_cells = n(),
    mean_contamination =
      mean(DecontX_contamination, na.rm = TRUE),
    median_contamination =
      median(DecontX_contamination, na.rm = TRUE),
    pct_gt_10 =
      mean(DecontX_contamination > 0.10, na.rm = TRUE) * 100,
    pct_gt_20 =
      mean(DecontX_contamination > 0.20, na.rm = TRUE) * 100,
    pct_gt_30 =
      mean(DecontX_contamination > 0.30, na.rm = TRUE) * 100,
    pct_gt_50 =
      mean(DecontX_contamination > 0.50, na.rm = TRUE) * 100,
    .groups = "drop"
  )
