library(Seurat)
library(Matrix)
library(dplyr)
library(edgeR)
library(limma)
################################################################################
cluster_use <- "0" #onlyclassical monocytes

data_sub <- subset(
  data2_mir, # data2_mir can be generated in miRNA_deconv_pipeline.R, it can be exported as a .rds file and imported before proceeding
  subset = seurat_clusters == cluster_use
)
DefaultAssay(data_sub) <- "RNA"
rna <- GetAssayData(
  data_sub,
  assay = "RNA",
  slot = "counts"
)
meta <- data_sub@meta.data

patient_info <- meta %>%
  distinct(ID, condition)
patient_ids <- patient_info$ID

rna_pb <- sapply(patient_ids,function(x){
  
  cells <- rownames(meta)[meta$ID==x]
  
  Matrix::rowSums(
    rna[,cells,drop=FALSE]
  )
  
})

colnames(rna_pb) <- patient_ids

dim(rna_pb)
group <- factor(
  patient_info$condition,
  levels=c("NON-ICU","ICU")
)

design <- model.matrix(~0+group)

colnames(design) <- c("NON_ICU","ICU")
dge <- DGEList(rna_pb)

dge <- calcNormFactors(dge)

v <- voom(dge,design)
fit <- lmFit(v,design)

cont <- makeContrasts(
  ICUvsNON = ICU-NON_ICU,
  levels=design
)

fit2 <- contrasts.fit(
  fit,
  cont
)

fit2 <- eBayes(fit2)
################################################################################


deg_patient <- topTable(
  fit2,
  coef = 1,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg_patient <- deg_patient %>%
  mutate(
    significance = case_when(
      adj.P.Val < 0.05 & logFC > 0.5 ~ "UP",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "DOWN",
      TRUE ~ "Not Significant"
    )
  )
################################################################################
################################################################################
# Create volcano plot
volcano_plot_deg <- ggplot(deg_patient, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = significance), alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "UP" = "red",
    "DOWN" = "blue", 
    "Not Significant" = "grey"
  )) +
  theme_minimal() +
  labs(
    title = "Volcano Plot - DEGs in Classical Monocytes",
    subtitle = "ICU vs NON-ICU",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Regulation"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave("./New_plots/volcano_deg_classical_monocytes.pdf",
       plot = volcano_plot_deg,
       width = 5,
       height = 8,
       dpi = 300)
################################################################################
library(pheatmap)
library(RColorBrewer)

#-----------------------------
# Obtain normalized expression
#-----------------------------
expr <- v$E      # voom normalized log2 CPM values

library(tibble)

deg_patient2 <- deg_patient %>%
  rownames_to_column("Gene")
#-----------------------------
# Select significant DEGs
#-----------------------------
# Significant genes
deg_sig <- deg_patient2 %>%
  filter(significance != "Not Significant")


# Top 25 upregulated
up_genes <- deg_sig %>%
  filter(significance == "UP") %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 25) %>%
  pull(Gene)

# Top 25 downregulated
down_genes <- deg_sig %>%
  filter(significance == "DOWN") %>%
  arrange(logFC) %>%
  slice_head(n = 25) %>%
  pull(Gene)

heatmap_genes <- c(up_genes, down_genes)

length(heatmap_genes)

#-----------------------------
# Expression matrix
#-----------------------------
heatmap_matrix <- expr[heatmap_genes, ]

# Z-score by gene
heatmap_matrix <- t(scale(t(heatmap_matrix)))

#-----------------------------
# Annotation
#-----------------------------
annotation_col <- data.frame(
  Condition = patient_info$condition
)

rownames(annotation_col) <- patient_info$ID

annotation_colors <- list(
  Condition = c(
    "ICU" = "#D73027",
    "NON-ICU" = "#4575B4"
  )
)

#-----------------------------
# Heatmap
#-----------------------------
pdf("./New_plots/heatmap_deg_classical_monocytes.pdf",
    width = 8,
    height = 10)

pheatmap(
  heatmap_matrix,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 7,
  fontsize_col = 9,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  main = "Top Differentially Expressed Genes\nClassical Monocytes"
)

dev.off()
################################################################################
DefaultAssay(data_sub) <- "miRNA"

mir <- GetAssayData(
  data_sub,
  assay = "miRNA",
  slot = "data"
)

meta <- data_sub@meta.data

patient_ids <- unique(meta$ID)

mir_pb <- sapply(patient_ids, function(x){
  
  cells <- rownames(meta)[meta$ID == x]
  
  Matrix::rowMeans(
    mir[, cells, drop = FALSE]
  )
  
})

colnames(mir_pb) <- patient_ids
mir_sd <- apply(mir_pb, 1, sd)

summary(mir_sd)

mir_pb <- mir_pb[
  mir_sd > 0.05,
]
deg_sig <- deg_patient %>%
  filter(
    P.Value < 0.05,
    abs(logFC) > 0.5
  )

genes_use <- rownames(deg_sig)

rna_use <- v$E[genes_use, ]

common_patients <- intersect(
  colnames(rna_use),
  colnames(mir_pb)
)

rna_use <- rna_use[, common_patients]

mir_pb <- mir_pb[, common_patients]
mir_pb <- mir_pb[
  startsWith(rownames(mir_pb), "hsa-"),
]
################################################################################
results <- list()

k <- 1

for(i in 1:nrow(mir_pb)){
  
  mir_name <- rownames(mir_pb)[i]
  mir_exp <- as.numeric(mir_pb[i, ])
  
  for(j in 1:nrow(rna_use)){
    
    gene_name <- rownames(rna_use)[j]
    gene_exp <- as.numeric(rna_use[j, ])
    
    ct <- cor.test(
      mir_exp,
      gene_exp,
      method = "spearman",
      exact = FALSE
    )
    
    rho <- unname(ct$estimate)
    
    # 95% CI using Fisher transformation
    z <- atanh(rho)
    se <- 1 / sqrt(length(mir_exp) - 3)
    
    CI_lower <- tanh(z - 1.96 * se)
    CI_upper <- tanh(z + 1.96 * se)
    
    results[[k]] <- data.frame(
      miRNA = mir_name,
      Gene = gene_name,
      rho = rho,
      CI_lower = CI_lower,
      CI_upper = CI_upper,
      P.Value = ct$p.value
    )
    
    k <- k + 1
  }
}

corr_table <- bind_rows(results)

corr_table$adj.P.Val <- p.adjust(
  corr_table$P.Value,
  method = "BH"
)

corr_all <- corr_table %>%
  filter(
    grepl("^hsa-", miRNA),
    P.Value < 0.05,
    abs(rho) > 0.5
  ) %>%
  arrange(P.Value)

dim(corr_all)
head(corr_all)

corr_negative <- corr_table %>%
  filter(
    grepl("^hsa-", miRNA),
    rho < -0.8,
    P.Value < 0.05
  ) %>%
  arrange(P.Value)

dim(corr_negative)
head(corr_negative)
write.csv(corr_all,
          "All_miRNA_DEG_correlations.csv",
          row.names = FALSE)

write.csv(corr_negative,
          "Negative_miRNA_DEG_correlations.csv",
          row.names = FALSE)
################################################################################
library(multiMiR)
significant_miRNAs <- rownames(mir_pb)
targets <- get_multimir(org = "hsa", 
                        mirna = significant_miRNAs,
                        #table = "all",
                        table = c("validated"),
                        predicted.cutoff = 50)

validated_targets <- targets@data[targets@data$type == "validated", ]
validated_targets <- validated_targets %>%
  distinct(mature_mirna_id, target_symbol, .keep_all = TRUE)
################################################################################
validated_negative_pairs <- corr_negative %>%
  inner_join(validated_targets, 
             by = c("miRNA" = "mature_mirna_id", "Gene" = "target_symbol"))

print(paste("Validated negative correlation pairs:", nrow(validated_negative_pairs)))
write.csv(validated_negative_pairs,
          "Negative_miRNA_DEG_correlations_validated.csv",
          row.names = FALSE)
################################################################################