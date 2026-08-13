##############################
## External Validation #######
library(dplyr)
counts <- read.delim(
  "GSE176290_raw_counts_GRCh38.p13_NCBI.tsv",
  header = TRUE,
  check.names = FALSE
)
################################################################################
meta <- read.csv(
  "GSE176290_SraRunTable.csv", # Download metadata from NCBI SRA
  stringsAsFactors = FALSE
)
meta <- meta %>%
  dplyr::select(Run, covid_severity, GEO_Accession)
meta <- meta %>%
  filter(covid_severity != "Non-COVID")
meta$covid_severity[meta$covid_severity == "Mild COVID without hospitalization"] <- "NON-ICU"

meta$covid_severity[meta$covid_severity == "Severe COVID with ICU"] <- "ICU"
################################################################################
# Sample names to keep
samples_keep <- meta$GEO_Accession

# Keep GeneID + selected samples
counts_sub <- counts %>%
  dplyr::select(GeneID, all_of(samples_keep))

dim(counts_sub)
################################################################################
gene_ids <- as.character(counts_sub$GeneID)

count_matrix <- as.matrix(counts_sub[, -1])

rownames(count_matrix) <- gene_ids

storage.mode(count_matrix) <- "integer"

library(org.Hs.eg.db)
library(AnnotationDbi)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(count_matrix),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)
rownames(count_matrix) <- gene_symbols
count_matrix <- count_matrix[
  !is.na(rownames(count_matrix)) &
    rownames(count_matrix) != "",
]
count_matrix <- count_matrix[
  !duplicated(rownames(count_matrix)),
]
count_matrix <- count_matrix[, meta$GEO_Accession]

all(colnames(count_matrix) == meta$GEO_Accession)
################################################################################
library(edgeR)

dge <- DGEList(counts = count_matrix)

dge <- calcNormFactors(dge, method = "TMM")
expr <- edgeR::cpm(
  dge,
  log = TRUE,
  prior.count = 1
)
################################################################################
class_labels <- factor(
  meta$covid_severity,
  levels = c("NON-ICU", "ICU")
)

names(class_labels) <- meta$GEO_Accession

table(class_labels)
identical(colnames(expr), names(class_labels))
################################################################################
library(pROC)

genes <- c(
  "ATAD3A", "FUNDC1", "HLA-DQA1", "KLF6", "LEPROTL1",
  "MAN2A1", "NENF", "RHOB", "UHMK1", "ZMAT3"
)

results <- lapply(genes, function(g) {
  
  roc.obj <- roc(
    response = class_labels,
    predictor = expr[g, ],
    levels = c("NON-ICU", "ICU"),
    direction = "auto"
  )
  
  auc_val <- as.numeric(auc(roc.obj))
  ci_val <- ci.auc(roc.obj)
  
  # Test AUC against 0.5
  test_auc <- roc.test(
    roc.obj,
    auc = 0.5,
    method = "delong"
  )
  
  data.frame(
    Gene = g,
    AUC = auc_val,
    Lower95 = ci_val[1],
    Upper95 = ci_val[3],
    P.Value = test_auc$p.value
  )
})

roc_results <- bind_rows(results)

# BH correction across the 10 genes
roc_results$FDR <- p.adjust(
  roc_results$P.Value,
  method = "BH"
)

roc_results