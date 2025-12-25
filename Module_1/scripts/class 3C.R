BiocManager::install(c("limma", "AnnotationDbi"))
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))

library(AnnotationDbi)   # Handles annotation and probe–gene mapping
library(limma)           # Performs linear modeling and differential expression
library(dplyr)           # Simplifies data manipulation tasks
library(tibble)          



setwd("C:\\Bioinformatics\\AI_Omics_Internship_2025\\Module_1")

# Load preprocessed expression and phenotype data
load("ShahdZiadMohamedRashad_Class_3B_Assignment.RData")


# check annotation slot of your dataset
annotation(raw_data)
# The auhors used a Brainarray custom annotation CDF instead of the default Affymetrix annotation 


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs_needed <- c("AnnotationDbi","org.Hs.eg.db")
for (p in pkgs_needed) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask=FALSE)

#We'll try several BrainArray versions & mirrors.

ba_versions <- c("25.0.0","24.0.0","23.0.0","22.0.0","21.0.0","20.0.0","19.0.0","18.0.0","15.1.0")
#Package base names for HuGene-1_1-st-v1 ENTREZG
bapkgs <- c(
  "hugene11stv1hsentrezg.db",   # annotation DB (PROBEID -> ENTREZID)
  "hugene11stv1hsentrezgcdf",   # CDF (if you process CEL with affy)
  "hugene11stv1hsentrezgprobe"  # probe-level (some pipelines need it)
)

mirrors <- c(
  # official
  "https://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/%s/entrezg.download/%s%s.tar.gz",
  # widely used mirror/alias (many tools use this)
  "http://mbni.org/customcdf/%s/entrezg.download/%s_%s.tar.gz"
)

try_install_ba <- function(pkg, ver) {
  ok <- FALSE
  for (tpl in mirrors) {
    url <- sprintf(tpl, ver, pkg, ver)
    message("Trying: ", url)
    # Use utils::install.packages so no devtools needed
    tf <- tempfile(fileext = ".tar.gz")
    res <- try(utils::download.file(url, tf, mode="wb", quiet=TRUE), silent=TRUE)
    if (!inherits(res, "try-error")) {
      ins <- try(install.packages(tf, repos=NULL, type="source"), silent=TRUE)
      if (!inherits(ins, "try-error")) { ok <- TRUE; break }
    }
  }
  ok
}
installed <- loadedNamespaces()
for (ver in ba_versions) {
  # if DB already installed we can stop early
  if (all(bapkgs %in% rownames(installed.packages()))) break
  for (pkg in bapkgs) {
    if (!(pkg %in% rownames(installed.packages()))) {
      ok <- try_install_ba(pkg, ver)
      if (ok) message("Installed ", pkg, " (", ver, ")")
    }
  }
}
annotate_hugene11 <- function(ids, prefer_geo_fallback = TRUE) {
  stopifnot(is.character(ids))
  out <- data.frame(ID = ids, stringsAsFactors = FALSE)

  is_ba_entrezg <- mean(grepl("^[0-9]+_at$", ids)) > 0.9
  is_transcriptcluster <- mean(grepl("^[0-9]{7}$", ids)) > 0.9

  # Always keep org.Hs.eg.db around for SYMBOL lookups
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
  })

  if (is_ba_entrezg) {
    # BrainArray ENTREZG IDs -> strip "_at" and map via Entrez
    entrez <- sub("_at$", "", ids)
    dfE <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = unique(entrez),
                                 columns = c("SYMBOL","GENENAME"),
                                 keytype = "ENTREZID")
    dfE$ID <- paste0(dfE$ENTREZID, "_at")
    ann <- dfE[match(ids, dfE$ID), c("SYMBOL","GENENAME")]
    out$SYMBOL   <- ann$SYMBOL
    out$GENENAME <- ann$GENENAME

  } else if (is_transcriptcluster) {
    # Standard Affy transcript-cluster IDs -> use Bioconductor's hugene11sttranscriptcluster.db
if (!requireNamespace("hugene11sttranscriptcluster.db", quietly = TRUE)) {
  BiocManager::install("hugene11sttranscriptcluster.db", ask = FALSE)
}
suppressPackageStartupMessages(library(hugene11sttranscriptcluster.db))
df <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
                            keys = unique(ids),
                            columns = c("SYMBOL","GENENAME","ENTREZID"),
                            keytype = "PROBEID")
df <- df[!duplicated(df$PROBEID), ]
ann <- df[match(ids, df$PROBEID), c("SYMBOL","GENENAME","ENTREZID")]
out$SYMBOL   <- ann$SYMBOL
out$GENENAME <- ann$GENENAME
out$ENTREZID <- ann$ENTREZID

} else if (prefer_geo_fallback) {
  # Fallback for BrainArray HuGene11stv1 ENTREZG using GEO platform GPL22654 (v15.1.0)
  if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery", ask=FALSE)
  suppressPackageStartupMessages(library(GEOquery))
  gpl <- GEOquery::getGEO("GPL22654", AnnotGPL = TRUE)
  tab <- Biobase::pData(gpl)
  tab <- tab[, intersect(c("ID","ENTREZ_GENE_ID","Gene_Symbol","Description"), colnames(tab))]
  # Merge by ID directly (these IDs are BrainArray ENTREZG like '1_at')
  ann <- tab[match(ids, tab$ID), ]
  out$SYMBOL   <- ann$Gene_Symbol
  out$GENENAME <- ann$Description
  out$ENTREZID <- ann$ENTREZ_GENE_ID
  
  # If still missing symbols but have ENTrez IDs, upgrade from org.Hs.eg.db
  need <- which(is.na(out$SYMBOL) & !is.na(out$ENTREZID))
  if (length(need)) {
    e2s <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = as.character(out$ENTREZID[need]),
                                 keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
    out$SYMBOL[need] <- unname(e2s)
  }
} else {
  stop("Unrecognized ID format and GEO fallback disabled.")
}

out
}


# Extract probe IDs
ids <- rownames(processed_data)  # or use rownames(processed_data) if exprs() fails

# Run the annotation function
annotation_table <- annotate_hugene11(ids)

# Summarize number of probes per gene symbol
duplicate_summary <- annotation_table %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))


# Identify genes associated with multiple probes
duplicate_genes <- duplicate_summary %>%
  filter(probes_per_gene > 1)

sum(duplicate_genes$probes_per_gene)

## It showed that about 7619 total probes mapped to duplicated genes 
7619 - 5104  # 2515 total probes

# check no. of total probes and % of duplicates
total_probes <- nrow(annotation_table)
multi_probe_symbols <- duplicate_genes$SYMBOL
duplicate_probes <- annotation_table %>%
  filter(SYMBOL %in% multi_probe_symbols)

num_duplicate_probes <- nrow(duplicate_probes)

percent_duplicate <- (num_duplicate_probes / total_probes) * 100

cat("Total probes:", total_probes, "\n")
cat("Probes in multi-probe genes:", num_duplicate_probes, "\n")
cat("Percentage of duplicates:", round(percent_duplicate, 2), "%\n")


# Total probes: 20233 
# Probes in multi-probe genes: 7619 
# Percentage of duplicates: 37.66 %


# Verify if probe IDs in mapping correspond to expression data
all(annotation_table$ID == row.names(processed_data))


# Merge annotation (SYMBOL) with expression matrix
processed_data_df <- processed_data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::mutate(SYMBOL = annotation_table$SYMBOL) %>%
  dplyr::relocate(SYMBOL, .after = PROBEID)


# Check if the annotation is aligned
identical(processed_data_df$PROBEID, annotation_table$ID)

# Remove probes without valid gene symbol annotation
processed_data_df <- processed_data_df %>%
  dplyr::filter(!is.na(SYMBOL))

# Select only numeric expression columns
expr_only <- processed_data_df %>%
  dplyr::select(-PROBEID, -SYMBOL)

# Collapse multiple probes per gene using average expression using the limma package
averaged_data <- limma::avereps(expr_only, ID = processed_data_df$SYMBOL)

dim(averaged_data)


# Convert averaged expression data to matrix format
data <- as.data.frame(averaged_data)
data <- data.matrix(data)
str(data)        # Structure check
is.numeric(data) # Confirm numeric matrix

# Define experimental groups 
groups <- factor(sdrf_file$Factor.Value.compound.,
                 levels = c("atorvastatin", "none"),
                 label = c("treated", "diseased"))

class(groups)
levels(groups)

# Create design matrix for linear modeling
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)


groups <- factor(groups)   # where group = your sample condition vector
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
design

# Fit linear model to expression data
fit_1 <- lmFit(data, design)

contrast_matrix <- makeContrasts(treated_vs_diseased = treated - diseased,
                                 levels = design)

# Apply contrasts and compute moderated statistics
fit_contrast <- contrasts.fit(fit_1, contrast_matrix)

fit_2 <- eBayes(fit_contrast)

# Extract list of differentially expressed genes (DEGs)
deg_results <- topTable(fit_2,
                        coef = "treated_vs_diseased",  # Specify contrast of interest
                        number = Inf,               # Return all genes
                        adjust.method = "BH") # Benjamini-Hochberg correction


# Classify DEGs into Upregulated, Downregulated, or Not Significant
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 0.3, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -0.3, "Downregulated",
         "No")
))

# Subset genes by regulation direction
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")

## This dataset didn't show any up- or downreguated genes for 
# an adj. p-value of 0.05 which means it didn't survive multiple testing
# Atorvastatin "drug" showed no significant change on patients with clinically isolated syndrome



