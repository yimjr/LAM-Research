#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(affxparser)
  library(limma)
  library(AnnotationDbi)
  library(hta20transcriptcluster.db)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
chp_dir <- file.path(project_root, "data/raw/GSE104335/extracted_gene_expression")
out_dir <- file.path(project_root, "results/mechanisms")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chp_files <- list.files(chp_dir, pattern = "sst-rma-gene-full\\.chp$",
                        full.names = TRUE, ignore.case = TRUE)
if (length(chp_files) != 9L) {
  stop("Expected 9 processed gene-level CHP files, found ", length(chp_files), " in ", chp_dir)
}

gsm <- sub("_.*$", "", basename(chp_files))
gsm_num <- as.integer(sub("GSM", "", gsm))
ord <- order(gsm_num)
chp_files <- chp_files[ord]
gsm <- gsm[ord]
gsm_num <- gsm_num[ord]
group <- ifelse(gsm_num %in% 2795553:2795555, "shGFP_vehicle",
         ifelse(gsm_num %in% 2795556:2795558, "shGFP_rapamycin",
         ifelse(gsm_num %in% 2795559:2795561, "shSRPK2_vehicle", NA)))
if (anyNA(group)) stop("Could not assign a group to every CHP file")

sample_design <- data.frame(sample = gsm, gsm = gsm, group = group,
                            stringsAsFactors = FALSE)
write.csv(sample_design,
          file.path(out_dir, "GSE104335_CHP_sample_design.csv"), row.names = FALSE)

message("Reading Expression Console sst-rma-gene-full CHP files")
chp <- lapply(chp_files, affxparser::readChp)
probe_ids <- chp[[1]]$QuantificationEntries$ProbeSetName
probe_expr <- do.call(cbind, lapply(chp, function(x) x$QuantificationEntries$QuantificationValue))
colnames(probe_expr) <- gsm
rownames(probe_expr) <- probe_ids

algorithm <- unique(vapply(chp, function(x) x$AlgorithmName, character(1)))
array_type <- unique(vapply(chp, function(x) x$ArrayType, character(1)))
if (!all(vapply(chp, function(x) identical(x$QuantificationEntries$ProbeSetName, probe_ids), logical(1)))) {
  stop("CHP files do not have identical probe-set ordering")
}

message("Mapping HTA 2.0 transcript clusters to gene symbols")
probe_map <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = rownames(probe_expr),
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)
probe_map <- probe_map[!is.na(probe_map$SYMBOL) & probe_map$SYMBOL != "", ]
probe_map <- probe_map[!duplicated(probe_map$PROBEID), c("PROBEID", "SYMBOL", "GENENAME")]
map_idx <- match(rownames(probe_expr), probe_map$PROBEID)
keep <- !is.na(map_idx)
probe_expr <- probe_expr[keep, , drop = FALSE]
symbols <- probe_map$SYMBOL[map_idx[keep]]

gene_expr <- rowsum(probe_expr, group = symbols, reorder = FALSE)
gene_counts <- table(symbols)
gene_expr <- gene_expr / as.numeric(gene_counts[rownames(gene_expr)])
gene_expr <- as.matrix(gene_expr)
gene_expr <- gene_expr[order(rownames(gene_expr)), , drop = FALSE]

gene_out <- data.frame(gene_symbol = rownames(gene_expr), gene_expr,
                       check.names = FALSE, row.names = NULL)
write.csv(gene_out,
          file.path(out_dir, "GSE104335_gene_expression_log2_chp_rma.csv"),
          row.names = FALSE)

group_factor <- factor(group,
                       levels = c("shGFP_vehicle", "shGFP_rapamycin", "shSRPK2_vehicle"))
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- levels(group_factor)
contrasts <- makeContrasts(
  shGFP_rapamycin_minus_shGFP_vehicle = shGFP_rapamycin - shGFP_vehicle,
  shSRPK2_vehicle_minus_shGFP_vehicle = shSRPK2_vehicle - shGFP_vehicle,
  shSRPK2_vehicle_minus_shGFP_rapamycin = shSRPK2_vehicle - shGFP_rapamycin,
  levels = design
)
fit <- eBayes(contrasts.fit(lmFit(gene_expr, design), contrasts))

contrast_tables <- lapply(seq_len(ncol(contrasts)), function(i) {
  tab <- topTable(fit, coef = i, number = Inf, sort.by = "none")
  tab$gene_symbol <- rownames(tab)
  tab$contrast <- colnames(contrasts)[i]
  tab[, c("gene_symbol", "contrast", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
})
contrast_out <- do.call(rbind, contrast_tables)
rownames(contrast_out) <- NULL
write.csv(contrast_out,
          file.path(out_dir, "GSE104335_gene_level_contrasts.csv"),
          row.names = FALSE)

write.csv(data.frame(
  array_type = array_type,
  algorithm = algorithm,
  n_transcript_clusters = nrow(probe_expr),
  n_genes_after_mapping = nrow(gene_expr),
  n_samples = ncol(gene_expr),
  stringsAsFactors = FALSE
), file.path(out_dir, "GSE104335_CHP_processing_summary.csv"), row.names = FALSE)

message("Wrote CHP-derived gene-level expression and moderated contrasts to ", out_dir)
