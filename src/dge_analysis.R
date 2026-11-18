#!/usr/bin/env -S Rscript --no-init-file
# Run Differential Gene Expression (DESeq2/edgeR), GEO ingestion, GO enrichment, and plots for Lectures 2 & 5.

parse_cli_args <- function(args) {
    prog_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    prog_name <- if (length(prog_arg) > 0) basename(sub("^--file=", "", prog_arg[1])) else "dge_analysis.R"

    print_help <- function() {
        cat(sprintf(paste0(
            "usage: %s [OPTIONS]\n\n",
            "Run DESeq2 or edgeR differential expression analysis, Top-10 DEG Markdown reporting, and enrichment plots.\n\n",
            "options:\n",
            "  -c, --counts PATH                  Gene count matrix (TSV/CSV with gene IDs in column 1)\n",
            "  -m, --metadata PATH                Sample metadata TSV/CSV (columns: sample, condition)\n",
            "  -g, --geo ACCESSION                NCBI GEO series accession (e.g. GSE154132)\n",
            "  --case {cs1,cs2,cs3,cs4a,cs4b,cs5,cs6}\n",
            "                                     Predefined Lecture 5 case study preset\n",
            "  --method {deseq2,edger}            DGE statistical method (default: deseq2)\n",
            "  --lfc FLOAT                        Absolute log2 fold-change threshold (default: 1.0)\n",
            "  --padj FLOAT                       Adjusted p-value significance cutoff (default: 0.05)\n",
            "  --species {hs,mm,ec,at,none}       Organism annotation DB for clusterProfiler GO (default: none)\n",
            "  --synthetic-demo                   Generate and analyze an offline synthetic RNA-seq dataset\n",
            "  -o, --output-dir DIR               Directory for tables and plots (default: lecture-05)\n",
            "  -n, --dry-run                      Preview planned analysis without executing\n",
            "  -h, --help                         Show this help message and exit\n\n",
            "examples:\n",
            "  # Run DESeq2 on Lecture 5 Case Study 2 (E. coli ciprofloxacin response, GSE220559)\n",
            "  %s --case cs2 -o lecture-05/cs2_ecoli\n\n",
            "  # Run DESeq2 on a local featureCounts matrix from rnaseq_pipeline.py\n",
            "  %s -c lecture-02/g1/counts/counts_matrix.tsv -m lecture-02/g1/metadata.tsv -o lecture-02/g1/dge\n\n",
            "  # Verify environment end-to-end with an offline synthetic dataset\n",
            "  %s --synthetic-demo -o .scratch/dge-demo\n"
        ), prog_name, prog_name, prog_name, prog_name))
    }

    opts <- list(
        counts = NULL,
        metadata = NULL,
        geo = NULL,
        case_id = NULL,
        method = "deseq2",
        lfc = 1.0,
        padj = 0.05,
        species = "none",
        synthetic_demo = FALSE,
        output_dir = "lecture-05",
        dry_run = FALSE
    )

    i <- 1
    while (i <= length(args)) {
        flag <- args[i]
        need_val <- function() {
            if (i >= length(args)) {
                cat(sprintf("Error: Option '%s' requires an argument.\nRun '%s --help' for usage.\n", flag, prog_name), file = stderr())
                quit(status = 1)
            }
            args[i + 1]
        }

        if (flag %in% c("-h", "--help")) {
            print_help()
            quit(status = 0)
        } else if (flag %in% c("-c", "--counts")) {
            opts$counts <- need_val()
            i <- i + 2
        } else if (flag %in% c("-m", "--metadata")) {
            opts$metadata <- need_val()
            i <- i + 2
        } else if (flag %in% c("-g", "--geo")) {
            opts$geo <- need_val()
            i <- i + 2
        } else if (flag == "--case") {
            val <- need_val()
            if (!val %in% c("cs1", "cs2", "cs3", "cs4a", "cs4b", "cs5", "cs6")) {
                cat(sprintf("Error: Invalid --case '%s'.\nRun '%s --help' for usage.\n", val, prog_name), file = stderr())
                quit(status = 1)
            }
            opts$case_id <- val
            i <- i + 2
        } else if (flag == "--method") {
            val <- need_val()
            if (!val %in% c("deseq2", "edger")) {
                cat(sprintf("Error: Invalid --method '%s'.\nRun '%s --help' for usage.\n", val, prog_name), file = stderr())
                quit(status = 1)
            }
            opts$method <- val
            i <- i + 2
        } else if (flag == "--lfc") {
            opts$lfc <- as.numeric(need_val())
            i <- i + 2
        } else if (flag == "--padj") {
            opts$padj <- as.numeric(need_val())
            i <- i + 2
        } else if (flag == "--species") {
            val <- need_val()
            if (!val %in% c("hs", "mm", "ec", "at", "none")) {
                cat(sprintf("Error: Invalid --species '%s'.\nRun '%s --help' for usage.\n", val, prog_name), file = stderr())
                quit(status = 1)
            }
            opts$species <- val
            i <- i + 2
        } else if (flag == "--synthetic-demo") {
            opts$synthetic_demo <- TRUE
            i <- i + 1
        } else if (flag %in% c("-o", "--output-dir")) {
            opts$output_dir <- need_val()
            i <- i + 2
        } else if (flag %in% c("-n", "--dry-run")) {
            opts$dry_run <- TRUE
            i <- i + 1
        } else {
            cat(sprintf("Error: Unknown option '%s'.\nRun '%s --help' for usage.\n", flag, prog_name), file = stderr())
            quit(status = 1)
        }
    }

    if (!opts$synthetic_demo && is.null(opts$counts) && is.null(opts$geo) && is.null(opts$case_id)) {
        cat(sprintf("Error: Provide --counts + --metadata, --geo, --case, or --synthetic-demo.\nRun '%s --help' for usage.\n", prog_name), file = stderr())
        quit(status = 1)
    }
    opts
}

CASE_REGISTRY <- list(
    cs1 = list(geo = "GSE33110", title = "Liver cell-to-iPSC reprogramming", species = "mm"),
    cs2 = list(geo = "GSE220559", title = "E. coli purA ciprofloxacin response", species = "ec"),
    cs3 = list(geo = "GSE80448", title = "Arabidopsis thaliana RPS4 immune response", species = "at"),
    cs4a = list(geo = "GSE120602", title = "Pseudomonas aeruginosa clinical isolates", species = "none"),
    cs4b = list(geo = "GSE154132", title = "Mus musculus diabetic wound healing", species = "mm"),
    cs5 = list(geo = "GSE193336", title = "Human macrophage LPS immune activation (note: manual also mentions GSE19420)", species = "hs"),
    cs6 = list(geo = "GSE229613", title = "Curcumin treatment in HCT116 colon cancer cells", species = "hs")
)

load_geo_dataset <- function(gse_id, cache_dir) {
    suppressPackageStartupMessages(library(GEOquery))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    supp_dir <- file.path(cache_dir, gse_id)
    if (!dir.exists(supp_dir)) {
        tryCatch(
            GEOquery::getGEOSuppFiles(gse_id, baseDir = cache_dir, makeDirectory = TRUE),
            error = function(e) message("[WARN] Supplementary download note: ", e$message)
        )
    }
    supp_files <- list.files(supp_dir, pattern = "\\.(txt|tsv|csv)(\\.gz)?$", full.names = TRUE, recursive = TRUE)
    if (length(supp_files) > 0) {
        count_file <- supp_files[1]
        sep <- if (grepl("\\.csv", count_file, ignore.case = TRUE)) "," else "\t"
        df <- read.table(count_file, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
        rownames(df) <- make.unique(as.character(df[[1]]))
        num_cols <- sapply(df[, -1, drop = FALSE], is.numeric)
        counts <- as.matrix(round(df[, -1, drop = FALSE][, num_cols, drop = FALSE]))
        gse <- GEOquery::getGEO(gse_id, destdir = cache_dir, GSEMatrix = TRUE)
        pdata <- Biobase::pData(gse[[1]])
        cond <- if ("source_name_ch1" %in% colnames(pdata)) pdata$source_name_ch1 else rep(c("Control", "Treated"), length.out = ncol(counts))
        meta <- data.frame(sample = colnames(counts), condition = make.names(cond[seq_len(ncol(counts))]), stringsAsFactors = FALSE)
        return(list(counts = counts, metadata = meta))
    }
    stop(sprintf("No tabular supplementary count matrix automatically detected for %s in %s; inspect downloaded files and pass -c/-m explicitly.", gse_id, supp_dir))
}

build_synthetic_data <- function() {
    set.seed(42)
    n_genes <- 400
    n_samples <- 6
    genes <- sprintf("GENE_%04d", seq_len(n_genes))
    base_means <- round(rexp(n_genes, rate = 1 / 250)) + 20
    counts <- matrix(rpois(n_genes * n_samples, lambda = rep(base_means, n_samples)), nrow = n_genes)
    # Inject strong biological signal into first 30 genes (15 up, 15 down)
    counts[1:15, 4:6] <- counts[1:15, 4:6] * 5L
    counts[16:30, 1:3] <- counts[16:30, 1:3] * 5L
    rownames(counts) <- genes
    colnames(counts) <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Treat_1", "Treat_2", "Treat_3")
    meta <- data.frame(
        sample = colnames(counts),
        condition = factor(c("Control", "Control", "Control", "Treated", "Treated", "Treated"), levels = c("Control", "Treated")),
        stringsAsFactors = FALSE
    )
    list(counts = counts, metadata = meta)
}

read_count_and_meta <- function(counts_path, meta_path) {
    sep_c <- if (grepl("\\.csv$", counts_path, ignore.case = TRUE)) "," else "\t"
    raw <- read.table(counts_path, header = TRUE, sep = sep_c, check.names = FALSE, stringsAsFactors = FALSE)
    rownames(raw) <- make.unique(as.character(raw[[1]]))
    counts <- as.matrix(raw[, -1, drop = FALSE])
    storage.mode(counts) <- "integer"

    sep_m <- if (grepl("\\.csv$", meta_path, ignore.case = TRUE)) "," else "\t"
    meta <- read.table(meta_path, header = TRUE, sep = sep_m, stringsAsFactors = FALSE)
    if (!all(c("sample", "condition") %in% colnames(meta))) {
        colnames(meta)[1:2] <- c("sample", "condition")
    }
    counts <- counts[, meta$sample, drop = FALSE]
    meta$condition <- factor(make.names(meta$condition))
    list(counts = counts, metadata = meta)
}

run_dge <- function(counts, meta, method, lfc_cut, padj_cut) {
    counts <- counts[rowSums(counts >= 10) >= 2, , drop = FALSE]
    if (method == "deseq2") {
        suppressPackageStartupMessages(library(DESeq2))
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
        dds <- DESeq2::DESeq(dds, quiet = TRUE)
        res <- as.data.frame(DESeq2::results(dds))
        res$gene_id <- rownames(res)
        norm_mat <- DESeq2::counts(dds, normalized = TRUE)
        return(list(res = res, norm = norm_mat))
    } else {
        suppressPackageStartupMessages(library(edgeR))
        y <- edgeR::DGEList(counts = counts, group = meta$condition)
        y <- edgeR::calcNormFactors(y)
        design <- model.matrix(~ meta$condition)
        y <- edgeR::estimateDisp(y, design)
        fit <- edgeR::glmQLFit(y, design)
        qlf <- edgeR::glmQLFTest(fit, coef = 2)
        tab <- edgeR::topTags(qlf, n = Inf)$table
        res <- data.frame(
            gene_id = rownames(tab),
            baseMean = 2^tab$logCPM,
            log2FoldChange = tab$logFC,
            pvalue = tab$PValue,
            padj = tab$FDR,
            stringsAsFactors = FALSE
        )
        rownames(res) <- res$gene_id
        norm_mat <- edgeR::cpm(y, log = TRUE)
        return(list(res = res, norm = norm_mat))
    }
}

run_go_enrichment <- function(sig_genes, species, out_csv, out_png) {
    org_map <- list(hs = "org.Hs.eg.db", mm = "org.Mm.eg.db", ec = "org.EcK12.eg.db", at = "org.At.tair.db")
    pkg <- org_map[[species]]
    if (is.null(pkg) || !requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace(pkg, quietly = TRUE)) {
        return(invisible(NULL))
    }
    suppressPackageStartupMessages({
        library(clusterProfiler)
        library(enrichplot)
        library(pkg, character.only = TRUE)
    })
    key_type <- if (all(grepl("^ENS", sig_genes))) "ENSEMBL" else "SYMBOL"
    ego <- tryCatch(
        clusterProfiler::enrichGO(gene = sig_genes, OrgDb = get(pkg), keyType = key_type, ont = "BP", pAdjustMethod = "BH", readable = TRUE),
        error = function(e) NULL
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        write.csv(as.data.frame(ego), out_csv, row.names = FALSE)
        p <- enrichplot::dotplot(ego, showCategory = 15)
        ggplot2::ggsave(out_png, p, width = 8, height = 6, dpi = 300)
    }
}

write_markdown_top10 <- function(res_df, run_id, method, lfc_cut, padj_cut, out_md) {
    clean <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
    sig <- clean[clean$padj <= padj_cut & abs(clean$log2FoldChange) >= lfc_cut, ]
    up <- head(sig[sig$log2FoldChange > 0, ][order(-sig$log2FoldChange[sig$log2FoldChange > 0]), ], 10)
    down <- head(sig[sig$log2FoldChange < 0, ][order(sig$log2FoldChange[sig$log2FoldChange < 0]), ], 10)

    fmt_rows <- function(df, dir_label) {
        if (nrow(df) == 0) return(sprintf("| *(none)* | %s | - | - | - |", dir_label))
        vapply(seq_len(nrow(df)), function(i) {
            sprintf("| `%s` | %s | %.3f | %.2f | %.2e |", df$gene_id[i], dir_label, df$log2FoldChange[i], df$baseMean[i], df$padj[i])
        }, character(1))
    }

    lines <- c(
        sprintf("# Differential Gene Expression Summary: `%s` (%s)", run_id, toupper(method)),
        "",
        sprintf("- **Thresholds:** `|log2FoldChange| >= %.2f`, `padj <= %.3g`", lfc_cut, padj_cut),
        sprintf("- **Total Tested Genes:** %d | **Significant DEGs:** %d (Up: %d, Down: %d)",
                nrow(clean), nrow(sig), sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)),
        "",
        "## Top 10 Upregulated & Top 10 Downregulated Genes",
        "",
        "| Gene ID | Direction | log2FoldChange | BaseMean | Adjusted P-value |",
        "| :--- | :--- | ---: | ---: | ---: |",
        fmt_rows(up, "Upregulated"),
        fmt_rows(down, "Downregulated")
    )
    writeLines(lines, out_md)
}

generate_plots <- function(res_df, norm_mat, meta, lfc_cut, padj_cut, fig_dir, run_id) {
    suppressPackageStartupMessages({
        library(ggplot2)
        library(pheatmap)
    })
    clean <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
    clean$status <- ifelse(
        clean$padj <= padj_cut & clean$log2FoldChange >= lfc_cut, "Up",
        ifelse(clean$padj <= padj_cut & clean$log2FoldChange <= -lfc_cut, "Down", "NS")
    )
    clean$neg_log10_padj <- -log10(pmax(clean$padj, 1e-300))

    p_vol <- ggplot(clean, aes(x = log2FoldChange, y = neg_log10_padj, color = status)) +
        geom_point(alpha = 0.7, size = 1.6) +
        scale_color_manual(values = c(Up = "#D73027", Down = "#4575B4", NS = "grey70")) +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey40") +
        geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey40") +
        theme_minimal(base_size = 12) +
        labs(title = sprintf("Volcano Plot: %s", run_id), x = "log2 Fold Change", y = "-log10(adjusted p-value)")
    ggsave(file.path(fig_dir, sprintf("%s_volcano.png", run_id)), p_vol, width = 7.5, height = 5.5, dpi = 300)

    log_norm <- log2(norm_mat + 1)
    pca <- prcomp(t(log_norm), scale. = FALSE)
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], condition = meta$condition, sample = colnames(norm_mat))
    pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample)) +
        geom_point(size = 3.5) +
        geom_text(vjust = -0.8, size = 3) +
        theme_minimal(base_size = 12) +
        labs(title = sprintf("PCA Plot: %s", run_id), x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2]))
    ggsave(file.path(fig_dir, sprintf("%s_pca.png", run_id)), p_pca, width = 7, height = 5, dpi = 300)

    top_genes <- head(clean$gene_id[order(clean$padj)], 20)
    if (length(top_genes) >= 2) {
        ann_col <- data.frame(Condition = meta$condition, row.names = meta$sample)
        pheatmap::pheatmap(
            log_norm[top_genes, , drop = FALSE],
            scale = "row",
            annotation_col = ann_col,
            main = sprintf("Top 20 DEGs Heatmap: %s", run_id),
            filename = file.path(fig_dir, sprintf("%s_top-degs-heatmap.png", run_id)),
            width = 7,
            height = 7
        )
    }
}

opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))

run_id <- if (opts$synthetic_demo) {
    "synthetic-demo"
} else if (!is.null(opts$case_id)) {
    CASE_REGISTRY[[opts$case_id]]$geo
} else if (!is.null(opts$geo)) {
    opts$geo
} else {
    tools::file_path_sans_ext(basename(opts$counts))
}

if (!is.null(opts$case_id) && opts$species == "none") {
    opts$species <- CASE_REGISTRY[[opts$case_id]]$species
}

if (opts$dry_run) {
    cat(sprintf("[DRY-RUN] Run ID      : %s\n", run_id))
    cat(sprintf("[DRY-RUN] Method      : %s (|log2FC| >= %.2f, padj <= %.3g)\n", opts$method, opts$lfc, opts$padj))
    cat(sprintf("[DRY-RUN] Species DB  : %s\n", opts$species))
    cat(sprintf("[DRY-RUN] Output dir  : %s\n", opts$output_dir))
    quit(status = 0)
}

ngs_lib <- "/home/ht2905/miniforge3/envs/ngs/lib/R/library"
if (dir.exists(ngs_lib)) {
    .libPaths(unique(c(ngs_lib, .libPaths())))
}

table_dir <- file.path(opts$output_dir, "tables")
fig_dir <- file.path(opts$output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

dataset <- if (opts$synthetic_demo) {
    build_synthetic_data()
} else if (!is.null(opts$counts)) {
    if (is.null(opts$metadata)) {
        stop("Error: --metadata PATH is required when passing --counts.")
    }
    read_count_and_meta(opts$counts, opts$metadata)
} else {
    load_geo_dataset(run_id, file.path(opts$output_dir, "geo_cache"))
}

dge_out <- run_dge(dataset$counts, dataset$metadata, opts$method, opts$lfc, opts$padj)
res_df <- dge_out$res
res_sorted <- res_df[order(res_df$padj), ]

full_csv <- file.path(table_dir, sprintf("%s_dge-full.csv", run_id))
top10_md <- file.path(table_dir, sprintf("%s_top10-degs.md", run_id))
write.csv(res_sorted, full_csv, row.names = FALSE)
write_markdown_top10(res_sorted, run_id, opts$method, opts$lfc, opts$padj, top10_md)
generate_plots(res_sorted, dge_out$norm, dataset$metadata, opts$lfc, opts$padj, fig_dir, run_id)

sig_genes <- res_sorted$gene_id[!is.na(res_sorted$padj) & res_sorted$padj <= opts$padj & abs(res_sorted$log2FoldChange) >= opts$lfc]
run_go_enrichment(
    sig_genes,
    opts$species,
    file.path(table_dir, sprintf("%s_go-enrichment.csv", run_id)),
    file.path(fig_dir, sprintf("%s_go-dotplot.png", run_id))
)

message(sprintf("[OK] DGE complete: %s and %s", full_csv, top10_md))
