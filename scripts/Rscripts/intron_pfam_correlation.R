#!/usr/bin/env Rscript
# intron_pfam_correlation.R
#
# Tests for correlations between intron length and Pfam domain
# presence/absence or copy number across the Fungi5K dataset.
#
# Two complementary analyses per stratum:
#   1. Species-level: Spearman correlation between species mean intron length
#      and per-species Pfam domain copy number (and Wilcoxon for presence/absence).
#   2. Gene-level: Wilcoxon rank-sum test comparing per-gene mean intron length
#      between proteins with vs. without each Pfam domain.
#
# Usage:
#   Rscript intron_pfam_correlation.R [options]
#
#   --stratify none       Run on all taxa together (default)
#   --stratify phylum     Repeat analysis within each phylum
#   --stratify subphylum  Repeat analysis within each subphylum
#   --taxon <name>        Restrict to a single phylum/subphylum value
#                         (requires --stratify phylum or subphylum)
#
# Outputs land in:
#   results/pfam_intron_corr/<stratum>/  and  plots/pfam_intron_corr/<stratum>/
# where <stratum> is "all", "phylum/Ascomycota", "subphylum/Pezizomycotina", etc.

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(duckdb)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
  library(ggrepel)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-s", "--stratify"), type = "character", default = "none",
              metavar = "LEVEL",
              help = "Stratification level: none, phylum, or subphylum [default: %default]"),
  make_option(c("-t", "--taxon"), type = "character", default = NULL,
              metavar = "NAME",
              help = "Restrict to one taxon value within the chosen stratum level"),
  make_option(c("--db-dir"), type = "character", default = "functionalDB",
              metavar = "DIR",
              help = "Directory containing function.duckdb [default: %default]"),
  make_option(c("--db-name"), type = "character", default = "function.duckdb",
              metavar = "FILE",
              help = "DuckDB database filename [default: %default]"),
  make_option(c("--evalue"), type = "double", default = 1e-5,
              metavar = "FLOAT",
              help = "Pfam domain_i_evalue cutoff [default: %default]"),
  make_option(c("--min-species"), type = "integer", default = 50,
              metavar = "N",
              help = "Min species carrying a domain to test (all-taxa run) [default: %default]"),
  make_option(c("--min-species-strat"), type = "integer", default = 10,
              metavar = "N",
              help = "Min species carrying a domain when stratified [default: %default]"),
  make_option(c("--min-genes"), type = "integer", default = 200,
              metavar = "N",
              help = "Min genes with a domain to test (all-taxa gene-level run) [default: %default]"),
  make_option(c("--min-genes-strat"), type = "integer", default = 30,
              metavar = "N",
              help = "Min genes with a domain when stratified [default: %default]"),
  make_option(c("--fdr"), type = "double", default = 0.05,
              metavar = "FLOAT",
              help = "BH-FDR significance threshold [default: %default]"),
  make_option(c("--top-n"), type = "integer", default = 20,
              metavar = "N",
              help = "Number of top hits to plot in detail [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

STRATIFY_BY   <- tolower(opt$stratify)
TAXON_FILTER  <- opt$taxon
DBDIR         <- opt[["db-dir"]]
DBNAME        <- opt[["db-name"]]
EVALUE_CUTOFF <- opt$evalue
MIN_SPECIES   <- opt[["min-species"]]
MIN_SP_STRAT  <- opt[["min-species-strat"]]
MIN_GENES     <- opt[["min-genes"]]
MIN_GENES_STRAT <- opt[["min-genes-strat"]]
FDR_THRESHOLD <- opt$fdr
TOP_N_HITS    <- opt[["top-n"]]

if (!STRATIFY_BY %in% c("none", "phylum", "subphylum")) {
  stop("--stratify must be one of: none, phylum, subphylum")
}
if (!is.null(TAXON_FILTER) && STRATIFY_BY == "none") {
  stop("--taxon requires --stratify phylum or subphylum")
}

BASE_PLOTDIR    <- file.path("plots",   "pfam_intron_corr")
BASE_RESULTSDIR <- file.path("results", "pfam_intron_corr")

# ── Connect ───────────────────────────────────────────────────────────────────
con <- dbConnect(duckdb(), dbdir = file.path(DBDIR, DBNAME), read_only = TRUE)
message("Connected to ", file.path(DBDIR, DBNAME))

# ── Fetch all base data ───────────────────────────────────────────────────────
message("Fetching base data ...")

# Taxonomy: one row per locustag
taxonomy <- dbGetQuery(con, "
  SELECT LOCUSTAG, PHYLUM, SUBPHYLUM
  FROM species")
message(sprintf("  %d species in taxonomy table", nrow(taxonomy)))

# Mean/median intron length per species
intron_sp_all <- dbGetQuery(con, '
  SELECT locustag,
         avg("end" - "start" + 1)    AS mean_intron_len,
         median("end" - "start" + 1) AS median_intron_len,
         count(*)                AS n_introns
  FROM gene_introns
  GROUP BY locustag
  HAVING count(*) >= 10') %>%
  left_join(taxonomy, by = c("locustag" = "LOCUSTAG"))
message(sprintf("  %d species with intron data", nrow(intron_sp_all)))

# Pfam copy number per species
pfam_sp_all <- dbGetQuery(con, sprintf("
  SELECT species_prefix AS locustag,
         pfam_id,
         count(DISTINCT protein_id) AS copy_number
  FROM pfam
  WHERE domain_i_evalue < %g
  GROUP BY species_prefix, pfam_id", EVALUE_CUTOFF))
message(sprintf("  %d species × domain combinations", nrow(pfam_sp_all)))

# Mean intron length per transcript
intron_tx_all <- dbGetQuery(con, '
  SELECT transcript_id,
         locustag,
         avg("end" - "start" + 1)    AS mean_intron_len,
         count(*)                AS n_introns
  FROM gene_introns
  GROUP BY transcript_id, locustag') %>%
  left_join(taxonomy, by = c("locustag" = "LOCUSTAG"))
message(sprintf("  %d transcripts with intron data", nrow(intron_tx_all)))

# Pfam domains per transcript (via gene_proteins join)
pfam_tx_all <- dbGetQuery(con, sprintf("
  SELECT DISTINCT gp.transcript_id, p.pfam_id
  FROM gene_proteins AS gp
  JOIN pfam AS p ON gp.transcript_id = p.protein_id
  WHERE p.domain_i_evalue < %g", EVALUE_CUTOFF))
message(sprintf("  %d transcript × domain pairs", nrow(pfam_tx_all)))

dbDisconnect(con, shutdown = TRUE)
message("Disconnected.\n")

# ═══════════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════════

make_volcano <- function(results, title, fdr_thresh) {
  if (nrow(results) == 0) return(NULL)
  has_sig <- any(results$wilcox_fdr < fdr_thresh, na.rm = TRUE)
  top_rows <- results %>% slice_min(wilcox_fdr, n = 20, with_ties = FALSE)
  dat <- results %>%
    mutate(
      sig       = wilcox_fdr < fdr_thresh,
      log10p    = -log10(wilcox_p + 1e-300),
      label_txt = if_else(pfam_id %in% top_rows$pfam_id, pfam_id, NA_character_)
    )
  hline_y <- if (has_sig)
    -log10(max(results$wilcox_p[results$wilcox_fdr < fdr_thresh], na.rm = TRUE))
  else NA_real_

  p <- ggplot(dat, aes(x = rb_r, y = log10p, color = sig)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_label_repel(aes(label = label_txt), size = 2.4, na.rm = TRUE,
                     max.overlaps = 20, box.padding = 0.3) +
    { if (!is.na(hline_y))
        geom_hline(yintercept = hline_y, linetype = "dashed", color = "grey40") } +
    scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
      labels = c("not significant", sprintf("FDR < %g", fdr_thresh)),
      name = NULL) +
    xlab("Rank-biserial r  (positive = longer introns in group with domain)") +
    ylab(expression(-log[10](p))) +
    ggtitle(title) +
    theme_cowplot(12) +
    theme(legend.position = "bottom")
  p
}

# ── Species-level analysis ────────────────────────────────────────────────────
run_species_analysis <- function(intron_sp, pfam_sp, min_sp,
                                  label, plotdir, resultsdir) {
  dir.create(plotdir,    showWarnings = FALSE, recursive = TRUE)
  dir.create(resultsdir, showWarnings = FALSE, recursive = TRUE)

  # Domains present in enough of these species
  domain_freq <- pfam_sp %>%
    filter(locustag %in% intron_sp$locustag) %>%
    count(pfam_id, name = "n_species") %>%
    filter(n_species >= min_sp)

  if (nrow(domain_freq) == 0) {
    message(sprintf("    [%s] No domains pass min_species=%d – skipping species level", label, min_sp))
    return(invisible(NULL))
  }
  message(sprintf("    [%s] Testing %d domains at species level (min_species=%d)",
                  label, nrow(domain_freq), min_sp))

  pfam_sp_filt <- pfam_sp %>% semi_join(domain_freq, by = "pfam_id")

  run_one_domain <- function(domain) {
    with_dom <- pfam_sp_filt %>%
      filter(pfam_id == domain) %>%
      select(locustag, copy_number)

    joined <- intron_sp %>%
      left_join(with_dom, by = "locustag") %>%
      mutate(has_domain  = !is.na(copy_number),
             copy_number = replace_na(copy_number, 0))

    n_with    <- sum(joined$has_domain)
    n_without <- sum(!joined$has_domain)
    if (n_with < 5 || n_without < 5) return(NULL)

    ct <- suppressWarnings(
      cor.test(joined$mean_intron_len, joined$copy_number,
               method = "spearman", exact = FALSE))
    wt <- wilcox.test(
      joined$mean_intron_len[joined$has_domain],
      joined$mean_intron_len[!joined$has_domain],
      exact = FALSE)
    rb_r <- as.numeric(1 - 2 * wt$statistic / (n_with * n_without))

    tibble(
      pfam_id               = domain,
      n_species_with        = n_with,
      n_species_total       = nrow(joined),
      spearman_rho          = ct$estimate,
      spearman_p            = ct$p.value,
      wilcox_p              = wt$p.value,
      rb_r                  = rb_r,
      median_intron_with    = median(joined$mean_intron_len[joined$has_domain]),
      median_intron_without = median(joined$mean_intron_len[!joined$has_domain])
    )
  }

  results <- map_dfr(domain_freq$pfam_id, run_one_domain) %>%
    mutate(
      spearman_fdr = p.adjust(spearman_p, method = "BH"),
      wilcox_fdr   = p.adjust(wilcox_p,   method = "BH")
    ) %>%
    arrange(wilcox_fdr, spearman_fdr)

  write_tsv(results, file.path(resultsdir, "species_level_results.tsv"))
  n_sig <- sum(results$wilcox_fdr < FDR_THRESHOLD, na.rm = TRUE)
  message(sprintf("    [%s] Species-level significant (FDR<%.2f): %d",
                  label, FDR_THRESHOLD, n_sig))

  # Volcano
  p_vol <- make_volcano(results,
                        sprintf("Species-level: intron length ~ Pfam  [%s]", label),
                        FDR_THRESHOLD)
  if (!is.null(p_vol))
    ggsave(file.path(plotdir, "volcano_species.pdf"), p_vol, width = 10, height = 8)

  # Top-hits violin plots
  top_domains <- results %>%
    filter(wilcox_fdr < FDR_THRESHOLD) %>%
    slice_max(abs(rb_r), n = min(TOP_N_HITS, n_sig), with_ties = FALSE) %>%
    pull(pfam_id)

  if (length(top_domains) > 0) {
    pdf(file.path(plotdir, "top_hits_species.pdf"), width = 6, height = 5)
    for (dom in top_domains) {
      dat <- intron_sp %>%
        left_join(pfam_sp_filt %>% filter(pfam_id == dom), by = "locustag") %>%
        mutate(has_domain = factor(
          if_else(!is.na(copy_number), "Present", "Absent"),
          levels = c("Absent", "Present")))
      stat_row <- results %>% filter(pfam_id == dom)
      subtitle  <- sprintf("FDR=%.2e | rb_r=%.3f | n_with=%d",
                           stat_row$wilcox_fdr, stat_row$rb_r, stat_row$n_species_with)
      p <- ggplot(dat, aes(x = has_domain, y = mean_intron_len, fill = has_domain)) +
        geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
        geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +
        scale_fill_manual(values = c("Absent" = "#4575b4", "Present" = "#d73027")) +
        xlab(dom) + ylab("Mean intron length (bp)") +
        ggtitle(sprintf("[%s]  %s", label, dom), subtitle = subtitle) +
        theme_cowplot(10) + theme(legend.position = "none")
      print(p)
    }
    dev.off()
  }

  invisible(results)
}

# ── Gene-level analysis ───────────────────────────────────────────────────────
run_gene_analysis <- function(intron_tx, pfam_tx, min_genes,
                               label, plotdir, resultsdir) {
  dir.create(plotdir,    showWarnings = FALSE, recursive = TRUE)
  dir.create(resultsdir, showWarnings = FALSE, recursive = TRUE)

  # Restrict pfam_tx to transcripts in this stratum
  pfam_tx_sub <- pfam_tx %>%
    filter(transcript_id %in% intron_tx$transcript_id)

  domain_freq <- pfam_tx_sub %>%
    count(pfam_id, name = "n_genes") %>%
    filter(n_genes >= min_genes)

  if (nrow(domain_freq) == 0) {
    message(sprintf("    [%s] No domains pass min_genes=%d – skipping gene level", label, min_genes))
    return(invisible(NULL))
  }
  message(sprintf("    [%s] Testing %d domains at gene level (min_genes=%d)",
                  label, nrow(domain_freq), min_genes))

  pfam_tx_filt <- pfam_tx_sub %>% semi_join(domain_freq, by = "pfam_id")

  # Build index: transcript_id → has_domain per pfam_id tested later
  tx_ids <- intron_tx$transcript_id

  run_one_domain <- function(domain) {
    tx_with_dom  <- pfam_tx_filt %>% filter(pfam_id == domain) %>% pull(transcript_id)
    len_with     <- intron_tx$mean_intron_len[tx_ids %in% tx_with_dom]
    len_without  <- intron_tx$mean_intron_len[!(tx_ids %in% tx_with_dom)]
    n_with    <- length(len_with)
    n_without <- length(len_without)
    if (n_with < 10 || n_without < 10) return(NULL)

    wt   <- wilcox.test(len_with, len_without, exact = FALSE)
    rb_r <- as.numeric(1 - 2 * wt$statistic / (n_with * n_without))

    tibble(
      pfam_id               = domain,
      n_genes_with          = n_with,
      n_genes_without       = n_without,
      median_intron_with    = median(len_with),
      median_intron_without = median(len_without),
      log2FC_median         = log2((median(len_with) + 1) / (median(len_without) + 1)),
      wilcox_p              = wt$p.value,
      rb_r                  = rb_r
    )
  }

  results <- map_dfr(domain_freq$pfam_id, run_one_domain) %>%
    mutate(wilcox_fdr = p.adjust(wilcox_p, method = "BH")) %>%
    arrange(wilcox_fdr, desc(abs(rb_r)))

  write_tsv(results, file.path(resultsdir, "gene_level_results.tsv"))
  n_sig <- sum(results$wilcox_fdr < FDR_THRESHOLD, na.rm = TRUE)
  message(sprintf("    [%s] Gene-level significant (FDR<%.2f): %d",
                  label, FDR_THRESHOLD, n_sig))

  # Volcano
  p_vol <- make_volcano(results,
                        sprintf("Gene-level: intron length ~ Pfam  [%s]", label),
                        FDR_THRESHOLD)
  if (!is.null(p_vol))
    ggsave(file.path(plotdir, "volcano_gene.pdf"), p_vol, width = 10, height = 8)

  # Top-hits violin plots
  top_domains <- results %>%
    filter(wilcox_fdr < FDR_THRESHOLD) %>%
    slice_max(abs(rb_r), n = min(TOP_N_HITS, n_sig), with_ties = FALSE) %>%
    pull(pfam_id)

  if (length(top_domains) > 0) {
    pdf(file.path(plotdir, "top_hits_gene.pdf"), width = 6, height = 5)
    for (dom in top_domains) {
      tx_with_dom <- pfam_tx_filt %>% filter(pfam_id == dom) %>% pull(transcript_id)
      dat <- intron_tx %>%
        mutate(has_domain = factor(
          if_else(transcript_id %in% tx_with_dom, "Present", "Absent"),
          levels = c("Absent", "Present")))
      stat_row <- results %>% filter(pfam_id == dom)
      subtitle  <- sprintf("FDR=%.2e | rb_r=%.3f | n_with=%d",
                           stat_row$wilcox_fdr, stat_row$rb_r, stat_row$n_genes_with)
      p <- ggplot(dat, aes(x = has_domain, y = mean_intron_len, fill = has_domain)) +
        geom_violin(alpha = 0.7, draw_quantiles = c(0.25, 0.5, 0.75)) +
        scale_fill_manual(values = c("Absent" = "#4575b4", "Present" = "#d73027")) +
        scale_y_log10() +
        xlab(dom) + ylab("Mean intron length (bp, log scale)") +
        ggtitle(sprintf("[%s]  %s", label, dom), subtitle = subtitle) +
        theme_cowplot(10) + theme(legend.position = "none")
      print(p)
    }
    dev.off()
  }

  invisible(results)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Build the list of strata to run
# ═══════════════════════════════════════════════════════════════════════════════

# Each element of `strata` is a list(label, locustags, min_sp, min_genes, subdir)
strata <- list()

if (STRATIFY_BY == "none") {
  strata[[1]] <- list(
    label    = "all",
    locustags = intron_sp_all$locustag,
    min_sp   = MIN_SPECIES,
    min_genes = MIN_GENES,
    subdir   = "all"
  )
} else {
  col <- if (STRATIFY_BY == "phylum") "PHYLUM" else "SUBPHYLUM"

  taxon_values <- taxonomy %>%
    filter(!is.na(.data[[col]]), .data[[col]] != "") %>%
    pull(.data[[col]]) %>%
    unique() %>%
    sort()

  if (!is.null(TAXON_FILTER)) {
    taxon_values <- taxon_values[taxon_values == TAXON_FILTER]
    if (length(taxon_values) == 0)
      stop(sprintf("--taxon '%s' not found in %s column", TAXON_FILTER, col))
  }

  for (tv in taxon_values) {
    loctags <- taxonomy %>%
      filter(.data[[col]] == tv) %>%
      pull(LOCUSTAG)
    strata[[length(strata) + 1]] <- list(
      label     = tv,
      locustags = loctags,
      min_sp    = MIN_SP_STRAT,
      min_genes = MIN_GENES_STRAT,
      subdir    = file.path(STRATIFY_BY, tv)
    )
  }
}

message(sprintf("\nRunning %d stratum/strata  (stratify=%s)\n", length(strata), STRATIFY_BY))

# ═══════════════════════════════════════════════════════════════════════════════
# Main loop over strata
# ═══════════════════════════════════════════════════════════════════════════════

summary_rows <- list()

for (st in strata) {
  message(sprintf("── Stratum: %s  (%d locustags) ──", st$label, length(st$locustags)))

  intron_sp_sub <- intron_sp_all %>% filter(locustag %in% st$locustags)
  pfam_sp_sub   <- pfam_sp_all   %>% filter(locustag %in% st$locustags)
  intron_tx_sub <- intron_tx_all %>% filter(locustag %in% st$locustags)

  if (nrow(intron_sp_sub) < st$min_sp) {
    message(sprintf("  Skipping: only %d species with intron data (need %d)",
                    nrow(intron_sp_sub), st$min_sp))
    next
  }

  pdir <- file.path(BASE_PLOTDIR,    st$subdir)
  rdir <- file.path(BASE_RESULTSDIR, st$subdir)

  sp_res   <- run_species_analysis(intron_sp_sub, pfam_sp_sub,
                                   st$min_sp, st$label, pdir, rdir)
  gene_res <- run_gene_analysis(intron_tx_sub, pfam_tx_all,
                                st$min_genes, st$label, pdir, rdir)

  # Collect summary
  n_sp_sig   <- if (!is.null(sp_res))
    sum(sp_res$wilcox_fdr   < FDR_THRESHOLD, na.rm = TRUE) else NA_integer_
  n_gene_sig <- if (!is.null(gene_res))
    sum(gene_res$wilcox_fdr < FDR_THRESHOLD, na.rm = TRUE) else NA_integer_

  summary_rows[[length(summary_rows) + 1]] <- tibble(
    stratum           = st$label,
    n_species         = nrow(intron_sp_sub),
    n_transcripts     = nrow(intron_tx_sub),
    sp_sig_domains    = n_sp_sig,
    gene_sig_domains  = n_gene_sig
  )
}

# ── Cross-stratum summary ─────────────────────────────────────────────────────
if (length(summary_rows) > 0) {
  summary_tbl <- bind_rows(summary_rows)
  write_tsv(summary_tbl, file.path(BASE_RESULTSDIR, "stratum_summary.tsv"))
  message("\n── Stratum summary ──")
  print(summary_tbl, n = Inf)
}

message(sprintf("\nDone. Results in %s/  |  Plots in %s/",
                BASE_RESULTSDIR, BASE_PLOTDIR))
