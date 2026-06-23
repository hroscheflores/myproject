#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GO.db)
  library(AnnotationDbi)
})

base_dir <- "path/to/project_directory"
out_dir <- file.path(base_dir, "go_enrichment_output")

gene2go_file <- file.path(base_dir, "gene2go.tsv")
genes_sf_file <- file.path(base_dir, "genes_by_superfamily.tsv")

target_superfamily_pattern <- "Sola"
comparison_label <- "other TE-associated genes"

out_full <- file.path(out_dir, "target_vs_other_superfamilies_GO_full.tsv")
out_slim <- file.path(out_dir, "target_vs_other_superfamilies_GO_slim.tsv")
out_pdf <- file.path(out_dir, "target_vs_other_GO_slim_dotplot.pdf")

min_total <- 30
min_target <- 3
top_n <- 25

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

gene2go <- read.delim(gene2go_file, stringsAsFactors = FALSE)
genes_sf <- read.delim(genes_sf_file, stringsAsFactors = FALSE)

required_gene2go_cols <- c("gene_id", "go_id")
required_genes_sf_cols <- c("gene_id", "superfamily")

stopifnot(all(required_gene2go_cols %in% colnames(gene2go)))
stopifnot(all(required_genes_sf_cols %in% colnames(genes_sf)))

gene2go <- gene2go %>%
  mutate(
    gene_id = as.character(gene_id),
    go_id = sub("^[PFC]:", "", go_id)
  ) %>%
  filter(go_id != "", !is.na(go_id)) %>%
  distinct(gene_id, go_id)

genes_sf <- genes_sf %>%
  mutate(gene_id = as.character(gene_id))

genes_te <- intersect(unique(genes_sf$gene_id), unique(gene2go$gene_id))

gene2go <- gene2go %>%
  filter(gene_id %in% genes_te)

genes_sf <- genes_sf %>%
  filter(gene_id %in% genes_te)

target_genes <- unique(
  genes_sf$gene_id[
    grepl(target_superfamily_pattern, genes_sf$superfamily, ignore.case = TRUE)
  ]
)

other_genes <- setdiff(genes_te, target_genes)

go_counts <- gene2go %>%
  distinct(gene_id, go_id) %>%
  group_by(go_id) %>%
  summarise(
    target_with = sum(gene_id %in% target_genes),
    other_with = sum(gene_id %in% other_genes),
    total_with = n(),
    .groups = "drop"
  )

n_target <- length(target_genes)
n_other <- length(other_genes)

pvals <- numeric(nrow(go_counts))
odds <- numeric(nrow(go_counts))

for (i in seq_len(nrow(go_counts))) {
  a <- go_counts$target_with[i]
  b <- go_counts$other_with[i]
  c <- n_target - a
  d <- n_other - b
  
  ft <- suppressWarnings(fisher.test(matrix(c(a, b, c, d), nrow = 2)))
  
  pvals[i] <- ft$p.value
  odds[i] <- if (!is.null(ft$estimate)) as.numeric(ft$estimate) else NA_real_
}

go_counts <- go_counts %>%
  mutate(
    pvalue = pvals,
    OR = odds,
    padj = p.adjust(pvalue, method = "BH")
  )

term_info <- suppressMessages(
  AnnotationDbi::select(
    GO.db,
    keys = unique(go_counts$go_id),
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
  )
) %>%
  filter(!is.na(TERM)) %>%
  distinct(GOID, TERM, ONTOLOGY)

res_full <- go_counts %>%
  left_join(term_info, by = c("go_id" = "GOID")) %>%
  mutate(TERM = ifelse(is.na(TERM), go_id, TERM)) %>%
  arrange(pvalue)

write.table(
  res_full,
  file = out_full,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

res_slim <- res_full %>%
  filter(
    !is.na(ONTOLOGY),
    ONTOLOGY == "BP",
    total_with >= min_total,
    target_with >= min_target
  ) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(padj, pvalue)

write.table(
  res_slim,
  file = out_slim,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_df <- res_slim %>%
  mutate(
    log2_OR = log2(OR),
    log2_OR = ifelse(is.finite(log2_OR), log2_OR, NA_real_),
    TERM = ifelse(TERM == "" | is.na(TERM), go_id, TERM)
  ) %>%
  filter(!is.na(log2_OR)) %>%
  arrange(padj, pvalue) %>%
  slice_head(n = top_n) %>%
  mutate(TERM = factor(TERM, levels = rev(unique(TERM))))

p <- ggplot(
  plot_df,
  aes(x = log2_OR, y = TERM, size = target_with, color = padj)
) +
  geom_point() +
  scale_color_viridis_c(option = "D", direction = -1, name = "adj. p") +
  scale_size_continuous(name = "# target genes") +
  labs(
    title = paste(
      "GO terms enriched in",
      target_superfamily_pattern,
      "associated genes"
    ),
    subtitle = paste("Compared against", comparison_label),
    x = paste0("log2 odds ratio (", target_superfamily_pattern, " vs other TE genes)"),
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

ggsave(out_pdf, p, width = 8, height = 6)