#' load packages for RNAseq
#' @export
load_rnaseq <- function() {
  pacman::p_load(
    'tidyverse',
    'biomaRt',
    'tximport',
    'edgeR',
    'limma',
    'DESeq2',
    'glmGamPoi',
    'Glimma',
    'patchwork',
    'plotly',
    'DT',
    'gplots',
    'RColorBrewer',
    'Biobase',
    'GSEABase',
    'GSVA',
    'gprofiler2',
    'clusterProfiler',
    'msigdbr',
    'enrichplot',
    'CAMERA',
    'ggokabeito'
  )
}


#' PCA PVE Plot
#' @export
plot_pca_pve <- function(pc_pve) {
  par(mfrow = c(1, 2))
  plot(
    pc_pve, ylim = c(0,1), type = "b", col = 2,
    ylab = "Proportion of Variance Explained (PVE)", xlab = "Principal Component"
  )
  plot(
    cumsum(pc_pve), ylim = c(0,1), type = "b", col = 4,
    ylab = "Cummulative PVE", xlab = "Principal Component"
  )
}


#' PCA Plot
#' @export
plot_pca <- function(pc_wide, pca_x, pca_y, pve_x, pve_y, label, shape, color){
  ggplot(pc_wide) +
    aes(
      x = get(pca_x), y = get(pca_y), 
      label = label,  
      shape = shape,
      color = color
    ) +
    geom_point(size = 4) +
    # geom_label() +
    # stat_ellipse() +
    xlab(paste0(pca_x, " (", pve_x * 100, "%", ")")) + 
    ylab(paste0(pca_y, " (", pve_y * 100, "%", ")")) +
    labs(
      title = "PCA plot",
      caption = paste0("CQWU ", Sys.time())
    ) +
    coord_fixed() +
    theme_classic()
}


#' PCA Small Multiples Plot
#' @export
plot_pca_multiple <- function(pc_wide, sample, group){
  pc_small <- pc_wide[, 1:4] %>%
    add_column(
      sample = sample_info$sample_id,
      group = factor(sample_info$treatment)                           ## INPUT
    )
  pc_small_tidy <- pivot_longer(
    pc_small,
    cols = PC1:PC4, 
    names_to = "PC",
    values_to = "loadings"
  )
  ggplot(pc_small_tidy) +              
    aes(
      x = sample, y = loadings,
      fill = group                                                    ## INPUT
    ) +
    geom_bar(stat = "identity") +
    facet_wrap(~PC) +
    labs(
      title = "PCA 'small multiples' plot",
      caption = paste0("CQWU ", Sys.time())
    ) +
    theme_bw() +
    coord_flip()
}


#' Convert tibble to matrix retaining rownames
#' @export
tibble_to_matrix <- function(x) {
  y <- as.matrix(x[, -1])
  rownames(y) <- x[[1]]
  return (y)
}


#' Add SRA meta data into sample_info
#' @export
add_SRA_meta <- function(sample_info) {
  rename(sample_info, seq_id = SRA_Accession)  ## seq_id is SRA Accession
  SRA_meta <- read_csv("../meta/SRA_meta.csv") %>%
    rename(seq_id = Run) %>%
    select(seq_id, Model)
  sample_info <- sample_info %>%
    left_join(SRA_meta, by = "seq_id")
}


#' Get human or mouse transcript to gene annotation with biomaRt
#' @export
get_t2g <- function(x) {
  library(tidyverse)
  gene_ensembl <- case_when(
    x %in% c("human", "Human") ~ "hsapiens_gene_ensembl",
    x %in% c("mouse", "Mouse") ~ "mmusculus_gene_ensembl"
  )
  t2g <- biomaRt::getBM(
    attributes = c(
      'ensembl_transcript_id',
      'external_gene_name',
      'description'
    ),
    mart = biomaRt::useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = gene_ensembl
    )
  ) %>%
    as_tibble() %>%
    rename(
      gene_id = ensembl_transcript_id,
      gene_name = external_gene_name
    ) %>% 
    dplyr::distinct(gene_id, gene_name, .keep_all = TRUE)
  write_rds(t2g, "t2g_biomaRt.rds")

  mart <- biomaRt::useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = gene_ensembl
  )
  write_tsv(biomaRt::listEnsembl(mart)[1, ], "t2g_biomaRt_version.txt")

  return(t2g)
}

#' Get human or mouse gene_ID to gene_name annotation with biomaRt
#' @export
get_g2g <- function(x) {
  library(tidyverse)
  gene_ensembl <- case_when(
    x %in% c("human", "Human") ~ "hsapiens_gene_ensembl",
    x %in% c("mouse", "Mouse") ~ "mmusculus_gene_ensembl"
  )
  g2g <- biomaRt::getBM(
    attributes = c(
      'ensembl_gene_id',
      'external_gene_name',
      'description'
    ),
    mart = biomaRt::useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = gene_ensembl
    )
  ) %>%
    as_tibble() %>%
    rename(
      gene_id = ensembl_gene_id,
      gene_name = external_gene_name
    ) %>% 
    dplyr::distinct(gene_id, gene_name, .keep_all = TRUE)
  write_rds(g2g, "g2g_biomaRt.rds")

  mart <- biomaRt::useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = gene_ensembl
  )
  write_tsv(biomaRt::listEnsembl(mart)[1, ], "g2g_biomaRt_version.txt")

  return(g2g)
}

#' Get human transcript to gene annotation with EnsDb.Hsapiens.v86
#' @export
get_t2g_human_v86 <- function() {
  library(ensembldb)
  library(EnsDb.Hsapiens.v86) 
  t2g <- ensembldb::transcripts(
    EnsDb.Hsapiens.v86,                                  
    columns = c("tx_id", "gene_name")
  ) %>%
    as_tibble() %>%
    dplyr::rename(target_id = tx_id) %>%
    dplyr::select(target_id, gene_name) %>% 
  detach("package:EnsDb.Hsapiens.v86", unload = TRUE)
  detach("package:ensembldb", unload = TRUE)  
  return(t2g)
}

#' Get mouse transcript to gene annotation with EnsDb.Mmusculus.v79
#' @export
get_t2g_mouse_v79 <- function() {
  library(ensembldb)
  library(EnsDb.Mmusculus.v79)
  t2g <- ensembldb::transcripts(
    EnsDb.Mmusculus.v79,                                  
    columns = c("tx_id", "gene_name")
  ) %>%
    as_tibble() %>%
    dplyr::rename(target_id = tx_id) %>%
    dplyr::select(target_id, gene_name) %>% 
  detach("package:EnsDb.Mmusculus.v79", unload = TRUE)
  detach("package:ensembldb", unload = TRUE)
  return(t2g)
}


#' Get cpm_log2_tidy 
#' @export
get_cpm_log2_tidy <- function (cpm_log2, sample_info) {
  cpm_log2_wide <- as_tibble(cpm_log2, rownames = "gene")
  colnames(cpm_log2_wide) <- c("gene", sample_info$sample_id)
  cpm_log2_tidy <- pivot_longer(
    cpm_log2_wide,
    cols = sample_info$sample_id,
    names_to = "sample_id",
    values_to = "log2_expr"
  )
  return(cpm_log2_tidy)
}


#' Make TPM table
#' @export
table_tpm <- function(tx_gene, sample_info, cpm_log2) {
  tpm <- tx_gene$abundance %>% as_tibble(rownames = "gene")
  colnames(tpm) <- c("gene", sample_info$sample_id)
  cpm_log2_wide <- as_tibble(cpm_log2, rownames = "gene")
  tpm <- tpm[tpm$gene %in% cpm_log2_wide$gene, ]
  write_tsv(tpm, "tpm.tsv")
  DT::datatable(
    tpm,
    extensions = c('KeyTable', "FixedHeader"), 
    filter = 'top',
    options = list(
      keys = TRUE, 
      searchHighlight = TRUE, 
      pageLength = 10, 
      lengthMenu = c("10", "25", "50", "100")
    )
  ) %>%
    formatRound(2 : (length(sample_info$sample_id) + 1), digits = 1) %>%
    saveWidget("tpm.html")
}


#' CPM Violin Plot
#' @export
plot_cpm_violin <- function(cpm_tidy, subtitle) {
  ggplot(cpm_tidy) +
  aes(x = sample_id, y = log2_expr, fill = sample_id) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(
    fun = "median", geom = "point", shape = 95, 
    size = 5, color = "black", show.legend = FALSE
  ) +
  labs(
    x = "sample_id", y = "log2_expr", 
    title = " Log2 counts per million (CPM)",
    subtitle = subtitle,
    caption = paste0("By CQWU on ", Sys.time())
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
}


#' DGE Volcano plot
#' @export
plot_volcano <- function(dge, subtitle, adj_p = 0.05) {
  ggplot(dge) +
    aes(x = logFC, y = -log10(adj.P.Val), text = paste("Symbol:", gene)) +
    geom_point(size = 2) +
    geom_hline(
      yintercept=-log10(adj_p),linetype="longdash",colour="grey",size=1
    ) +
    geom_vline(xintercept=1, linetype="longdash", colour="#BE684D", size=1) +
    geom_vline(xintercept=-1, linetype="longdash", colour="#2C467A", size=1) +
    #annotate(
    #  "rect",xmin=1,xmax=12,ymin=-log10(0.01),ymax=7.5,alpha=.2,fill="#BE684D"
    #) +
    #annotate(
    #  "rect",xmin=-1,xmax=-12,ymin=-log10(0.01),ymax=7.5,alpha=.2,fill="#2C467A"
    #) +
    labs(
      subtitle = subtitle,
      caption = paste0("by CQWU on ", Sys.time())) +
    theme_classic()
}


#' Join top hit
#' @export
table_top_hit_to_join <- function(top_hits, name){
  top_hits_to_join <- top_hits %>%  
    rename(
      !! paste0(str_replace(name, "top_hits_", ""), "_logFC") := logFC,
      !! paste0(str_replace(name, "top_hits_", ""), "_AveExpr") := AveExpr,        
      !! paste0(str_replace(name, "top_hits_", ""), "_P.Value") := P.Value,    
      !! paste0(str_replace(name, "top_hits_", ""), "_adj.P.Val") := adj.P.Val,
    )
  return (top_hits_to_join)
}

######################################################################## The END
