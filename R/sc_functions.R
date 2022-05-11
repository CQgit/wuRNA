#' Load packages
#' @export
load_sc <- function() {
  pacman::p_load(
    'BiocParallel',
    'tidyverse',
    'Matrix',
    'DropletUtils',
    'Seurat',  
    'sctransform',
    'glmGamPoi',  
    'patchwork',
    'scales',
    'rjson',
    'R2HTML',
    'DT',
    'viridis',
    'scater',
    'scran',
    'SingleR',
    'celldex',
    'pheatmap'
  )
}


#' Seurat QC Plot
#' @export
plot_qc_seurat <- function(neo, nCount_intercept = 20000, nFeature_intercept = 200) {
  p1 <- ggplot(neo@meta.data, aes(nCount_RNA, nFeature_RNA, color = pct_MT)) +
    geom_point(alpha = 0.6, size = 0.5) +
    #labs(x = "Total UMI counts per cell", y = "Number of genes detected") +
    scale_color_viridis()
  p2 <- ggplot(neo@meta.data, aes(nCount_RNA)) +
    geom_histogram(binwidth = 200) +
    geom_vline(xintercept = nCount_intercept, color = "Red", linetype = "dashed")
  p3 <- ggplot(neo@meta.data, aes(nFeature_RNA)) +
    geom_histogram(binwidth = 20) +
      geom_vline(xintercept = nFeature_intercept, color = "Red", linetype = "dashed")
  p4 <- VlnPlot(neo, features = c("nCount_RNA", "nFeature_RNA", "pct_MT"), pt.size = 0.01, ncol = 3)
  p1 / p2 / p3 / p4
}


#' plain 
#' @export
plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE, big.mark = ",")
}


#' Drop empty barcods and Plot
#' @export
drop_empty_barcodes <- function (path_cellranger) {
## -----------------------------------------------------------------------------
  bc_rank_plot <- function(stats, raw_cells, filt_cells, save){
    cells <- raw_cells %in% filt_cells
    keep <- !duplicated(stats$total)
    plot_df <- data.frame(Rx=stats$rank, Tx=stats$total, cell=cells)
    plot_df <- plot_df[keep, ]
    png(save, width = 700, height = 500)
    print({
      ggplot(subset(plot_df, plot_df$Tx >0), aes(x=Rx, y=Tx, col=cell, alpha=cell))+ 
        geom_point(size=3) + 
        geom_hline(yintercept = stats@metadata$knee, lty = 2, col = '#0972D5', size=1.5) +
        annotate("text", x=max(plot_df$Tx), y=stats@metadata$knee+10000, label="Knee", color = "#0972D5", size=5.5) +
        geom_hline(yintercept = stats@metadata$inflection, lty = 2, col = '#09CFD5', size = 1.5) +
        annotate("text", x=max(plot_df$Tx), y=stats@metadata$inflection+500, label="Inflection", color = "#09CFD5", size=5.5) +
        scale_x_log10(labels = plain, breaks = trans_breaks("log10", function(x) round(10^x, 0))) + 
        scale_y_log10(breaks = trans_breaks('log10', function(x) floor(10^x)), labels = plain)+
        scale_color_manual(values = c('#8595A8', '#6406B6'), name = NULL, labels = c("Background", "Cells")) +
        scale_alpha_manual(values = c(0.5,1)) +
        labs(x = 'Barcodes', y = 'UMI counts', title = 'Barcode Rank Plot') +
        guides(alpha = "none", colour = guide_legend(reverse = TRUE, override.aes=list(size = 5))) +
        theme_linedraw() + 
        theme(plot.title = element_text(size=20, hjust=0.5, face = 'bold'),
              axis.title = element_text(size = 18),
              axis.text = element_text(size = 15),
              legend.text = element_text(size=19),
              legend.background = element_rect(fill = 'transparent'),
              legend.position = c(0.15,0.15))
    })
    dev.off()
  }
  
  print_HTML <- function(seq_stats, cell_stats, dir, sample_id){
    system(paste0('base64 ', dir, '/barcode_rank.png > ', dir, '/barcode_rank.txt'))
    b64_bc <- readChar(paste0(dir, '/barcode_rank.txt'), file.info(paste0(dir, '/barcode_rank.txt'))$size)
    target <- HTMLInitFile(dir, filename=paste0(sample_id, '_summary'))
    HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Roboto">', file=target)
    HTML("<div class='title'>", file=target)
    HTML.title(' Pre-Processing Summary', HR=1, file = target)
    HTML("</div>", file = target)
    HTML.title(sample_id, HR=2, file = target)
    HTML("<div id='wrapper'>", file=target)
    HTML("<div class='boxed' id='left' align='center'>", file=target)
    HTML.title('Sequencing/Alignment Stats', HR=3, file=target)
    HTML('<table style="width:100%">', file=target)
    HTML(paste('<tr> <td>', seq_stats$stat, '</td> <td align="right">', seq_stats$value, '</td> </tr>'), file=target)
    HTML('</table> <hr>', file=target)
    HTML.title('Cell Stats', HR=3, file=target)
    HTML('<table style="width:100%">', file=target)
    HTML(paste('<tr> <td>', cell_stats$stat, '</td> <td align="right">', cell_stats$value, '</td> </tr>'), file=target)
    HTML('</table>', file=target)
    HTML("</div>", file = target)
    HTML("<div class='boxed' id='right' align='center'>", file=target)
    HTML(paste0("<img src='data:image/png;base64,", b64_bc, "' width=90%>"), file=target)
    HTML("</div>", file = target)
    HTML("</div>", file = target)
    HTML('<style type="text/css">
  		.title {
      			background-color: #0972D5;
      			padding: 8px;
      			color: white;
      			position: fixed;
      			top: 0;
      			left: 0;
      			z-index: 999;
      			width: 100%;
  		}
  		.boxed {
    			border: 1px solid #868D96;
    			padding: 10px;
    			margin: 20px;
  		}
  		h1 {
  			font-family: "Roboto";
  			font-size: 33px;
  		}
  		h2 {
  			font-family: "Roboto";
  			font-size: 26px;
  		}
  		h3 {
  			font-family: "Roboto";
  			font-size: 18px;
  		}
  		#wrapper {
    			display: flex;
  		}
  		#left {
    			width: 50%;
  		}
  		#right {
    			width: 50%;
  		}
  		table tr:nth-child(even) {
    			background-color: #eee;
  		}
  		table tr:nth-child(odd) {
    			background-color: #fff;
  		}
  		table {
    			font-family: "Roboto";
  			font-size: 20px;
  			border: 1px solid #868D96;
  		}
  		#mathplayer{
    			height: 80px;
  		}
  		</style> </head>', file=target)
    HTMLEndFile()
  }

## ---------------------------------------------------------- Read kb-cellranger  
  path_mtx        <- file.path(path_cellranger, "matrix.mtx")
  path_genes      <- file.path(path_cellranger, "genes.tsv")
  path_barcodes   <- file.path(path_cellranger, "barcodes.tsv")
  
  mtx <- Matrix::readMM(path_mtx)
  genes <- read_tsv(path_genes, col_names = FALSE) %>% as.data.frame()
  barcodes <- read_tsv(path_barcodes, col_names = FALSE) %>% as.data.frame()
  rownames(mtx) <- genes[, 1]
  colnames(mtx) <- barcodes[, 1]
  
## ---------------------------------------------------------- Drop empty barcode
  keep <- DropletUtils::emptyDrops(mtx, BPPARAM = MulticoreParam())$FDR <= 0.05
  keep <- replace_na(keep, FALSE)
  mtx_d <- mtx[, keep]
  DropletUtils::write10xCounts(
    "r/kb_cellranger_d", gene.symbol = genes[, 2], mtx_d, overwrite = TRUE
  ) 
  
## ------------------------------------------------------------------- seq_stats
  kb_stats <- c(
    fromJSON(file = "kallisto/inspect.json"),
    fromJSON(file = "kallisto/run_info.json")
  ) 
  chem_10X <- grep('10X(.*)', strsplit(kb_stats$call, '\\s')[[1]], value = TRUE)
  seq_stats <- data.frame(
    stat = c(
      'Sequencing technology', 'Number of reads processed', 
      '% reads pseudoaligned', '% reads on whitelist'
    ),
    value = prettyNum(
      c(
        chem_10X, kb_stats$n_processed, kb_stats$p_pseudoaligned, 
        round(kb_stats$percentageReadsOnWhitelist,2)
      ),
      big.mark = ','
    )
  )
  
## ------------------------------------------------------------------ cell_stats
  per_counts_cell <- round((sum(mtx_d)/sum(mtx)) * 100, 2) 
  med_counts <- median(colSums(mtx_d))
  med_genes <- median(apply(mtx_d, 2, function(x) sum(x >= 1)))
  total_genes <- sum(rowSums(mtx_d) >= 1)
  cell_stats <- data.frame(
    stat = c(
      'Estimated number of cells', '% counts in cells', 
      'Median counts per cell', 'Median genes per cell', 'Total genes detected'
    ), 
    value = prettyNum(
      c(
        ncol(mtx_d), per_counts_cell, med_counts, med_genes, total_genes
      ),
      big.mark = ','
    )
  )
  
## ----------------------------------------------------------- Barcode QA report
  stats <- DropletUtils::barcodeRanks(mtx)
  barcodes_d <- read_tsv("r/kb_cellranger_d/barcodes.tsv", col_names = FALSE) %>%
   as.data.frame()
  bc_rank_plot(
    stats = stats, raw_cells = barcodes, filt_cells = barcodes_d, 
    save = "r/barcode_rank.png"
  ) 
  print_HTML(
    seq_stats = seq_stats, cell_stats = cell_stats, 
    dir = "r/", sample_id = NULL                                        ## INPUT
  )
}
