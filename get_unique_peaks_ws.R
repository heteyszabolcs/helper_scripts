# packages
suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("ggpubr")
  library("cowplot")
  library("Seurat")
  library("wigglescout")
  library("GenomicRanges")
  library("ChIPseeker")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
})

# get unique peaks by wigglescout
get_unique_ws = function(bw,
                         bw_backgr,
                         subset,
                         thr = 4) {
  label = "fold_change"
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  
  read_cov = bw_loci(
    bwfiles = bw,
    bg_bwfiles = bw_backgr,
    labels = label,
    loci = subset
  )
  
  annot = annotatePeak(
    read_cov,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db"
  )
  annot = as.data.frame(annot)
  annot = annot %>% dplyr::select(seqnames,
                                  start,
                                  end,
                                  starts_with("fold_change"),
                                  distanceToTSS,
                                  gene_symbol = SYMBOL) %>% dplyr::filter(abs(fold_change) > thr) %>%
    dplyr::filter(!fold_change == Inf)
  
  return(annot)
  
}