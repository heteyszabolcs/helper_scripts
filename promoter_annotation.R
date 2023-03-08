suppressPackageStartupMessages({
  library("GenomicRanges")
  library("tidyverse")
  library("data.table")
  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("org.Hs.eg.db")
  library("org.Mm.eg.db")
  library("biomaRt")
  library("glue")
})


assign_promoter_region_mm10 = function(gene) {
  # promoter regions
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  genes = transcriptsBy(txdb, "gene")
  proms = promoters(genes, upstream = 1500, downstream = 500)
  proms = as_tibble(proms)
  proms = proms %>% separate(tx_name, into = c("tx_name", "id"), sep = "\\.") %>% dplyr::select(-id)
  
  # annotate promoters
  regions = GRanges(
    seqnames = proms$seqnames,
    ranges = IRanges(start = proms$start,
                     end = proms$end)
  )
  
  annot = annotatePeak(
    regions,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db"
  )
  
  annot = as.data.frame(annot)
  annot = annot %>% dplyr::select(seqnames, start, end, annotation, distanceToTSS, gene_symbol = SYMBOL) %>%
    dplyr::filter(gene_symbol %in% gene)
  
  return(annot)
  
}

assign_promoter_region_hg38 = function(gene) {
  # promoter regions
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  genes = transcriptsBy(txdb, "gene")
  proms = promoters(genes, upstream = 1500, downstream = 500)
  proms = as_tibble(proms)
  proms = proms %>% separate(tx_name, into = c("tx_name", "id"), sep = "\\.") %>% dplyr::select(-id)
  
  # annotate promoters
  regions = GRanges(
    seqnames = proms$seqnames,
    ranges = IRanges(start = proms$start,
                     end = proms$end)
  )
  
  annot = annotatePeak(
    regions,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )
  
  annot = as.data.frame(annot)
  annot = annot %>% dplyr::select(seqnames, start, end, annotation, distanceToTSS, gene_symbol = SYMBOL) %>%
    dplyr::filter(gene_symbol %in% gene)
  
  return(annot)
  
}
