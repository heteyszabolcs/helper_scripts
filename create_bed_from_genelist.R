library("biomaRt")
library("tidyverse")
library("data.table")

# feed gene symbols and return bed file (tibble) with mm10 regions
# aug2020 archive returns mm10
mm10_gene_regions = function(gene_list) {
  ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host =
                      "https://aug2020.archive.ensembl.org")
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  location = getBM(
    attributes = c(
      'mgi_symbol',
      "chromosome_name",
      "start_position",
      "end_position"
    ),
    mart = ensembl,
    filters = 'mgi_symbol',
    values = gene_list
  )
  bed = location %>% mutate(seqname = paste0("chr", as.character(chromosome_name))) %>%
    dplyr::select(seqname, start_position, end_position)
  
  return(bed)
  
}
