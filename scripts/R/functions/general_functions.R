# function to read GTF files
read_gtf <- function(x){
  # read GTF
  gtf <- rtracklayer::import(x) %>% 
    as_tibble()
  
  # remove genes with PAR_Y (in pseudoautosomal region)
  gtf <- gtf %>% filter(!gene_id %>% str_detect("PAR_Y$"))
  
  # remove ensembl version number from gtf annotation
  gtf$gene_id <- gtf$gene_id %>% str_remove("\\..*$")
  
  # extract protein coding genes only from gtf
  gtf <- gtf %>% filter(type == "gene", gene_type == "protein_coding")
  
  # return 
  return(gtf)
}
