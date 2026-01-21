library(GEOfastq)


gse_name <- 'GSE52778'
gse_text <- crawl_gse(gse_name)

gsm_names <- extract_gsms(gse_text)
gsm_name <- gsm_names[182]
srp_meta <- crawl_gsms(gsm_name)
