## 202601 eutils seems to be having some trouble.

if (FALSE) {
  test_accession <- "AE009948"
  gbk_file <- download_gbk(accessions = test_accession, write = TRUE)

  sagalacticae_genbank_annot <- load_genbank_annotations(accession = test_accession)
  dim(as.data.frame(sagalacticae_genbank_annot$cds))
}
