## 202601 eutils seems to be having some trouble.

test_accession <- "AE009948"
gbk_file <- try(download_gbk(accessions = test_accession))
if (! "try-error" %in% class(gbk_file)) {
  removed <- file.remove(gbk_file[["written_file"]])
}

sagalacticae_genbank_annot <- load_genbank_annotations(accession = test_accession, type = "CDS")
if (! "try-error" %in% class(sagalacticae_genbank_annot)) {
  dim(sagalacticae_genbank_annot)
}
