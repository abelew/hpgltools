## On those relatively rare occasions when the internet is being
## cranky, these tests will most certainly fail.

## Get the numeric ID from microbesonline for a species of interest.
## When microbesonline is under maintenance/etc this fails.
ecoli_id <- try(get_microbesonline_taxid(species = "S88"))
ecoli_id

if (! "try-error" %in% class(ecoli_id)) {
  ## Download a genbank/fasta/etc file for that species.
  ecoli_files <- try(download_microbesonline_files(id = ecoli_id, type = "gbk"))
  ecoli_files
  if (! "try-error" %in% class(ecoli_files)) {
    removed <- file.remove(ecoli_files[["gbk"]])
  }

  ## Download a species annotations to a dataframe
  ecoli_annot <- try(load_microbesonline_annotations(id = ecoli_id))
  if (! "try-error" %in% class(ecoli_annot)) {
    dim(ecoli_annot)
    colnames(ecoli_annot)
  }

  ## Get gene ontology data from microbesonline
  ecoli_go <- try(load_microbesonline_go(id = ecoli_id))
  if (! "try-error" %in% class(ecoli_go)) {
    dim(ecoli_go)
  }
}
