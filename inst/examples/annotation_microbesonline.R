## On those relatively rare occasions when the internet is being
## cranky, these tests will most certainly fail.

## Get the numeric ID from microbesonline for a species of interest.
ecoli_id <- try(get_microbesonline_taxid(species = "coli S88"))
ceoli_id

## Download a genbank/fasta/etc file for that species.
ecoli_files <- try(download_microbesonline_files(id = ecoli_id, type = "gbk"))
ecoli_files
file.exists(ecoli_files[["gbk"]])
removed <- file.remove(ecoli_files[["gbk"]])

## Download a species annotations to a dataframe
ecoli_annot <- try(load_microbesonline_annotations(id = ecoli_id))
dim(ecoli_annot)
colnames(ecoli_annot)

## Get gene ontology data from microbesonline
ecoli_go <- try(load_microbesonline_go(id = ecoli_id))
dim(ecoli_go)
