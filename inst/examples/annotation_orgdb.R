hs_annot <- load_orgdb_annotations("org.Hs.eg.db")
dim(hs_annot[["genes"]])

## I changed this fairly significantly to simplify the select statements.
hs_go <- load_orgdb_go("org.Mm.eg.db")
summary(hs_go)

## extract_eupath_orthologs should be moved to EuPathDB

## This is similar to clusterProfiler's bitr I think.
dm_symbols <- map_orgdb_ids("org.Dm.eg.db", keytype = "entrezid", mapto = "symbol")
head(dm_symbols)

## Guess the most appropriate keytype given some IDs
## Essentially: what column has this?
wanted_ids <- c("G9a", "cin", "ewg")
keytype_guess <- guess_orgdb_keytype(ids = wanted_ids, "org.Dm.eg.db")
keytype_guess

## Extract The actual data of interest from AnnotationHub
if (interactive()) {
  wanted <- orgdb_from_ah(type = "OrgDb", title = "Gorilla")
  class(wanted)
}
