start <- as.POSIXlt(Sys.time())
context("470se_spyogenes.R: Does a small bacterial RNAseq experiment load?")

mgas_data <- new.env()
cdm_data <- system.file("share/cdm_se.rda", package = "hpgldata")
load(cdm_data, envir = mgas_data)
rm(cdm_data)

mgas_se <- sm(create_se(count_dataframe = mgas_data[["cdm_counts"]],
                        metadata = mgas_data[["cdm_metadata"]],
                        gene_info = mgas_data[["gene_info"]]))

expected <- c("dnaA", "dnaN", "M5005_Spy_0003", "ychF", "pth", "trcF")
actual <- head(fData(mgas_se)[["Name"]])
test_that("Did the gene information load?", {
  expect_equal(expected, actual)
})

mgas_norm <- normalize(mgas_se, transform = "log2", convert = "cbcbcpm", filter = TRUE)
mgas_state <- state(mgas_norm)
test_that("Is the filter state maintained?", {
  expect_equal("cbcb", mgas_state[["filter"]])
})
test_that("Is the normalization state maintained?", {
  expect_equal("raw", mgas_state[["normalization"]])
})
test_that("Is the conversion state maintained?", {
  expect_equal("cbcbcpm", mgas_state[["conversion"]])
})
test_that("Is the transformation state maintained?", {
  expect_equal("log2", mgas_state[["transform"]])
})

mgas_norm <- normalize(mgas_norm, batch = "combat_scale")
mgas_state <- state(mgas_norm)
test_that("Is the batch state maintained?", {
  expect_equal("combat_scale", mgas_state[["batch"]])
})

mgas_pairwise <- all_pairwise(mgas_se)
expected <- 0.50
actual <- min(mgas_pairwise[["comparison"]][["comp"]])
test_that("Do we get reasonably high similarities among the various DE tools?", {
  expect_gt(actual, expected)
})

mgas_combined <- combine_de_tables(mgas_pairwise, excel = FALSE)
mgas_sig <- extract_significant_genes(mgas_combined, excel = FALSE)
expected <- 150
actual <- nrow(mgas_sig[["deseq"]][["ups"]][["wt_ll_cf_vs_mga1_ll_cf"]])
test_that("Do we find some significant genes in the mga/wt fructose analysis?", {
  expect_gt(actual, expected)
})

mgas_data <- load_genbank_annotations(accession = "AE009949")
expected <- 1895017
actual <- Biostrings::width(mgas_data[["sequence"]])
actual_width <- actual
test_that("Can I extract the chromosome sequence from a genbank file? (widths)", {
  expect_equal(expected, actual)
})

cds_features <- mgas_data[["feature_list"]][["CDS"]]

expected <- c(1845, 10)
actual <- dim(cds_features)
test_that("Can I extract the chromosome sequence from a genbank file? (exons)", {
  expect_equal(expected, actual)
})

expected <- c("spyM18_0001", "spyM18_0002", "spyM18_0004",
              "spyM18_0005", "spyM18_0007", "spyM18_0008")
actual <- head(cds_features[["locus_tag"]])
test_that("Can I extract the chromosome sequence from a genbank file? (gene names)", {
  expect_equal(expected, actual)
})

taxon <- "293653"
mgas_df <- load_microbesonline_annotations(id = taxon)
mgas_df[["sysName"]] <- gsub(pattern = "Spy_", replacement = "Spy", x = mgas_df[["sysName"]])
expected <- c("dnaA","dnaN","M5005_Spy_0003","M5005_Spy_0004","pth","trcF")
actual <- as.character(head(mgas_df[["name"]]))
test_that("Did the mgas annotations download?", {
  expect_equal(expected, actual)
})

mgas_go <- load_microbesonline_go(taxon)
colnames(mgas_go) <- c("ID", "GO")
mgas_go <- unique(mgas_go)
expected <- c(4161, 2)
actual <- dim(mgas_go)
test_that("Do we get expected gene ontology information?", {
  expect_equal(expected, actual)
})

circos_annot_df <- as.data.frame(mgas_df)
circos_annot_df <- circos_annot_df[, c("start", "stop", "strand", "COGFun")]
circos_annot_df[["chromosome"]] <- "chr1"
rownames(circos_annot_df) <- make.names(gsub(x = mgas_df[["sysName"]],
                                             pattern = "Spy_", replacement = "Spy"),
                                        unique = TRUE)

## There is no way circos will work on travis, lets be realistic.
##if (identical(Sys.getenv("HAS_CIRCOS"), "true")) {
## Plot the coefficients of latelog glucose
glucose_table <- mgas_pairwise[["limma"]][["identity_tables"]][["mga1_ll_cg"]]
wtvmga_glucose <- mgas_pairwise[["limma"]][["all_tables"]][["wt_ll_cg_vs_mga1_ll_cg"]]
relevant_widths <- merge(glucose_table, mgas_df, by.x = "row.names",
                         by.y = "sysName", all.x = TRUE)
## Since genbankr died, get the gene lengths from microbesonline
relevant_widths <- suppressWarnings(as.numeric(relevant_widths[["width"]]))
na_widths <- is.na(relevant_widths)
relevant_widths[na_widths] <- 0

na_cog <- is.na(circos_annot_df[["COGFun"]])
circos_annot_df[na_cog, "COGFun"] <- "Z"

circos_test <- circos_prefix(circos_annot_df, name = "mgas",
                             chr_column = "chromosome", stop_column = "stop")
lengths <- 1835600
names(lengths) <- "chr1"
circos_kary <- circos_karyotype(circos_test, lengths = lengths)
circos_plus <- circos_plus_minus(circos_test)
circos_hist_ll_cg <- circos_hist(circos_test, glucose_table,
                                 colname = "logFC", outer = circos_plus)
circos_tile_wtmga <- circos_tile(circos_test, wtvmga_glucose,
                                 colname = "logFC", outer = circos_hist_ll_cg)
circos_suffix(circos_test)
circos_made <- circos_make(circos_test, target = "mgas")
expected <- "circos/mgas.svg"
test_that("Did circos run?", {
  expect_true(file.exists(expected))
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = (as.numeric(end - start)))
message("\nFinished 70se_spyogenes.R in ", elapsed,  " seconds.")
