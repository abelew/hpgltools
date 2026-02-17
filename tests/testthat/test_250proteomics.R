start <- as.POSIXlt(Sys.time())
library(hpgldata)
context("250proteomics.R")

## Available functions:  add_conditional_nas(), extract_mayu_pps_fdr(),
## extract_scan_data(), extract_mzXML_scans(), extract_mzML_scans(),
## extract_msraw_data(), extract_peprophet_data(), extract_pyprophet_data(),
## impute_se(), mean_by_bioreplicate(), read_thermo_xlsx(), s2s_all_filters(),
## subset_pyprophet_data(), and gather_masses().

## Plotting functions:
## plot_intensity_mz(), plot_mzxml_boxplot(), plot_pyprophet_counts(),
## plot_pyprophet_xy(), plot_pyprophet_distribution(), plot_pyprophet_protein(),
## plot_pyprophet_points(), plot_peprophet_data(), plot_cleaved(),
## cleavage_histogram()


## I would like to test out my various proteomics functions.  Unfortunately,
## many of those functions require rather large input files.  So, I made a
## tarfile of the smallest of them in order to play with some of the
## functionality.

## As a reminder to myself, the general process followed is: (check out my
## 02_preprocessing*.Rmd for more details))
##  1.  Convert the thermo-fisher proprietary output files to mzXML/mzML using
##      MSreader on a windows virtual machine.
##  2.  Create a spectral library in the pqp format.  This may either be created
##      from an existing DDA experiment, or from some downloaded data.  Let us
##      assume the latter for the moment, then the command is
##      'TargetedFileConverter'
##  3.  Check over the mzXML/mzML files to ensure that there are no problems so
##      significant that they will cause later steps to fail. (That is what the
##      first few functions help with)
##  4.  Write appropriate window files for OpenSwathWorkflow, the function
##      'extract_msraw_data()' will do that for you, but it requires the rather
##      large mzXML files which I don't want to include here.  So unless I find
##      a data package containing some, I guess you will have to trust me that
##      those parsers actually work.
##  5.  Pass the transition library and raw data files to OpenSwathWorkFlow.
##  6.  Run pyprophet on the openswath outputs to normalize the error rates.
##  7.  Examine the pyprophet outputs to see if there are problems.
##      extract_pyprophet_data() and related functions help with that.
##  8.  Invoke SWATH2stats to filter the data and create matrices in the formats
##      expected by MSstats etc.
##  9.  Conversely, take those matrices and pass them to hpgltools.
##  10. Conversely, pass the mzXML data directly to encylopeDIA and extract the
##      output matrices.

## For most of the above choices, I have functions in the hpgltools to help.

meta <- system.file("share/mtb_prot/dia_samples.ods", package = "hpgldata")
untarred <- utils::untar(tarfile = system.file("share/mtb_prot/sb_prot.tar.xz",
                                               package = "hpgldata"))

## As the name implies, this function uses the dia_scored column in the metadata
## and reads the pyprophet output files indicated therein; it parses them into a
## list of data tables which are later used for plotting/examining.
## In that list is 'failed, colors, metadata, and sample_data'.  Failed is
## comprised of the files which failed to read properly, colors are for
## plotting, metadata is a copy of the metadata, and sample_data is the fun.
pyprophet_fun <- extract_pyprophet_data(metadata = meta,
                                        pyprophet_column = "dia_scored")

test_that("Did extract_pyprophet_data have failures?", {
  expect_equal(NULL, pyprophet_fun[["failed"]])
})
expected <- 35
test_that("Did extract_pyprophet_data provide the expected number of elements?", {
  expect_equal(expected, length(pyprophet_fun[["sample_data"]]))
})
## These are the columns provided by pyprophet.  They are columns provided by
## openswath with some extra scores appended by pyprophet.
expected <- c(
    "transition_group_id", "decoy", "run_id", "filename", "rt",  "assay_rt",
    "delta_rt", "irt", "assay_irt", "delta_irt", "id", "sequence",
    "fullpeptidename", "charge", "mz", "intensity", "aggr_prec_peak_area",
    "aggr_prec_peak_apex", "leftwidth", "rightwidth", "peak_group_rank",
    "d_score", "m_score", "aggr_peak_area", "aggr_peak_apex",
    "aggr_fragment_annotation", "proteinname", "mass", "seqlength")
test_that("Did extract_pyprotphet_data provide the expected columns?", {
  expect_equal(expected, colnames(pyprophet_fun[["sample_data"]][[1]]))
})
## The various plot_pyprophet functions may be used on pretty much any of the
## above columns, some of them are more useful than others for diagnosing
## problems.

plot_type <- "ggplot2::ggplot"
mass_plot <- plot_pyprophet_distribution(pyprophet_fun, column = "mass")
test_that("Does plot_pyprotphet_distribution return some plots?", {
  ## Yeah, so I couldn't decide which type of plot was best for representing this...
  expect_equal(class(mass_plot[["violin"]])[1], plot_type)
  expect_equal(class(mass_plot[["boxplot"]])[1], plot_type)
  expect_equal(class(mass_plot[["dotboxplot"]])[1], plot_type)
  expect_equal(class(mass_plot[["density"]])[1], plot_type)
})

## plot_pyprophet_counts has some parameters to ignore/include decoys
## as well as different types: 'protein_count' which counts how many proteins
## were found, 'intensity' which sums the intensities observed, 'count' which
## counts how many identifications were observed, otherwise it sums whatever
## columnname was provided, assuming it was numeric.
peptide_identifications <- plot_pyprophet_counts(pyprophet_fun, keep_decoys = FALSE,
                                                 type = "count")
test_that("Does plot_pyprophet_counts return some plots?", {
  expect_equal(class(peptide_identifications[["plot"]])[1], plot_type)
})

## widths with respect to counts are a surprisingly reliable way to find
## problematic samples.
pyprophet_lwidths <- plot_pyprophet_xy(pyprophet_fun, x_type = "count",
                                       y_type = "leftwidth")
test_that("Does plot_pyprophet_xy return a plot?", {
  expect_equal(class(pyprophet_lwidths)[1], plot_type)
})

## The above are global metrics, we can plot some metrics of individual
## proteins.  Any column from above may be used.
intensities_esxG <- plot_pyprophet_protein(pyprophet_fun, scale = "log",
                                           title = "esxG Intensities",
                                           column = "intensity", protein = "Rv0287")
test_that("Does plot_pyprophet_protein return a plot?", {
  expect_equal(class(intensities_esxG)[1], plot_type)
})

## The final step of pyprophet is to export the individual matrices into a
## single table, that is the input for SWATH2stats.
## I have one function which helps with that, s2s_all_filters() which just
## invokes all of the SWATH2stats filters with some hopefully sane defaults.
## The result of that set of filters is what I previously used to create an
## expressionset (now SE).

## In my previous notebook I did the following:
## 1.  extract_pyprophet_data() to get pyprophet_fun (done above)
## 2.  QC it with some plots (done above)
## 3.  Repeat with decoy data.
## 4.  Look at identified proteins/sample via plot_pyprophet_counts()
## (done above)
## 5.  Look for outliers with plot_pyprophet_xy and plot_pyprophet_points()
## 6.  Load the tric results into swath2stats from the tric-output csv
## file and the swath2stats function sample_annotation().
## This was done by just using read_tsv with 3 tric outputs and usin
## gsub to massage the protein IDs so that they follow the convention
## expected by swath2stats.
## I therefore copied those files to hpgldata as
## inst/share/cf_comet_HCD.tsv, wc_comet_HCD.csv, and
## previous_comet_HCD.csv
cf_tric_file <- system.file(file.path("share", "cf_comet_HCD.tsv"),
                            package = "hpgldata")
cf_tric_data <- readr::read_tsv(cf_tric_file)
cf_tric_data[["ProteinName"]] <- gsub(pattern="^(.*)_.*$", replacement="\\1",
                                      x=cf_tric_data[["ProteinName"]])
wc_tric_file <- system.file(file.path("share", "wc_comet_HCD.tsv"),
                            package = "hpgldata")
wc_tric_data <- readr::read_tsv(wc_tric_file)
wc_tric_data[["ProteinName"]] <- gsub(pattern="^(.*)_.*$", replacement="\\1",
                                      x=wc_tric_data[["ProteinName"]])
tric_data <- rbind(wc_tric_data, cf_tric_data)
prev_tric_file <- system.file(file.path("share", "previous_comet_HCD.tsv"),
                              package = "hpgldata")
prev_tric_data <- readr::read_tsv(prev_tric_file)
tric_data <- rbind(prev_tric_data, tric_data)

sample_annotation <- extract_metadata(meta)
## A weird thing is happening with column names, I will chase it down later.
colnames(sample_annotation) <- gsub(x = colnames(sample_annotation),
                               pattern = "[[:punct:]]", replace = "")
kept <- ! grepl(x = rownames(sample_annotation), pattern = "^s\\.\\.")
sample_annotation <- sample_annotation[kept, ]
rownames(sample_annotation) <- paste0("s", rownames(sample_annotation))

##devtools::load_all("~/scratch/git/SWATH2stats_myforked")
## This only works using my forked copy of SWATH2stats which allows
## one to specify the peptide column name.
##s2s_data <- sample_annotation(data = tric_data,
##                             sample_annotation = sample_annotation,
##                             fullpeptidename_column = "fullpeptidename")
s2s_file <- system.file(file.path("share", "s2s_data.rda"), package = "hpgldata")
## I therefore performed the above and copied the resulting csv file
## to hpgldata as s2s_data.csv
s2s_data <- load(s2s_file)

## With that in mind, this is one of my functions, but it just runs swath2stats filters.
s2s_filtered <- s2s_all_filters(s2s_data, target_fdr = 0.1, mscore = 0.1,
                                upper_fdr = 0.1, do_min = FALSE)
filtered_mtrx <- as.data.frame(write_matrix_proteins(s2s_filtered[["final"]]))
rownames(filtered_mtrx) <- filtered_mtrx[[1]]
filtered_mtrx[[1]] <- NULL
## Sadly, s2s leaves the weirdly formatted sample names in place, so I still need to do
## some manual massaging before I can use it.
reordered <- colnames(filtered_mtrx)
metadata <- sample_annotation[reordered, ]
colnames(filtered_mtrx) <- gsub(x = colnames(filtered_mtrx),
                                pattern = "^.*(2019.*$)", replacement="s\\1")
colnames(filtered_mtrx) <- gsub(x = colnames(filtered_mtrx),
                                pattern = "^.*(2018.*$)", replacement = "s\\1")

mtb_microbes <- load_microbesonline_annotations(id = 83332)
mtb_annotations <- as.data.frame(mtb_microbes)
rownames(mtb_annotations) <- make.names(mtb_annotations[["sysName"]], unique = TRUE)

existing_data <- rownames(sample_annotation) %in% colnames(filtered_mtrx)
sample_annot_subset <- sample_annotation[existing_data, ]
existing_data <- rownames(t(filtered_mtrx)) %in% rownames(sample_annotation)
t_filtered_mtrx <- t(filtered_mtrx)[existing_data, ]
filtered_mtrx <- t(t_filtered_mtrx)
reordered <- colnames(filtered_mtrx)
metadata <- sample_annot_subset[reordered, ]
metadata[["sampleid"]] <- rownames(metadata)
weird_ids <- grepl(x = rownames(filtered_mtrx), pattern = "_")
filtered_mtrx <- filtered_mtrx[!weird_ids, ]

rownames(metadata)
colnames(filtered_mtrx)

mtb_se <- create_se(metadata, count_dataframe = filtered_mtrx, gene_info = mtb_annotations)
plot_libsize(mtb_se)
plot_nonzero(mtb_se)

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 250proteomics.R in ", elapsed,  " seconds.")
