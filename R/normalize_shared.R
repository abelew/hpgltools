## normalize_shared.r: Functions which bring together everything in the other
## normalize_* files.

#' Actually runs the batch method, this more than anything shows that
#' hpgl_norm is too complicated.
#'
#' @param count_table The counts in their current state.
#' @param method Batch/SV method to employ.
#' @param design Experimental design, requiring columns named
#'  'condition' and 'batch'.
#' @param batch1 Column containing the primary batch factor
#' @param current_state State of the data before messing with it.
#' @param current_design Redundant with design above...
#' @param expt_state Redundant with current_state above.
#' @param surrogate_method Method to estimate the number of SVs.
#' @param num_surrogates Explicit number of surrogates.
#' @param low_to_zero Push values less than zero to zero.
#' @param cpus Number of cores to use, should be removed.
#' @param batch2 Secondary batch factor column.
#' @param noscale Scale the data (combat)?
#' @param adjust_method Method to use to modify the counts after
#'  finding the surrogates.
#' @param batch_step Choose when to perform this in the set of normalization tasks.
do_batch <- function(count_table, method = "raw", design, batch1 = "batch",
                     current_state = NULL, current_design = NULL, expt_state = NULL,
                     surrogate_method = "be", num_surrogates = NULL, low_to_zero = FALSE,
                     cpus = 4, batch2 = NULL, noscale = TRUE, adjust_method = "ruv",
                     batch_step = 4) {
  retlist <- list(
    "result" = NULL,
    "batched_counts" = count_table,
    "batch_performed" = "raw",
    "adjust_performed" = "none")

  if (is.null(method)) {
    method <- "raw"
  }

  if (method == "raw") {
    mesg("Step ", batch_step, ": not doing batch correction.")
    return(retlist)
  } else {
    mesg("Step ", batch_step, ": doing batch correction with ",
         method, ".")
    result <- try(batch_counts(count_table, method = method,
                               design = design, adjust_method = adjust_method,
                               current_state = current_state, num_surrogates = num_surrogates))
    ##new_counts <- batch_counts(count_table, method = method,
    ##                           design = design, adjust_method = adjust_method,
    ##                           current_state = current_state)
    if ("try-error" %in% class(result)) {
      warning("The batch_counts call failed.  Returning non-batch reduced data.")
      return(retlist)
    } else {
      retlist[["result"]] <- result
      retlist[["batched_counts"]] <- result[["count_table"]]
      retlist[["batch_performed"]] <- method
      retlist[["adjust_performed"]] <- adjust_method
    }
  }
  return(retlist)
}

#' Every time I think I have a handle on R dispatch I get some BS like this.
#'
#' @param exp Input to normalize
#' @param ... other arguments passed along to BiocGenerics::normalize().
#' @export
normalize <- function(exp, ...) {
  message("This function is intended to normalize an experimental dataset.")
  message("It was passed an object of type ", class(exp),
          " and does not know what to do.")
  standardGeneric("normalize")
  return(NULL)
}
setGeneric("normalize")

#' Normalization of an expt
#'
#' The primary definition for normalize resides in BiocGenerics.
#'
#' @param object Data structure to normalize
#' @param ... Other options.
#' @include 01_hpgltools.R
#' @importFrom BiocGenerics normalize
#' @export
setMethod(
  "normalize", signature = signature(exp = "expt"),
  definition = function(exp, ...) {
    message("Running normalize.")
    normalize_expt(exp, ...)
  })

#' Normalization of se, taking a hint from BiocGenerics::normalize()
#'
#' @param object Data structure to normalize
#' @param ... Other options.
#' @export
setMethod(
  "normalize", signature = signature(exp = "SummarizedExperiment"),
  definition = function(exp, ...) {
    message("Running normalize_se.")
    normalize_se(exp, ...)
  })

#' Normalize the data of an expt object.  Save the original data, and note what
#' was done.
#'
#' It is the responsibility of normalize() to perform any arbitrary
#' normalizations desired as well as to ensure that the data integrity is
#' maintained.  In order to do this, it writes the actions performed in
#' expt$state and saves the intermediate steps of the normalization in
#' expt$intermediate_counts.  Furthermore, it should tell you every step of the
#' normalization process, from count filtering, to normalization, conversion,
#' transformation, and batch correction.
#'
#' @param expt Original expt.
#' @param transform Transformation desired, usually log2.
#' @param norm How to normalize the data? (raw, quant, sf, upperquartile, tmm, rle)
#' @param convert Conversion to perform? (raw, cpm, rpkm, cp_seq_m)
#' @param batch Batch effect removal tool to use? (limma sva fsva ruv etc)
#' @param filter Filter out low/undesired features? (cbcb, pofa, kofa, others?)
#' @param annotations Used for rpkm -- probably not needed as this is in fData now.
#' @param fasta Fasta file for cp_seq_m counting of oligos.
#' @param entry_type For getting genelengths by feature type (rpkm or cp_seq_m).
#' @param batch1 Experimental factor to extract first.
#' @param batch2 Second factor to remove (only with limma's removebatcheffect()).
#' @param batch_step From step 1-5, when should batch correction be applied?
#' @param low_to_zero When log transforming, change low numbers (< 0) to 0 to avoid NaN?
#' @param thresh Used by cbcb_lowfilter().
#' @param min_samples Also used by cbcb_lowfilter().
#' @param p Used by genefilter's pofa().
#' @param A Also used by genefilter's pofa().
#' @param k Used by genefilter's kofa().
#' @param cv_min Used by genefilter's cv().
#' @param cv_max Also used by genefilter's cv().
#' @param na_to_zero Sometimes rpkm gives some NA values for very low numbers.
#' @param adjust_method Given a set of sv estimates, change the counts with this method.
#' @param verbose Print what is happening while the normalization is performed?
#'  I am not sure why, but I think they should be 0.
#' @param num_surrogates Explicit number of way to seek the number
#'  of surrogates.
#' @param length_column annotation column containing gene lengths for rpkm.
#' @param ... more options
#' @return Expt object with normalized data and the original data saved as
#'  'original_expressionset'
#' @seealso [convert_counts()] [normalize_counts()] [batch_counts()]
#'  [filter_counts()] [transform_counts()]
#' @examples
#' \dontrun{
#'  normed <- normalize(exp, transform='log2', norm='rle', convert='cpm',
#'                           batch='raw', filter='pofa')
#'  normed_batch <- normalize(exp, transform='log2', norm='rle', convert='cpm',
#'                                 batch='sva', filter='pofa')
#' }
#' @export
normalize_expt <- function(expt, ## The expt class passed to the normalizer
                           ## choose the normalization strategy
                           transform = "raw", norm = "raw", convert = "raw",
                           batch = "raw", filter = FALSE,
                           ## annotations used for rpkm/cpseqm, original may be
                           ## used to ensure double-normalization isn't
                           ## performed.
                           annotations = NULL, fasta = NULL, entry_type = "gene",
                           batch1 = "batch", batch2 = NULL, batch_step = 4,
                           low_to_zero = TRUE, ## extra parameters for batch correction
                           thresh = 2, min_samples = 2, p = 0.01, A = 1, k = 1,
                           cv_min = 0.01, cv_max = 1000,  ## extra parameters for low-count filtering
                           na_to_zero = FALSE, adjust_method = "ruv", verbose = TRUE,
                           num_surrogates = "be", length_column = NULL,
                           ...) {
  arglist <- list(...)
  expt_state <- state(expt)
  new_expt <- expt
  type <- ""
  current_exprs <- expt[["expressionset"]]
  if (is.null(filter) || isFALSE(filter)) {
    filter <- "raw"
  } else if (isTRUE(filter)) {
    filter <- "cbcb"
  }
  if (is.null(convert) || isFALSE(convert)) {
    convert <- "raw"
  } else if (isTRUE(convert)) {
    convert <- "cbcbcpm"
  }
  if (is.null(norm) || isFALSE(norm)) {
    norm <- "raw"
  } else if (isTRUE(norm)) {
    norm <- "tmm"
  }
  if (is.null(transform) || isFALSE(transform)) {
    transform <- "raw"
  } else if (isTRUE(transform)) {
    transform <- "log2"
  }
  if (is.null(batch) || isFALSE(batch)) {
    batch <- "raw"
  } else if (isTRUE(batch)) {
    batch <- "sva"
  }

  mesg("This function will replace the expressionset with:")
  operations <- what_happened(transform = transform, batch = batch, convert = convert,
                              norm = norm, filter = filter)
  mesg(operations)
  mesg("It will save copies of each step along the way
 in expt$normalized with the corresponding libsizes. Keep libsizes in mind
 when invoking limma.  The appropriate libsize is non-log(cpm(normalized)).
 This is most likely kept at:
 'new_expt$normalized$intermediate_counts$normalization$libsizes'
 A copy of this may also be found at:
 new_expt$best_libsize
")

  if (filter == "raw") {
    mesg("Filter is false, this should likely be set to something, good
 choices include cbcb, kofa, pofa (anything but FALSE).  If you want this to
 stay FALSE, keep in mind that if other normalizations are performed, then the
 resulting libsizes are likely to be strange (potentially negative!)
")
  }
  if (transform == "raw") {
    mesg("Leaving the data in its current base format, keep in mind that
 some metrics are easier to see when the data is log2 transformed, but
 EdgeR/DESeq do not accept transformed data.
")
  }
  if (convert == "raw") {
    mesg("Leaving the data unconverted.  It is often advisable to cpm/rpkm
 the data to normalize for sampling differences, keep in mind though that rpkm
 has some annoying biases, and voom() by default does a cpm (though hpgl_voom()
 will try to detect this).
")
  }
  if (norm == "raw") {
    mesg("Leaving the data unnormalized.  This is necessary for DESeq, but
 EdgeR/limma might benefit from normalization.  Good choices include quantile,
 size-factor, tmm, etc.
")
  }
  if (batch == "raw") {
    mesg("Not correcting the count-data for batch effects.  If batch is
 included in EdgerR/limma's model, then this is probably wise; but in extreme
 batch effects this is a good parameter to play with.
")
  }
  if (convert == "cpm" && transform == "tmm") {
    warning("Cpm and tmm perform similar purposes. They should not be applied to the same data.")
  }
  if (norm == "quant" && isTRUE(grepl(x = batch, pattern = "sva"))) {
    warning("Quantile normalization and sva do not always play well together.")
  }

  data <- exprs(current_exprs)
  design <- pData(expt)
  if (is.null(annotations)) {
    annotations <- fData(current_exprs)
  }
  ## A bunch of these options should be moved into ...
  ## Having them as options to maintain is foolish
  if (isTRUE(verbose)) {
    normalized <- hpgl_norm(
      data, expt_state = expt_state, design = design, transform = transform, norm = norm,
      convert = convert, batch = batch, batch1 = batch1, batch2 = batch2, low_to_zero = low_to_zero,
      filter = filter, annotations = annotations, fasta = fasta, thresh = thresh,
      batch_step = batch_step, min_samples = min_samples, p = p, A = A, k = k, cv_min = cv_min,
      cv_max = cv_max, entry_type = entry_type, adjust_method = adjust_method,
      num_surrogates = num_surrogates, length_column = length_column,
      ...)
  } else {
    normalized <- sm(hpgl_norm(
      data, expt_state = expt_state, design = design, transform = transform, norm = norm,
      convert = convert, batch = batch, batch1 = batch1, batch2=batch2, low_to_zero = low_to_zero,
      filter = filter, annotations = annotations, fasta = fasta, thresh = thresh,
      batch_step = batch_step, min_samples = min_samples, p = p, A = A, k = k, cv_min = cv_min,
      cv_max = cv_max, entry_type = entry_type, adjust_method = adjust_method,
      num_surrogates = num_surrogates, length_column = length_column,
      ...))
  }

  final_libsize <- normalized[["libsize"]]
  final_data <- as.matrix(normalized[["count_table"]])

  ## A recent update to Biobase adds a test in the function
  ## assayDataElementReplace() which no longer allows one to just
  ## replace an expressionset with a smaller version (low-filtered).
  ## Instead, one must properly subset the object first, then replace.
  ## While this is annoying, I suppose it is smart.

  ## There is another peculiarity which has popped up in my proteomics data.
  ## Some rows are all NA!  How this happened I have not yet figured out,
  ## but these rows need to go.
  row_na_idx <- !is.na(rownames(final_data))
  final_data <- final_data[row_na_idx, ]
  unfiltered_genes <- rownames(exprs(current_exprs)) %in% rownames(final_data)
  current_exprs <- current_exprs[unfiltered_genes, ]
  ## This next line was added in response to an annoying occurance when
  ## combining technical replicates of a data set for which the samples started
  ## with numbers ('2018_0315' for example).  In that instance, my
  ## concatenate_runs() does not properly check to ensure that the column names
  ## do not change. I should therefore change that, but for the moment I will
  ## add a check here.
  colnames(final_data) <- sampleNames(current_exprs)
  exprs(current_exprs) <- final_data

  ## This state slot should match the information available in
  ## new_expt$normalized$actions
  ## I am hoping this will prove a more direct place to access it and provide a
  ## chance to double-check that things match
  new_state <- list(
    "filter" = normalized[["actions"]][["filter"]],
    "normalization" = normalized[["actions"]][["normalization"]],
    "conversion" = normalized[["actions"]][["conversion"]],
    "batch" = normalized[["actions"]][["batch"]],
    "transform" = normalized[["actions"]][["transform"]])

  ## Keep in mind that low_to_zero should be ignored if transform_state is not raw.
  if (new_state[["transform"]] == "raw" & isTRUE(low_to_zero)) {
    low_idx <- exprs(current_exprs) < 0
    ## Do not forget that na does not count when looking for numbers < x.
    na_idx <- is.na(exprs(current_exprs))
    low_idx[na_idx] <- FALSE
    low_num <- sum(low_idx)
    if (low_num > 0) {
      message("Setting ", low_num, " entries to zero.")
      exprs(current_exprs)[low_idx] <- 0
    }
  }

  if (isTRUE(na_to_zero)) {
    na_idx <- is.na(exprs(current_exprs))
    exprs(current_exprs)[na_idx] <- 0
  }

  ## The original data structure contains the following slots:
  ## colors, batches, convert, conditions, design, expressionset,
  ## filtered, initial_metadata, norm, original_expressionset,
  ## original_libsize, samplenames, stages, types, transform
  ## batches, convert, transform, filtered, normalization are being
  ## replaced with 'state' containing all of them.
  ## Similarly, the multiple libsizes maintained are going to be put into
  ## libsizes = list(original, norm, best, etc)
  ## Finally, I want to remove 'stages' and 'types' those are data in design.

  ## The structure of the 'normalized' data is fairly complex and includes the following:
  ## actions_performed -- a list of the actions done by hpgl_norm()
  ## intermediate_counts -- counts after each step in the hpgl_norm() process
  ## count_table -- a dataframe count table from hpgl_norm()
  ## libsize -- a final libsize from hpgl_norm()
  new_expt[["normalized"]] <- normalized
  new_expt[["state"]] <- new_state
  new_expt[["SVs"]] <- NULL
  if (!is.null(normalized[["sv_df"]])) {
    new_expt[["sv_df"]] <- normalized[["sv_df"]]
  }

  ## My concept of the 'best library size' comes from Kwame's work where the
  ## libsize was kept after performing quantile normalization, but before doing
  ## a log2(cpm()) The problem with this is that, if one does a normalize() then
  ## another normalize() then the assumptions used  may get violated.
  if (!is.null(normalized[["intermediate_counts"]][["normalization"]][["libsize"]])) {
    new_expt[["best_libsize"]] <-
      normalized[["intermediate_counts"]][["normalization"]][["libsize"]]
  } else if (!is.null(normalized[["intermediate_counts"]][["filter"]][["libsize"]])) {
    new_expt[["best_libsize"]] <-
      normalized[["intermediate_counts"]][["filter"]][["libsize"]]
  } else {
    new_expt[["best_libsize"]] <- NULL
  }
  ## limma should probably use this
  new_expt[["norm_result"]] <- normalized
  new_expt[["expressionset"]] <- current_exprs
  current_notes <- glue("{new_expt[['notes']]} Normalized with {type} at {date()}.
")
  new_expt[["notes"]] <- toString(current_notes)
  return(new_expt)
}
setGeneric("normalize")

#' If I call normalize on a SE, catch it and redispatch appropriately.
#'
#' @inheritParams normalize
#' @export
setMethod(
  "normalize", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, transform = "raw", norm = "raw", convert = "raw",
                        batch = "raw", filter = FALSE, annotations = NULL, fasta = NULL,
                        entry_type = "gene", batch1 = "batch",
                        batch2 = NULL, batch_step = 4, low_to_zero = TRUE,
                        thresh = 2, min_samples = 2, p = 0.01, A = 1, k = 1,
                        cv_min = 0.01, cv_max = 1000, na_to_zero = FALSE,
                        adjust_method = "ruv", verbose = FALSE, ...) {
    message("Running normalize_se")
    normalize_se(expt, transform = transform, norm = norm, convert = convert,
                 batch = batch, filter = filter, annotations = annotations,
                 fasta = fasta, entry_type = entry_type,
                 batch1 = batch1, batch2 = batch2, batch_step = batch_step,
                 low_to_zero = low_to_zero, thresh = thresh, min_samples = min_samples,
                 p = p, A = A, k = k, cv_min = cv_min, cv_max = cv_max,
                 na_to_zero = na_to_zero, adjust_method = adjust_method,
                 verbose = verbose, ...)
  })
setGeneric("normalize")

#' Normalize a SummarizedExperiment and think about how I want to reimplement some of this.
#'
#' @param se Summarized Experiment as input.
#' @param transform Transformation desired, usually log2.
#' @param norm How to normalize the data? (raw, quant, sf, upperquartile, tmm, rle)
#' @param convert Conversion to perform? (raw, cpm, rpkm, cp_seq_m)
#' @param batch Batch effect removal tool to use? (limma sva fsva ruv etc)
#' @param filter Filter out low/undesired features? (cbcb, pofa, kofa, others?)
#' @param annotations Used for rpkm -- probably not needed as this is in fData now.
#' @param fasta Fasta file for cp_seq_m counting of oligos.
#' @param entry_type For getting genelengths by feature type (rpkm or cp_seq_m).
#' @param batch1 Experimental factor to extract first.
#' @param batch2 Second factor to remove (only with limma's removebatcheffect()).
#' @param batch_step From step 1-5, when should batch correction be applied?
#' @param low_to_zero When log transforming, change low numbers (< 0)
#'  to 0 to avoid NaN?
#' @param thresh Used by cbcb_lowfilter().
#' @param min_samples Also used by cbcb_lowfilter().
#' @param p Used by genefilter's pofa().
#' @param A Also used by genefilter's pofa().
#' @param k Used by genefilter's kofa().
#' @param cv_min Used by genefilter's cv().
#' @param cv_max Also used by genefilter's cv().
#' @param na_to_zero Sometimes rpkm gives some NA values for very low numbers.
#' @param adjust_method Given a set of sv estimates, change the counts with this method.
#' @param verbose Print what is happening while the normalization is performed?
#'  I am not sure why, but I think they should be 0.
#' @param num_surrogates Explicit number of way to seek the number
#'  of surrogates.
#' @param surrogate_method Method to seek surrogates.
#' @param cpus Number of cpus to use.
#' @param noscale I think used only by combat -- scale the data?
#' @param length_column annotation column containing gene lengths for rpkm.
#' @param ... Arbitrary arguments.
#' @export
normalize_se <- function(se, ## The expt class passed to the normalizer
                         ## choose the normalization strategy
                         transform = "raw", norm = "raw", convert = "raw",
                         batch = "raw", filter = FALSE,
                         ## annotations used for rpkm/cpseqm, original may be
                         ## used to ensure double-normalization isn't
                         ## performed.
                         annotations = NULL, fasta = NULL, entry_type = "gene",
                         batch1 = "batch", batch2 = NULL, batch_step = 4,
                         low_to_zero = TRUE, ## extra parameters for batch correction
                         thresh = 2, min_samples = 2, p = 0.01, A = 1, k = 1,
                         cv_min = 0.01, cv_max = 1000,  ## extra parameters for low-count filtering
                         na_to_zero = FALSE, adjust_method = "ruv", verbose = TRUE,
                         num_surrogates = "be", surrogate_method = NULL, cpus = 4,
                         noscale = TRUE, length_column = NULL,
                         ...) {
  arglist <- list(...)

  ## I need to go through and rewrite how these get set:
  adjust_performed <- "none"
  batch_performed <- "none"
  convert_performed <- "none"
  filter_performed <- "none"
  transform_performed <- "none"

  ## A stopgap while I clean up these arguments
  if (is.null(surrogate_method) && "character" %in% class(num_surrogates)) {
    surrogate_method <- num_surrogates
  }

  count_table <- assay(se)
  meta <- S4Vectors::metadata(se)
  se_state <- state(se)
  current_state <- se_state
  original_libsize <- meta[["libsize"]]
  if (is.null(original_libsize)) {
    original_libsize <- colSums(count_table)
  }
  if (is.null(annotations)) {
    annotations <- rowData(se)
  }

  design <- colData(se)
  original_counts <- count_table

  type <- ""
  if (is.null(filter) || isFALSE(filter)) {
    filter <- "raw"
  } else if (isTRUE(filter)) {
    filter <- "cbcb"
  }
  if (is.null(convert) || isFALSE(convert)) {
    convert <- "raw"
  } else if (isTRUE(convert)) {
    convert <- "cbcbcpm"
  }
  if (is.null(norm) || isFALSE(norm)) {
    norm <- "raw"
  } else if (isTRUE(norm)) {
    norm <- "tmm"
  }
  if (is.null(transform) || isFALSE(transform)) {
    transform <- "raw"
  } else if (isTRUE(transform)) {
    transform <- "log2"
  }
  if (is.null(batch) || isFALSE(batch)) {
    batch <- "raw"
  } else if (isTRUE(batch)) {
    batch <- "sva"
  }

  mesg("This will normalize the data via:")
  operations <- what_happened(transform = transform, batch = batch, convert = convert,
                              norm = norm, filter = filter)
  mesg(operations)

  if (filter == "raw") {
    mesg("Filter is false, this should likely be set to something, good
 choices include cbcb, kofa, pofa (anything but FALSE).  If you want this to
 stay FALSE, keep in mind that if other normalizations are performed, then the
 resulting libsizes are likely to be strange (potentially negative!)
")
  }
  if (transform == "raw") {
    mesg("Leaving the data in its current base format, keep in mind that
 some metrics are easier to see when the data is log2 transformed, but
 EdgeR/DESeq do not accept transformed data.
")
  }
  if (convert == "raw") {
    mesg("Leaving the data unconverted.  It is often advisable to cpm/rpkm
 the data to normalize for sampling differences, keep in mind though that rpkm
 has some annoying biases, and voom() by default does a cpm (though hpgl_voom()
 will try to detect this).
")
  }
  if (norm == "raw") {
    mesg("Leaving the data unnormalized.  This is necessary for DESeq, but
 EdgeR/limma might benefit from normalization.  Good choices include quantile,
 size-factor, tmm, etc.
")
  }
  if (batch == "raw") {
    mesg("Not correcting the count-data for batch effects.  If batch is
 included in EdgerR/limma's model, then this is probably wise; but in extreme
 batch effects this is a good parameter to play with.
")
  }
  if (convert == "cpm" && transform == "tmm") {
    warning("Cpm and tmm perform similar purposes. They should not be applied to the same data.")
  }
  if (norm == "quant" && isTRUE(grepl(x = batch, pattern = "sva"))) {
    warning("Quantile normalization and sva do not always play well together.")
  }

  count_lst <- list(
    "count_table" = count_table,
    "libsize" = original_libsize)
  current_libsize <- original_libsize
  current_mtrx <- as.matrix(count_table)

  batched_counts <- NULL
  normalized_counts <- NULL
  filtered_counts <- NULL
  sv_df <- NULL
  if (batch_step == 1) {
    batch_data <- do_batch(
      count_lst, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = se_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    current_libsize <- batch_data[["libsize"]]
    count_table <- batch_data[["batched_counts"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
  }

  ## Step 1: count filtering
  mesg("Step 1: performing count filter with option: ", filter)
  filter_performed <- filter
  ## All the other intermediates have a libsize slot, perhaps this should too
  if (filter == "raw") {
    mesg("Step 1: not filtering the data.")
  } else {
    mesg("Step 1: filtering the data with ", filter, ".")
    filtered_counts <- filter_counts(count_table, method = filter, p = p,
                                     A = A, k = k, cv_min = cv_min, cv_max = cv_max,
                                     thresh = thresh, min_samples = min_samples,
                                     ...)
    ## filtered_counts <- filter_counts(count_table, method = filter)
    current_libsize <- filtered_counts[["libsize"]]
    count_table <- filtered_counts[["count_table"]]
    current_state[["filter"]] <- filter
  }

  if (batch_step == 2) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = se_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    current_libsize <- batch_data[["libsize"]]
    count_table <- batch_data[["batched_counts"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
  }

  normalized <- count_table
  norm_performed <- norm
  if (norm == "raw") {
    mesg("Step 2: not normalizing the data.")
  } else {
    mesg("Step 2: Normalizing the data with ", norm, ".")
    normalized_counts <- normalize_counts(data = count_table, method = norm,
                                          ...)
    current_libsize <- normalized_counts[["libsize"]]
    count_table <- normalized_counts[["count_table"]]
    current_state[["normalization"]] <- norm
  }

  ## Step 3: Convert the data to (likely) cpm
  ## The following stanza handles the three possible output types
  ## cpm and rpkm are both from edgeR
  ## They have nice ways of handling the log2 which I should consider
  if (batch_step == 3) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = se_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    current_libsize <- batch_data[["libsize"]]
    count_table <- batch_data[["count_table"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
  }
  converted_counts <- count_table
  convert_performed <- convert
  if (convert == "raw") {
    mesg("Step 3: not converting the data.")
  } else {
    mesg("Step 3: converting the data with ", convert, ".")
    converted_counts <- convert_counts(count_table, method = convert, annotations = annotations,
                                       length_column = length_column,
                                       ...)
    ## converted_counts <- convert_counts(count_table, method = convert, annotations = annotations)
    current_libsize <- converted_counts[["libsize"]]
    count_table <- converted_counts[["count_table"]]
    convert_performed <- convert
    current_state[["conversion"]] <- convert
  }

  ## Step 4: Transformation
  ## Finally, this considers whether to log2 the data or no
  if (batch_step == 4) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = se_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    current_libsize <- batch_data[["libsize"]]
    count_table <- batch_data[["batched_counts"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
    ## count_table <- do_batch(count_table, method = batch,
    ## design = design, current_state = current_state,
    ## num_surrogates = num_surrogates)
  }

  transformed_counts <- count_table
  transform_performed <- transform
  if (transform == "raw") {
    mesg("Step 4: not transforming the data.")
  } else {
    mesg("Step 4: transforming the data with ", transform, ".")
    transformed_counts <- transform_counts(count_table, method = transform,
                                           ...)
    ## transformed_counts <- transform_counts(count_table, transform = transform)
    current_libsize <- transformed_counts[["libsize"]]
    count_table <- transformed_counts[["count_table"]]
    if (transform == "round") {
      transform_performed <- "raw"
      transform <- "raw"
      current_state[["rounded"]] <- TRUE
    } else {
      transform_performed <- transform
    }
    current_state[["transform"]] <- transform
  }

  if (batch_step == 5) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = se_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    current_libsize <- batch_data[["libsize"]]
    count_table <- batch_data[["count_table"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
  }

  impute <- "raw"
  imputed_counts <- count_table
  if (!is.null(arglist[["impute"]])) {
    impute <- arglist[["impute"]]
  }
  impute_performed <- impute
  if (impute != "raw") {
    imputed_counts <- impute_expt(count_table)
    current_libsize <- imputed_counts[["libsize"]]
    count_table <- imputed_counts[["count_table"]]
    current_state[["impute"]] <- impute
  }

  if (isTRUE(low_to_zero)) {
    low_idx <- count_table < 0
    ## Do not forget that na does not count when looking for numbers < x.
    na_idx <- is.na(count_table)
    low_idx[na_idx] <- FALSE
    low_num <- sum(low_idx)
    if (low_num > 0) {
      message("Setting ", low_num, " entries to zero.")
      count_table[low_idx] <- 0
    }
  }

  if (isTRUE(na_to_zero)) {
    na_idx <- is.na(count_table)
    count_table[na_idx] <- 0
  }

  ## This list provides the list of operations performed on the data in order
  ## they were done.
  actions <- list(
    "num_surrogates" = num_surrogates,
    "filter" = filter_performed,
    "normalization" = norm_performed,
    "conversion" = convert_performed,
    "batch" = batch_performed,
    "adjust" = adjust_performed,
    "transform" = transform_performed)

  meta[["libsize"]] <- current_libsize
  ## meta[["state"]] <- actions_performed
  meta[["SVs"]] <- NULL
  if (!is.null(sv_df)) {
    meta[["SVs"]] <- sv_df
  }

  kept_genes <- rownames(annotations) %in% rownames(count_table)
  annotations <- annotations[kept_genes, ]

  se <- SummarizedExperiment(assays = count_table,
                             rowData = annotations,
                             colData = design)
  metadata(se) <- meta
  state(se) <- actions
  return(se)
}
setGeneric("normalize_se")

#' Normalize a dataframe/expt, express it, and/or transform it
#'
#' There are many possible options to this function.  Refer to normalize()
#' for a more complete list.
#'
#' FIXME: This function is defunct and should be deleted in favor of normalize()
#'
#' @param data Some data as a df/expt/whatever.
#' @param expt_state State of the input data.
#' @param design experimental design.
#' @param transform Desired transformation.
#' @param norm Desired normalization.
#' @param convert Desired conversion.
#' @param batch Desired batch estimator/adjuster.
#' @param batch1 Factor containing the primary batch factor.
#' @param batch2 Factor containing secondary batch factor.
#' @param low_to_zero Push values < 0 to zero.
#' @param filter Filter the data before SV searching.
#' @param annotations Annotation source.
#' @param fasta Input fasta file (for rpkm I assume)
#' @param thresh threshold number of reads/gene.
#' @param batch_step When to perform batch correction.
#' @param min_samples Requisite number of samples for collapsing.
#' @param p p-value cutoff.
#' @param k K!
#' @param A A
#' @param cv_min Minimum coefficient of variance.
#' @param cv_max Maximum allowed coefficient of variance.
#' @param entry_type Specific annotation type for matching rpkm lengths.
#' @param adjust_method Adjustment method when seeking SVs
#' @param num_surrogates Explicit number or character describing how  many SVs
#' @param length_column annotation column containing gene lengths.
#' @param cpus deprecated: Number of cpus to use.
#' @param noscale Used by combat.
#' @param ... I should put all those other options here
#' @return edgeR's DGEList expression of a count table.  This seems to
#'  me to be the easiest to deal with.
#' @seealso [edgeR] [DESeq2] [edgeR::cpm()] [filter_counts()] [batch_counts()]
#'  [convert_counts()] [transform_counts()]
#' @export
#' @examples
#' \dontrun{
#'  df_raw = hpgl_norm(expt = expt)  ## Only performs low-count filtering
#'  df_raw = hpgl_norm(df = a_df, design = a_design) ## Same, but using a df
#'  df_ql2rpkm = hpgl_norm(expt = expt, norm='quant', transform='log2',
#'                         convert='rpkm')  ## Quantile, log2, rpkm
#'  count_table = df_ql2rpkm$counts
#' }
hpgl_norm <- function(data, expt_state = NULL, design = NULL, transform = "raw",
                      norm = "raw", convert = "raw", batch = "raw", batch1 = "batch",
                      batch2 = NULL, low_to_zero = TRUE, filter = "raw", annotations = NULL,
                      fasta = NULL, thresh = 2, batch_step = 4, min_samples = 2, p = 0.05,
                      k = 0.01, A = 1000, cv_min = 0.001, cv_max = 1, entry_type = "gene",
                      adjust_method = "ruv", num_surrogates = "be", length_column = "length",
                      cpus = 4, noscale = FALSE, ...) {
  arglist <- list(...)
  filter_performed <- "raw"
  norm_performed <- "raw"
  convert_performed <- "raw"
  transform_performed <- "raw"
  batch_performed <- "raw"
  adjust_performed <- "none"
  expt_state <- list(
    "low_filter" = filter_performed,
    "normalization" = norm_performed,
    "conversion" = convert_performed,
    "batch" = batch_performed,
    "adjust" = adjust_performed,
    "transform" = transform_performed)
  data_class <- class(data)[1]
  original_counts <- NULL
  original_libsize <- NULL
  annot <- NULL
  counts <- NULL
  surrogate_method <- num_surrogates
  ## I never quite realized just how nice data.tables are.  To what extent can I refactor
  ## all of my data frame usage to them?
  if (data_class == "expt") {
    original_counts <- data[["original_counts"]]
    original_libsizes <- data[["original_libsize"]]
    design <- pData(data)
    annot <- fData(data)
    counts <- exprs(data)
    expt_state <- state(expt)
  } else if (data_class == "ExpressionSet") {
    counts <- exprs(data)
    design <- pData(data)
    annot <- fData(data)
  } else if (data_class == "list") {
    counts <- data[["count_table"]]
    if (is.null(data)) {
      stop("The list provided contains no count_table.")
    }
  } else if (data_class == "matrix" ||
               data_class == "data.frame" ||
               data_class == "data.table") {
    counts <- as.data.frame(data)  ## some functions prefer matrix, so I am
    ## keeping this explicit for the moment. In the case of data.tables, even if
    ## you set the rownames, the first column might still be rowname characters
    ## I don't yet fully understand this, so I will add an explicit test here.
    if (data_class == "data.table" & class(counts[[1]]) == "character") {
      rownames(counts) <- make.names(counts[[1]], unique = TRUE)
      counts <- counts[-1]
    }
  } else {
    stop("This only understands types: expt, ExpressionSet, data.frame, and matrix.")
  }

  count_table <- as.matrix(counts)
  if (is.null(original_counts)) {
    original_counts <- counts
  }
  if (is.null(original_libsize)) {
    original_libsize <- colSums(count_table)
  }
  if (is.null(annotations) && !is.null(annot)) {
    annotations <- annot
  }

  ## Make changes to this as we go.
  current_state <- expt_state
  batched_counts <- NULL
  if (!is.numeric(batch_step)) {
    batch_step <- 4
  } else if (batch_step > 5 | batch_step < 0) {
    batch_step <- 4
  }

  sv_df <- NULL
  if (batch_step == 1) {
    batch_data <- do_batch(count_table, method = batch, design = design,
                           batch1 = batch1, current_state = current_state,
                           expt_state = expt_state, surrogate_method = surrogate_method,
                           num_surrogates = num_surrogates, low_to_zero = low_to_zero,
                           cpus = cpus, batch2 = batch2, noscale = noscale,
                           adjust_method = adjust_method)
    count_table <- batch_data[["batched_counts"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
  }

  ## Step 1: count filtering
  filtered_counts <- NULL
  if (filter == FALSE | filter == "raw") {
    mesg("Step 1: not doing count filtering.")
  } else {
    if (isTRUE(filter)) {
      filter <- "cbcb"
    }
    mesg("Step 1: performing count filter with option: ", filter, " and ", thresh, ".")
    ## All the other intermediates have a libsize slot, perhaps this should too
    filtered_counts <- filter_counts(count_table, method = filter,
                                     ...)
    ## filtered_counts <- filter_counts(count_table, method = filter)
    count_table <- filtered_counts[["count_table"]]
    filter_performed <- filter
    current_state[["filter"]] <- filter
  }

  if (batch_step == 2) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = expt_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    count_table <- batch_data[["count_table"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
  }
  ## Step 2: Normalization
  ## This section handles the various normalization strategies
  normalized_counts <- NULL
  if (norm == "raw") {
    mesg("Step 2: not normalizing the data.")
  } else {
    mesg("Step 2: normalizing the data with ", norm, ".")
    if (is.null(design)) {
      message("The experimental design is null.  Some normalizations will fail.")
      message("If you get an error about 'no dimensions', that is likely why.")
    }
    normalized_counts <- normalize_counts(data = count_table, design = design,
                                          method = norm)
    count_table <- normalized_counts[["count_table"]]
    norm_performed <- norm
    current_state[["normalization"]] <- norm
  }

  ## Step 3: Convert the data to (likely) cpm
  ## The following stanza handles the three possible output types
  ## cpm and rpkm are both from edgeR
  ## They have nice ways of handling the log2 which I should consider
  if (batch_step == 3) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = expt_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    count_table <- batch_data[["count_table"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
  }

  converted_counts <- NULL
  if (convert == "raw") {
    mesg("Step 3: not converting the data.")
  } else {
    mesg("Step 3: converting the data with ", convert, ".")
    converted_counts <- convert_counts(count_table, method = convert, annotations = annotations,
                                       length_column = length_column,
                                       ...)
    ## converted_counts <- convert_counts(count_table, method = convert, annotations = annotations)
    count_table <- converted_counts[["count_table"]]
    convert_performed <- convert
    current_state[["conversion"]] <- convert
  }

  ## Step 4: Transformation
  ## Finally, this considers whether to log2 the data or no
  if (batch_step == 4) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = expt_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    count_table <- batch_data[["batched_counts"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
    ## count_table <- do_batch(count_table, method = batch,
    ##                         design = design, current_state = current_state)
  }
  transformed_counts <- NULL
  if (transform == "raw") {
    mesg("Step 4: not transforming the data.")
  } else {
    mesg("Step 4: transforming the data with ", transform, ".")
    transformed_counts <- transform_counts(count_table, method = transform,
                                           ...)
    ## transformed_counts <- transform_counts(count_table, transform = transform)
    count_table <- transformed_counts[["count_table"]]
    if (transform == "round") {
      transform_performed <- "raw"
      transform <- "raw"
      current_state[["rounded"]] <- TRUE
    } else {
      transform_performed <- transform
    }
    current_state[["transform"]] <- transform
  }

  if (batch_step == 5) {
    batch_data <- do_batch(
      count_table, method = batch, design = design, batch1 = batch1,
      current_state = current_state, expt_state = expt_state,
      surrogate_method = surrogate_method, num_surrogates = num_surrogates,
      low_to_zero = low_to_zero, cpus = cpus, batch2 = batch2, noscale = noscale,
      adjust_method = adjust_method)
    count_table <- batch_data[["count_table"]]
    batched_counts <- count_table
    batch_performed <- batch_data[["batch_performed"]]
    sv_df <- batch_data[["result"]][["result"]][["model_adjust"]]
    ## count_table <- do_batch(count_table, arglist)
  }

  ## This list provides the list of operations performed on the data in order
  ## they were done.
  actions <- list(
    "filter" = filter_performed,
    "normalization" = norm_performed,
    "conversion" = convert_performed,
    "batch" = batch_performed,
    "adjust" = adjust_performed,
    "transform" = transform_performed,
    "num_surrogates" = num_surrogates)
  ## This list contains the intermediate count tables generated at each step
  ## This may be useful if there is a problem in this process.
  ## Each of them also contains the libsize at that point in the process.
  intermediate_counts <- list(
    "original" = original_counts, ## The original count table, should never
    ## change from iteration to iteration
    "input" = as.matrix(data),  ## The input provided to this function, this may
    ## diverge from original
    "filter" = filtered_counts,  ## After filtering
    "normalization" = normalized_counts,  ## and normalization
    "conversion" = converted_counts,  ## and conversion
    "batch" = batched_counts,  ## and batch correction
    "transform" = transformed_counts)  ## and finally, transformation.

  retlist <- list(
    "actions" = actions,
    "intermediate_counts" = intermediate_counts,
    "count_table" = count_table,  ## The final count table
    "final_state" = current_state,
    "sv_df" = sv_df,
    "libsize" = colSums(count_table, na.rm = TRUE)  ## The final libsizes
  )
  return(retlist)
}

#' Simplified and ideally improved normalization function
#'
#' This function is ideally should provide a simpler and more capable
#' version of normalize. I also want to move everything to using
#' summarizedExperiments and this simpler method provides an
#' opportunity.
#'
#' @param expt Input data
#' @param todo List of tasks to perform.
normalize_todo <- function(expt, todo = list()) {
  ## This expects a list like:
  ## list("norm" = "quant", "filter" = c(pofa, "A" = 1))
  possible_methods <- list(
    "transform" = "transform_counts",
    "norm" = "normalize_counts",
    "convert" = "convert_counts",
    "filter" = "filter_counts",
    "batch" = "batch_counts",
    "impute" = "impute_counts")
  annot <- fData(expt)
  counts <- exprs(expt)
  design <- pData(expt)
  count_table <- list("count_table" = counts,
                      "libsize" = expt[["libsize"]])
  for (i in seq_along(todo)) {
    type <- names(todo)[i]
    operation <- todo[[i]]
    method <- operation[[1]]
    args <- c()
    if (length(operation) > 1) {
      args <- operation[2:length(operation)]
    }
    if (! type %in% names(possible_methods)) {
      stop("This type of todo is not known: ", type, ".")
    }
    call <- possible_methods[[type]]
    arglist <- list()
    arglist[["count_table"]] <- count_table
    arglist[["method"]] <- method
    for (a in args) {
      name <- names(args)[a]
      arglist[[a]] <- args[[a]]
    }
    mesg("Invoking: ", call)
    count_table <- base::do.call(call, arglist)
  }
  return(count_table)
}

## EOF
