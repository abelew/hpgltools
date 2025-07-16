#' Sum the reads/gene for multiple sequencing runs of a single condition/batch.
#'
#' On occasion we have multiple technical replicates of a sequencing run.  This
#' can use a column in the experimental design to identify those replicates and
#' sum the counts into a single column in the count tables.
#'
#' Untested as of 2016-12-01, but used in a couple of projects where sequencing
#' runs got repeated.
#'
#' @param expt Experiment class containing the requisite metadata and count tables.
#' @param column Column of the design matrix used to specify which samples are replicates.
#' @return Expt with the concatenated counts, new design matrix, batches, conditions, etc.
#' @seealso [Biobase] [exprs()] [fData()] [pData()] [create_expt()]
#' @examples
#' \dontrun{
#'  compressed <- concatenate_runs(expt)
#' }
#' @export
concatenate_runs <- function(expt, column = "replicate") {
  design <- pData(expt)
  message("The original expressionset has ", nrow(design), " samples.")
  replicates <- levels(as.factor(design[[column]]))
  final_expt <- expt
  final_data <- NULL
  final_design <- NULL
  column_names <- list()
  colors <- list()
  conditions <- list()
  batches <- list()
  samplenames <- list()
  for (rep in replicates) {
    ## expression <- paste0(column, "=='", rep, "'")
    expression <- glue("{column} == '{rep}'")
    tmp_expt <- subset_expt(expt, expression)
    tmp_data <- rowSums(exprs(tmp_expt))
    tmp_design <- pData(tmp_expt)[1, ]
    final_data <- cbind(final_data, tmp_data)
    final_design <- rbind(final_design, tmp_design)
    column_names[[rep]] <- as.character(tmp_design[, "sampleid"])
    colors[[rep]] <- as.character(tmp_expt[["colors"]][1])
    batches[[rep]] <- as.character(pData(tmp_expt)[["batch"]][1])
    conditions[[rep]] <- as.character(pData(tmp_expt)[["condition"]][1])
    samplenames[[rep]] <- paste(conditions[[rep]], batches[[rep]], sep = "-")
    colnames(final_data) <- column_names
  }
  metadata <- new("AnnotatedDataFrame", final_design)
  sampleNames(metadata) <- colnames(final_data)
  feature_data <- new("AnnotatedDataFrame", fData(expt))
  featureNames(feature_data) <- rownames(final_data)
  experiment <- new("ExpressionSet", exprs = final_data,
                    phenoData = metadata, featureData = feature_data)
  final_expt[["expressionset"]] <- experiment
  final_expt[["samples"]] <- final_design
  final_expt[["colors"]] <- as.character(colors)
  final_expt[["samplenames"]] <- as.character(samplenames)
  message("The final expressionset has ", nrow(pData(final_expt)), " samples.")
  return(final_expt)
}
setGeneric("concatenate_runs")

#' Sum the reads/gene for multiple sequencing runs of a single condition/batch.
#'
#' On occasion we have multiple technical replicates of a sequencing run.  This
#' can use a column in the experimental design to identify those replicates and
#' sum the counts into a single column in the count tables.
#'
#' Untested as of 2016-12-01, but used in a couple of projects where sequencing
#' runs got repeated.
#'
#' @param expt Experiment class containing the requisite metadata and count tables.
#' @param column Column of the design matrix used to specify which samples are replicates.
#' @return Expt with the concatenated counts, new design matrix, batches, conditions, etc.
#' @seealso [Biobase] [exprs()] [fData()] [pData()] [create_expt()]
#' @examples
#' \dontrun{
#'  compressed <- concatenate_runs(expt)
#' }
#' @export
setMethod(
  "concatenate_runs", signature(expt = "SummarizedExperiment", column = "character"),
  definition = function(expt, column = "replicate") {
    design <- colData(expt)
    message("The original SE has ", nrow(design), " samples.")
    replicates <- levels(as.factor(design[[column]]))
    final_expt <- expt
    final_data <- NULL
    final_design <- NULL
    column_names <- list()
    colors <- list()
    conditions <- list()
    batches <- list()
    samplenames <- list()
    for (rep in replicates) {
      idx <- design[[column]] == rep
      sub <- expt[, idx]
      tmp_data <- rowSums(assay(sub))
      tmp_design <- colData(sub)[1, ]
      final_data <- cbind(final_data, tmp_data)
      final_design <- rbind(final_design, tmp_design)
      column_names[[rep]] <- as.character(tmp_design[, "sampleid"])
      colors[[rep]] <- as.character(get_colors(sub)[1])
      batches[[rep]] <- as.character(batches(sub)[1])
      conditions[[rep]] <- as.character(conditions(sub)[1])
      samplenames[[rep]] <- paste(conditions[[rep]], batches[[rep]], sep = "-")
      colnames(final_data) <- column_names
    }

    metadata <- new("AnnotatedDataFrame", as.data.frame(final_design))
    sampleNames(metadata) <- colnames(final_data)
    feature_data <- new("AnnotatedDataFrame", as.data.frame(rowData(expt)))
    featureNames(feature_data) <- rownames(final_data)
    experiment <- SummarizedExperiment(assays = final_data,
                                       colData = as.data.frame(final_design),
                                       rowData = as.data.frame(rowData(expt)))
    colors(experiment) <- colors
    message("The final SE has ", nrow(colData(experiment)), " samples.")
    return(experiment)
  })
