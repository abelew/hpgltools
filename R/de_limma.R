## de_limma.r: Simplify/standardize inputs/outputs when doing differential
## expression with limma.  The ideas here are identical to de_edger.r and
## de_deseq.r.  It is worth noting that limma was the first tool I learned to
## use and therefore some of the code in this file is some of the very first R
## code I learned and therefore potentially... bad.

#' A minor change to limma's voom with quality weights to attempt to address some corner cases.
#'
#' This copies the logic employed in hpgl_voom().  I suspect one should not use it.
#'
#' @param data Some data!
#' @param fun_model A model for voom() and arrayWeights()
#' @param libsize Library sizes passed to voom().
#' @param normalize.method Passed to voom()
#' @param plot Do the plot of mean variance?
#' @param span yes
#' @param var.design maybe
#' @param method kitty!
#' @param maxiter 50 is good
#' @param tol I have no tolerance.
#' @param trace no trace for you.
#' @param replace.weights  Replace the weights?
#' @param col yay columns!
#' @param ... more arguments!
#' @return a voom return
#' @seealso [limma::voom()]
#' @examples
#' \dontrun{
#' ## No seriously, dont run this, I think it is wiser to use the functions
#' ## provided by limma. But this provides a place to test stuff out.
#'  voom_result <- hpgl_voomweighted(dataset, model)
#' }
#' @export
hpgl_voomweighted <- function(data, fun_model, libsize = NULL, normalize.method = "none",
                              plot = TRUE, span = 0.5, var.design = NULL, method = "genebygene",
                              maxiter = 50, tol = 1E-10, trace = FALSE, replace.weights = TRUE,
                              col = NULL, ...) {

  if (isTRUE(plot)) {
    oldpar <- par(mfrow = c(1, 2))
    on.exit(par(oldpar))
  }
  v1 <- hpgl_voom(data, model = fun_model, libsize = libsize,
                  normalize.method = normalize.method,
                  plot = TRUE, span = span, ...)
  aw <- try(limma::arrayWeights(v1, design = fun_model, method = method, maxiter = maxiter,
                                tol = tol, var.design = var.design))
  if (class(aw) == "try-error") {
    message("arrayWeights failed, returning the voom result.")
    return(v1)
  }
  v <- hpgl_voom(data, model = fun_model, weights = aw, libsize = libsize,
                 normalize.method = normalize.method, plot = TRUE, span = span, ...)
  aw <- limma::arrayWeights(v, design = fun_model, method = method, maxiter = maxiter,
                            tol = tol, trace = trace, var.design = var.design)
  wts <- limma::asMatrixWeights(aw, dim(v)) * v[["weights"]]
  attr(wts, "arrayweights") <- NULL
  if (plot) {
    tmp_file <- tmpmd5file(pattern = "voom", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    barplot(aw, names = 1:length(aw), main = "Sample-specific weights",
            ylab = "Weight", xlab = "Sample", col = col)
    abline(h = 1, col = 2, lty = 2)
    voom_barplot <- grDevices::recordPlot()
    plotted <- dev.off()
    removed <- file.remove(tmp_file)
  }
  if (replace.weights) {
    v[["weights"]] <- wts
    v[["sample.weights"]] <- aw
    v[["barplot"]] <- voom_barplot
    v[["first_iter"]] <- v1
    return(v)
  } else {
    return(wts)
  }
}

#' A slight modification of limma's voom().
#'
#' Estimate mean-variance relationship between samples and generate
#' 'observational-level weights' in preparation for linear modeling RNAseq data.
#' This particular implementation was primarily scabbed from cbcbSEQ, but
#' changes the mean-variance plot slightly and attempts to handle corner cases
#' where the sample design is confounded by setting the coefficient to 1 for
#' those samples rather than throwing an unhelpful error.  Also, the Elist
#' output gets a 'plot' slot which contains the plot rather than just printing
#' it.
#'
#' @param dataframe Dataframe of sample counts which have been normalized and
#'  log transformed.
#' @param model Experimental model defining batches/conditions/etc.
#' @param libsize Size of the libraries (usually provided by edgeR).
#' @param normalize.method Normalization method used in voom().
#' @param span The span used in voom().
#' @param stupid Cheat when the resulting matrix is not solvable?
#' @param logged Is the input data is known to be logged?
#' @param converted Is the input data is known to be cpm converted?
#' @param ... Extra arguments are passed to arglist.
#' @return EList containing the following information:
#'  E = The normalized data
#'  weights = The weights of said data
#'  design = The resulting design
#'  lib.size = The size in pseudocounts of the library
#'  plot = A ggplot of the mean/variance trend with a blue loess fit and red trend fit
#' @seealso [limma::voom()]
#' @examples
#' \dontrun{
#'  funkytown = hpgl_voom(samples, model)
#' }
#' @export
hpgl_voom <- function(dataframe, model = NULL, libsize = NULL,
                      normalize.method = "none", span = 0.5,
                      stupid = FALSE, logged = FALSE, converted = FALSE, ...) {
  arglist <- list(...)
  ## Going to attempt to as closely as possible dovetail the original implementation.
  ## I think at this point, my implementation is the same as the original with the exception
  ## of a couple of tests to check that the data is not fubar and I think my plot is prettier.
  counts <- dataframe
  out <- list()
  if (is(counts, "DGEList")) {
    out[["genes"]] <- counts[["genes"]]
    out[["targets"]] <- counts[["samples"]]
    if (is.null(model) &
          diff(range(as.numeric(counts[["sample"]][["group"]]))) > 0) {
      model <- model.matrix(~group, data = counts[["samples"]])
    }
    if (is.null(libsize)) {
      ## libsize <- with(counts[["samples"]], libsize * norm.factors)
      ## This is a bit confusing.
      libsize <- with(counts[["samples"]], counts[["libsize"]] * counts[["norm.factors"]])
    }
    counts <- counts[["counts"]]
  } else {
    isExpressionSet <- sm(is(counts, "ExpressionSet"))
    if (isExpressionSet) {
      if (length(fData(counts))) {
        out[["genes"]] <- fData(counts)
      }
      if (length(pData(counts))) {
        out[["targets"]] <- pData(counts)
      }
      counts <- exprs(counts)
    } else {
      counts <- as.matrix(counts)
    }
  }
  if (is.null(model)) {
    model <- matrix(1, ncol(counts), 1)
    rownames(model) <- colnames(counts)
    colnames(model) <- "GrandMean"
  }
  if (is.null(libsize)) {
    libsize <- colSums(dataframe, na.rm = TRUE)
  }
  if (converted == "cpm") {
    converted <- TRUE
  }
  if (!isTRUE(converted)) {
    message("The voom input was not cpm, converting now.")
    posed <- t(dataframe + 0.5)
    dataframe <- t(posed / (libsize + 1) * 1e+06)
    ##y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1000000)) ## from voom()
  }
  if (logged == "log2") {
    logged <- TRUE
  }
  if (isTRUE(logged)) {
    if (max(dataframe, na.rm = TRUE) > 1000) {
      warning("This data appears to not be logged, the lmfit will do weird things.")
    }
  } else {
    if (max(dataframe, na.rm = TRUE) < 200) {
      warning("This data says it was not logged, but the maximum counts seem small.")
      warning("If it really was log2 transformed, ",
              "then we are about to double-log it and that would be very bad.")
    }
    message("The voom input was not log2, transforming now.")
    dataframe <- log2(dataframe)
  }
  dataframe <- as.matrix(dataframe)
  dataframe <- limma::normalizeBetweenArrays(dataframe, method = normalize.method)
  linear_fit <- limma::lmFit(dataframe, model, method = "ls")
  if (is.null(linear_fit[["Amean"]])) {
    linear_fit[["Amean"]] <- rowMeans(dataframe, na.rm = TRUE)
  }
  sx <- linear_fit[["Amean"]] + mean(log2(libsize + 1)) - log2(1e+06)
  sy <- sqrt(linear_fit[["sigma"]])

  ## First remove NAs
  sx_nas <- is.na(sx)
  sy_nas <- is.na(sy)
  or_nas <- sx_nas | sy_nas
  ## Then look for all-zero entries.
  all_zero <- rowSums(dataframe, na.rm = TRUE) == 0
  or_nas_zero <- or_nas | all_zero
  sx <- sx[!or_nas_zero]
  sy <- sy[!or_nas_zero]

  fitted <- gplots::lowess(sx, sy, f = 0.5)
  f <- stats::approxfun(fitted, rule = 2)
  mean_var_df <- data.frame(mean = sx, var = sy)
  mean_var_plot <- ggplot(mean_var_df,
                          aes(x = .data[["mean"]], y = .data[["var"]])) +
    ggplot2::geom_point() +
    ggplot2::xlab("Log2(count size + 0.5)") +
    ggplot2::ylab("Square root of the standard deviation.") +
    ggplot2::stat_density2d(geom = "tile", aes(fill = ggplot2::after_stat(density^0.25)),
                            contour = FALSE, show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = grDevices::colorRampPalette(c("white", "black"))(256)) +
    ggplot2::geom_smooth(method = "loess") +
    ggplot2::stat_function(fun = f, colour = "red") +
    ggplot2::theme(legend.position = "none")
  if (is.null(linear_fit[["rank"]])) {
    message("Some samples cannot be balanced across the experimental design.")
    if (isTRUE(stupid)) {
      ## I think this is telling me I have confounded data, and so
      ## for those replicates I will have no usable coefficients, so
      ## I say set them to 1 and leave them alone.
      linear_fit[["coefficients"]][is.na(linear_fit[["coefficients"]])] <- 1
      fitted.values <- linear_fit[["coefficients"]] %*%
        t(linear_fit[["design"]])
    }
  } else if (linear_fit[["rank"]] < ncol(linear_fit[["design"]])) {
    j <- linear_fit[["pivot"]][1:linear_fit[["rank"]]]
    fitted.values <- linear_fit[["coefficients"]][, j, drop = FALSE] %*%
      t(linear_fit[["design"]][, j, drop = FALSE])
  } else {
    fitted.values <- linear_fit[["coefficients"]] %*%
      t(linear_fit[["design"]])
  }
  fitted.cpm <- 2 ^ fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (libsize + 1.0))
  fitted.logcount <- log2(fitted.count)
  w <- 1 / f(fitted.logcount) ^ 4
  dim(w) <- dim(fitted.logcount)
  rownames(w) <- rownames(dataframe)
  colnames(w) <- colnames(dataframe)
  out[["E"]] <- dataframe
  out[["weights"]] <- w
  out[["design"]] <- model
  out[["lib.size"]] <- libsize
  out[["plot"]] <- mean_var_plot
  new("EList", out)
}

#' Set up a model matrix and set of contrasts for pairwise comparisons using
#' voom/limma.
#'
#' Creates the set of all possible contrasts and performs them using voom/limma.
#'
#' @param input Dataframe/vector or expt class containing count tables,
#'  normalization state, etc.
#' @param model_fstring Formula string describing the statistical model of interest.
#' @param null_fstring Formula string describing the null model.
#' @param model_svs Matrix of surrogates or method to seek them.
#' @param filter Filter the data before seeking SVs?
#' @param extra_contrasts Some extra contrasts to add to the list.
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param annot_df Data frame for annotations.
#' @param libsize I've recently figured out that libsize is far more important
#'  than I previously realized.  Play with it here.
#' @param adjust Use this p-value adjustment.
#' @param which_voom Try out different invocations of voom.
#' @param limma_method And different invocations of limma itself.
#' @param limma_robust Pass along the robust args for limma?
#' @param voom_norm Use a specific normalization for voom?
#' @param limma_trend Include a trendline in the limma plot?
#' @param force Force data which may not be appropriate for limma into it?
#' @param keep_underscore Sanitize away model underscores?
#' @param num_surrogates Explicit number of surrogates or way to find them.
#' @param keepers Choose a set of contrasts instead of all.
#' @param ... Use the elipsis parameter to feed options to write_limma().
#' @return List including the following information:
#'  macb = the mashing together of condition/batch so you can look at it
#'  macb_model = The result of calling model.matrix(~0 + macb)
#'  macb_fit = The result of calling lmFit(data, macb_model)
#'  voom_result = The result from voom()
#'  voom_design = The design from voom (redundant from voom_result, but convenient)
#'  macb_table = A table of the number of times each condition/batch pairing happens
#'  cond_table = A table of the number of times each condition appears (the
#'   denominator for the identities)
#'  batch_table = How many times each batch appears
#'  identities = The list of strings defining each condition by itself
#'  all_pairwise = The list of strings defining all the pairwise contrasts
#'  contrast_string = The string making up the makeContrasts() call
#'  pairwise_fits = The result from calling contrasts.fit()
#'  pairwise_comparisons = The result from eBayes()
#'  limma_result = The result from calling write_limma()
#' @seealso [limma] [Biobase] [deseq_pairwise()] [edger_pairwise()] [basic_pairwise()]
#'  DOI:10.1093/nar/gkv007
#' @examples
#' \dontrun{
#'  pretend <- limma_pairwise(expt)
#' }
#' @export
limma_pairwise <- function(input = NULL, model_fstring = "~ 0 + condition + batch",
                           null_fstring = "~", model_svs = NULL, filter = FALSE,
                           extra_contrasts = NULL, annot_df = NULL,
                           libsize = NULL, adjust = "BH",
                           which_voom = "limma", limma_method = "ls",
                           limma_robust = FALSE, voom_norm = "quantile",
                           limma_trend = TRUE, force = FALSE, keep_underscore = TRUE,
                           num_surrogates = "be", keepers = NULL,
                           ...) {
  arglist <- list(...)
  ## This is used in the invocation of a voom() implementation for normalization.
  ## This is for the eBayes() call.

  current <- state(input)
  filteredp <- current[["filter"]]
  if ((is.null(filteredp) || filteredp == "raw") &&
        isTRUE(filter)) {
    input <- sm(normalize(input, filter = filter))
  }
  input <- sanitize_expt(input, keep_underscore = keep_underscore)
  input_data <- choose_limma_dataset(input, force = force, which_voom = which_voom)
  count_mtrx <- as.matrix(input_data[["data"]])
  fctrs <- get_formula_factors(model_fstring)
  condition_column <- fctrs[["factors"]][1]
  design <- pData(input)
  conditions <- droplevels(as.factor(design[[condition_column]]))
  batches <- droplevels(as.factor(design[["batch"]]))
  libsize <- libsize(input)

  ## The following small piece of logic is intended to handle situations where we use
  ## tximport for limma (kallisto/sailfish/salmon).
  if (is.null(input[["tximport"]])) {
    count_mtrx <- input_data[["data"]]
  } else {
    data <- edgeR::DGEList(input[["tximport"]][["scaled"]][["counts"]])
    count_mtrx <- edgeR::calcNormFactors(data)
  }

  if (is.null(libsize)) {
    libsize <- libsize(input)
  }

  condition_table <- table(conditions)
  batch_table <- table(batches)
  mesg("This limma pairwise comparison should compare across:")
  print(condition_table)

  loggedp <- current[["transform"]]
  if (is.null(loggedp)) {
    mesg("I don't know if this data is logged, testing if it is integer.")
    if (is.integer(data)) {
      loggedp <- FALSE
    } else {
      loggedp <- TRUE
    }
  } else {
    if (grepl(pattern = "log", x = loggedp)) {
      loggedp <- TRUE
    } else {
      loggedp <- FALSE
    }
  }

  convertedp <- current[["conversion"]]
  if (is.null(convertedp)) {
    message("I cannot determine if this data has been converted, assuming no.")
    convertedp <- FALSE
  } else {
    if (convertedp == "raw") {
      convertedp <- FALSE
    } else {
      convertedp <- TRUE
    }
  }

  appended_fstring <- model_fstring
  if ("character" %in% class(model_svs)) {
    model_params <- adjuster_expt_svs(input, model_fstring = model_fstring,
                                      null_fstring = null_fstring,
                                      model_svs = model_svs,
                                      num_surrogates = num_surrogates,
                                      ...)
    estimate_type <- model_svs
    model_svs <- model_params[["model_adjust"]]
    null_model <- model_params[["null_model"]]
    appended_fstring <- model_params[["appended_fstring"]]
    design <- pData(model_params[["modified_input"]])
  }
  model_mtrx <- model.matrix(as.formula(appended_fstring), data = design)

  ## voom() it, taking into account whether the data has been log2 transformed.
  ## Leaving the following here for the moment, but I think it will no longer be needed.
  ## Instead, I am checking the data state before passing it to this function with the
  ## choose_limma_dataset() call above.
  na_sum <- sum(is.na(count_mtrx))
  if (na_sum > 0) {
    ## My version of voom can handle nas
    which_voom <- "hpgl"
  }
  voom_result <- NULL
  voom_plot <- NULL
  if (which_voom == "hpgl_weighted") {
    mesg("Limma step 1/5: running hpgl_voomweighted(), switch with the argument 'which_voom'.")
    voom_result <- hpgl_voomweighted(
      count_mtrx, model_mtrx, libsize = libsize, voom_norm = voom_norm,
      span = 0.5, var.design = NULL, method = "genebygene",
      maxiter = 50, tol = 1E-10, trace = FALSE, replace.weights = TRUE, col = NULL,
      logged = loggedp, converted = convertedp)
    voom_plot <- voom_result[["plot"]]
  } else if (which_voom == "hpgl") {
    mesg("Limma step 1/5: running hpgl_voom(), switch with the argument 'which_voom'.")
    voom_result <- hpgl_voom(
      count_mtrx, model_mtrx, libsize = libsize,
      logged = loggedp, converted = convertedp)
    voom_plot <- voom_result[["plot"]]
  } else if (which_voom == "limma_weighted") {
    mesg("Limma step 1/5: running limma::voomWithQualityWeights(), ",
         "switch with the argument 'which_voom'.")
    tmp_file <- tmpmd5file(pattern = "voom", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    voom_result <- try(limma::voomWithQualityWeights(
      counts = count_mtrx, design = model_mtrx, lib.size = libsize,
      normalize.method = voom_norm, plot = TRUE, span = 0.5,
      var.design = NULL, method = "genebygene", maxiter = 50,
      tol = 1E-10, trace = FALSE, replace.weights = TRUE, col = NULL))
    if (class(voom_result) == "try-error") {
      mesg("voomWithQualityWeights failed, falling back to voom.")
      voom_result <- limma::voom(
        counts = count_mtrx, design = model_mtrx, lib.size = libsize,
        normalize.method = voom_norm, span = 0.5, plot = TRUE, save.plot = TRUE)
    }
    voom_plot <- grDevices::recordPlot()
    plotted <- dev.off()
    removed <- file.remove(tmp_file)
  } else if (which_voom == "limma") {
    mesg("Limma step 1/5: running limma::voom(), switch with the argument 'which_voom'.")
    mesg("Using normalize.method = ", voom_norm, " for voom.")
    ## Note to self, the defaults are span = 0.5, plot = FALSE, save.plot = FALSE,
    ## normalize.method = "none", lib.size = NULL, design = NULL
    tmp_file <- tmpmd5file(pattern = "voom", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    voom_result <- limma::voom(
      counts = count_mtrx, design = model_mtrx, lib.size = libsize,
      normalize.method = voom_norm, span = 0.5, plot = TRUE, save.plot = TRUE)
    voom_plot <- grDevices::recordPlot()
    plotted <- dev.off()
    removed <- file.remove(tmp_file)
  } else if (which_voom == "cpm") {
    ## Reyy reminded me today that one does not necessarily need voom, but
    ## logcpm might be sufficient when the data distributions are nice and
    ## consistent.
    mesg("Limma step 1/5: using edgeR::cpm(), switch with the argument 'which_voom'.")
    voom_result <- edgeR::cpm(count_mtrx, log = TRUE, prior.count = 3)
  } else if (which_voom == "none") {
    ## Then this is microarray-ish data.
    mesg("Assuming this data is similar to a microarray and not performing voom.")
    voom_result <- count_mtrx
  } else {
    mesg("Limma step 1/5: running limma::voom(), switch with the argument 'which_voom'.")
    mesg("Using normalize.method = ", voom_norm, " for voom.")
    ## Note to self, the defaults are span = 0.5, plot = FALSE, save.plot = FALSE,
    ## normalize.method = "none", lib.size = NULL, design = NULL
    tmp_file <- tmpmd5file(pattern = "voom", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    voom_result <- limma::voom(
      counts = count_mtrx, design = model_mtrx, lib.size = libsize,
      normalize.method = voom_norm, span = 0.5, plot = TRUE, save.plot = TRUE)
    voom_plot <- grDevices::recordPlot()
    plotted <- dev.off()
    removed <- file.remove(tmp_file)
  }
  one_replicate <- FALSE
  if (is.null(voom_result)) {
    ## Apparently voom returns null where there is only 1 replicate.
    message("voom returned null, I am not sure what will happen.")
    one_replicate <- TRUE
    voom_result <- count_mtrx
  }

  ## Do the lmFit() using this model
  pairwise_fits <- NULL
  identity_fits <- NULL
  mesg("Limma step 2/5: running lmFit with method: ", limma_method, ".")
  fitted_data <- limma::lmFit(object = voom_result,
                              design = model_mtrx,
                              method = limma_method)
  all_tables <- NULL
  if (is.null(fctrs[["cellmeans_intercept"]])) {
    message("Limma step 3/5: making and fitting contrasts with a model intercept. (~ factors)")
    contrasts <- "nointercept"
    all_pairwise_contrasts <- NULL
    contrast_string <- "no intercept done"
    all_pairwise <- NULL
    pairwise_fits <- fitted_data
    identity_contrasts <- NULL
    identities <- NULL
    identity_fits <- fitted_data
    mesg("Limma step 4/5: Running eBayes with robust = ",
         limma_robust, " and trend = ", limma_trend, ".")
    all_pairwise_comparisons <- limma::eBayes(fitted_data,
                                              robust = limma_robust,
                                              trend = limma_trend)
    all_identity_comparisons <- NULL
    mesg("Limma step 5/5: Creating tables with an intercept.")
    pairwise_results <- make_limma_tables(fit = all_pairwise_comparisons, adjust = "BH",
                                          n = 0, coef = NULL, annot_df = NULL, intercept = TRUE)
    limma_tables <- pairwise_results[["contrasts"]]
    contrasts_performed <- names(limma_tables)
    limma_identities <- pairwise_results[["identities"]]
  } else {
    mesg("Limma step 3/5: making and fitting contrasts with no intercept. (~ 0 + factors)")
    contrasts <- make_pairwise_contrasts(
      model = model_mtrx, conditions = conditions, contrast_factor = condition_column,
      extra_contrasts = extra_contrasts, keepers = keepers, keep_underscore = keep_underscore,
      do_identities = FALSE)
    all_pairwise_contrasts <- contrasts[["all_pairwise_contrasts"]]
    contrast_string <- contrasts[["contrast_string"]]
    all_pairwise <- contrasts[["all_pairwise"]]
    ## Once all that is done, perform the fit
    ## This will first provide the relative abundances of each condition
    ## followed by the set of all pairwise comparisons.
    pairwise_fits <- limma::contrasts.fit(fit = fitted_data, contrasts = all_pairwise_contrasts)
    identity_contrasts <- sm(make_pairwise_contrasts(model = model_mtrx, conditions = conditions,
                                                     contrast_factor = condition_column,
                                                     do_identities = TRUE, do_pairwise = FALSE,
                                                     keep_underscore = keep_underscore))
    identities <- identity_contrasts[["all_pairwise_contrasts"]]
    identity_fits <- limma::contrasts.fit(fit = fitted_data, contrasts = identities)
    mesg("Limma step 4/5: Running eBayes with robust = ",
            limma_robust, " and trend = ", limma_trend, ".")
    if (isTRUE(one_replicate)) {
      all_pairwise_comparisons <- pairwise_fits[["coefficients"]]
      all_identity_comparisons <- identity_fits[["coefficients"]]
    } else {
      all_pairwise_comparisons <- limma::eBayes(pairwise_fits,
                                                robust = limma_robust,
                                                trend = limma_trend)
      all_identity_comparisons <- limma::eBayes(identity_fits,
                                                robust = limma_robust,
                                                trend = limma_trend)
    }
    mesg("Limma step 5/5: Creating tables without an intercept.")
    pairwise_results <- make_limma_tables(fit = all_pairwise_comparisons, adjust = adjust,
                                          n = 0, coef = NULL, annot_df = NULL)
    limma_tables <- pairwise_results[["contrasts"]]
    identity_results <- make_limma_tables(fit = all_identity_comparisons, adjust = adjust,
                                          n = 0, coef = NULL, annot_df = NULL)
    limma_identities <- identity_results[["identities"]]

    contrasts_performed <- names(limma_tables)
  }

  retlist <- list(
    "all_pairwise" = all_pairwise,
    "all_tables" = limma_tables,
    "batches" = batches,
    "batches_table" = batch_table,
    "conditions" = conditions,
    "condition_table" = condition_table,
    "contrast_string" = contrast_string,
    "contrasts_performed" = contrasts_performed,
    "design" = design,
    "dispersion_plot" = voom_plot,
    "fit" = fitted_data,
    "identities" = identities,
    "identity_tables" = limma_identities,
    "identity_comparisons" = all_identity_comparisons,
    "input_data" = input,
    "method" = "limma",
    "model" = model_mtrx,
    "model_fstring" = model_fstring,
    "pairwise_comparisons" = all_pairwise_comparisons,
    "single_table" = all_tables,
    "voom_result" = voom_result)
  class(retlist) <- c("limma_pairwise", "list")
  if (!is.null(arglist[["limma_excel"]])) {
    retlist[["limma_excel"]] <- write_limma(retlist, excel = arglist[["limma_excel"]])
  }
  return(retlist)
}

#' @export
print.limma_pairwise <- function(x, ...) {
  summary_string <- glue("The results from the limma pairwise analysis.")
  message(summary_string)
  return(invisible(x))
}

#' Writes out the results of a limma search using toptable().
#'
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the toptable() output in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param fit Result from lmFit()/eBayes()
#' @param adjust Pvalue adjustment chosen.
#' @param n Number of entries to report, 0 says do them all.
#' @param coef Which coefficients/contrasts to report, NULL says do them all.
#' @param annot_df Optional data frame including annotation information to
#'  include with the tables.
#' @param intercept Intercept model?
#' @return List of data frames comprising the toptable output for each
#'  coefficient, I also added a qvalue entry to these toptable() outputs.
#' @seealso [limma] [write_xlsx()]
#' @examples
#' \dontrun{
#'  finished_comparison = eBayes(limma_output)
#'  table = make_limma_tables(finished_comparison, adjust = "fdr")
#' }
make_limma_tables <- function(fit = NULL, adjust = "BH", n = 0, coef = NULL,
                              annot_df = NULL, intercept = FALSE) {
  ## Figure out the number of genes if not provided
  if (n == 0) {
    n <- nrow(fit[["coefficients"]])
  }

  ## If specific contrast(s) is/are not requested, get them all.
  if (is.null(coef)) {
    if (isTRUE(intercept)) {
      coef <- colnames(fit[["coefficients"]])
      coef <- coef[2:length(coef)]
    } else {
      coef <- colnames(fit[["contrasts"]])
    }
  } else {
    coef <- as.character(coef)
  }
  return_identities <- list()
  return_data <- list()
  end <- length(coef)
  data_tables <- list()
  if (isTRUE(intercept)) {
    ## If we do have an intercept model, then we get the data
    ## in a slightly different fashion.
    for (d in seq_len(ncol(fit[["coefficients"]]))) {
      data_table <-  limma::topTable(fit, adjust.method = adjust,
                                     n = n, coef = d, sort.by = "logFC")

      for (column in seq_len(ncol(data_table))) {
        data_table[[column]] <- signif(x = as.numeric(data_table[[column]]), digits = 4)
      }
      if (!is.null(annot_df)) {
        data_table <- merge(data_table, annot_df, by.x = "row.names", by.y = "row.names")
      }

      if (d == 1) {
        return_identities[[1]] <- data_table
      } else {
        comparison <- colnames(fit[["coefficients"]])[d]
        return_data[[comparison]] <- data_table
      }
    }
  } else {
    ## If we do not have an intercept (~ 0 + ...)
    ## Then extract the coefficients and identities separately.
    for (it in seq_len(end)) {
      comparison <- coef[it]
      comp_name <- strsplit(x = comparison, split = " = ")[[1]][1]
      mesg("Limma step 6/6: ", it, "/", end, ": Creating table: ",
           comp_name, ".  Adjust = ", adjust)
      data_tables[[comp_name]] <- limma::topTable(
        fit, adjust.method = adjust,
        n = n, coef = comparison, sort.by = "logFC")
    }
    ## Take a moment to prettily format the numbers in the tables
    ## and fill in the identity table.
    for (e in seq_along(data_tables)) {
      comparison <- coef[e]
      comp_name <- strsplit(x = comparison, split = " = ")[[1]][1]
      table <- data_tables[[e]]
      for (column in seq_len(ncol(table))) {
        table[[column]] <- signif(x = as.numeric(table[[column]]), digits = 4)
      }
      if (!is.null(annot_df)) {
        table <- merge(table, annot_df, by.x = "row.names", by.y = "row.names")
      }
      if (grepl(pattern = "_vs_", x = comparison)) {
        return_data[[comp_name]] <- table
      } else {
        return_identities[[comp_name]] <- table
      }
    }
  } ## End checking for an intercept/nointercept model.

  retlist <- list(
    "identities" = return_identities,
    "contrasts" = return_data)
  return(retlist)
}

#' Writes out the results of a limma search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from limma and friends.
#'
#' @param data Output from limma_pairwise()
#' @param ... Options for writing the xlsx file.
#' @seealso [write_de_table()]
#' @examples
#' \dontrun{
#'  finished_comparison = limma_pairwise(expressionset)
#'  data_list = write_limma(finished_comparison)
#' }
#' @export
write_limma <- function(data, ...) {
  result <- write_de_table(data, type = "limma", ...)
  return(result)
}

## EOF
