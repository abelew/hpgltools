## Time-stamp: <Mon Apr 25 15:07:10 2016 Ashton Trey Belew (abelew@gmail.com)>

## plot_scatter.r: Various scatter plots

#' Steal edgeR's plotBCV() and make it a ggplot2
#' This was written primarily to understand what that function is doing in edgeR.
#'
#' @param data  A dataframe/expt/exprs with count data
#' @return a plot! of the BCV a la ggplot2.
#' @seealso \pkg{edgeR} \link[edgeR]{plotBCV}
#' @examples
#' \dontrun{
#' bcv <- hpgl_bcv_plot(expt)
#' summary(bcv$data)
#' bcv$plot
#' }
#' @export
hpgl_bcv_plot <- function(data) {
    data_class <- class(data)[1]
    if (data_class == "expt") {
        data <- Biobase::exprs(data$expressionset)
    } else if (data_class == "ExpressionSet") {
        data <- Biobase::exprs(data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    data <- edgeR::DGEList(counts=data)
    edisp <- edgeR::estimateDisp(data)
    avg_log_cpm <- edisp$AveLogCPM
    if (is.null(avg_log_cpm)) {
        avg_log_cpm <- edgeR::aveLogCPM(edisp$counts, offset=edgeR::getOffset(edisp))
    }
    disper <- edgeR::getDispersion(edisp)
    if (is.null(disper)) {
        stop("No dispersions to plot")
    }
    if (attr(disper, "type") == "common") {
        disper <- rep(disper, length=length(avg_log_cpm))
    }
    disp_df <- data.frame("A" = avg_log_cpm,
                          "disp" = sqrt(disper))
    fitted_disp <- gplots::lowess(disp_df$A, disp_df$disp, f=0.5)
    f <- stats::approxfun(fitted_disp, rule=2)
    disp_plot <- ggplot2::ggplot(disp_df, ggplot2::aes_string(x="A", y="disp")) +
        ggplot2::geom_point() +
        ggplot2::xlab("Average log(CPM)") +
        ggplot2::ylab("Dispersion of Biological Variance") +
        ##ggplot2::stat_density2d(geom="tile", ggplot2::aes(fill=..density..^0.25), contour=FALSE, show_guide=FALSE) +
        ## ..density.. leads to no visible binding for global variable, but I don't fully understand that notation
        ## I remember looking at it a while ago and being confused
        ggplot2::stat_density2d(geom="tile", ggplot2::aes_string(fill="..density..^0.25"), contour=FALSE, show.legend=FALSE) +
        ggplot2::scale_fill_gradientn(colours=grDevices::colorRampPalette(c("white","black"))(256)) +
        ggplot2::geom_smooth(method="loess") +
        ggplot2::stat_function(fun=f, colour="red") +
        ggplot2::theme(legend.position="none")
    ret <- list("data"=disp_df, "plot"=disp_plot)
    return(ret)
}

#' Make a pretty scatter plot between two sets of numbers with a
#' cheesy distance metric and some statistics of the two sets.
#'
#' The distance metric should be codified and made more intelligent.
#' Currently it creates a dataframe of distances which are absolute
#' distances from each axis, multiplied by each other, summed by axis,
#' then normalized against the maximum.
#'
#' @param df  a dataframe likely containing two columns
#' @param tooltip_data   a df of tooltip information for gvis
#' graphs.
#' @param gvis_filename   a filename to write a fancy html graph.
#' Defaults to NULL in which case the following parameter isn't needed.
#' @param size size of the dots
#' @return a ggplot2 scatter plot.  This plot provides a "bird's eye"
#' view of two data sets.  This plot assumes the two data structures
#' are not correlated, and so it calculates the median/mad of each
#' axis and uses these to calculate a stupid, home-grown distance
#' metric away from both medians.  This distance metric is used to
#' color dots which are presumed the therefore be interesting because
#' they are far from 'normal.'  This will make a fun clicky googleVis
#' graph if requested.
#' @seealso \pkg{ggplot2} \link{hpgl_gvis_scatter} \link[ggplot2]{geom_point}
#' \link{hpgl_linear_scatter}
#' @examples
#' \dontrun{
#'  hpgl_dist_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe,
#'                    gvis_filename="html/fun_scatterplot.html")
#' }
#' @export
hpgl_dist_scatter <- function(df, tooltip_data=NULL, gvis_filename=NULL, size=2) {
    hpgl_env <- environment()
    df <- data.frame(df[, c(1,2)])
    df <- df[complete.cases(df) ,]
    df_columns <- colnames(df)
    df_x_axis <- df_columns[1]
    df_y_axis <- df_columns[2]
    colnames(df) <- c("first","second")
    first_median <- summary(df[, 1])["Median"]
    second_median <- summary(df[, 2])["Median"]
    first_mad <- stats::mad(df[, 1])
    second_mad <- stats::mad(df[, 2])
    mydist <- sillydist(df[, 1], df[, 2], first_median, second_median)
    mydist$x <- abs((mydist[, 1] - first_median) / abs(first_median))
    mydist$y <- abs((mydist[, 2] - second_median) / abs(second_median))
    mydist$x <- mydist$x / max(mydist$x)
    mydist$y <- mydist$y / max(mydist$y)
    mydist$dist <- mydist$x * mydist$y
    mydist$dist <- mydist$dist / max(mydist$dist)
    line_size <- size / 2
    first_vs_second <- ggplot2::ggplot(df, ggplot2::aes_string(x="first", y="second"), environment=hpgl_env) +
        ggplot2::xlab(paste("Expression of", df_x_axis)) +
        ggplot2::ylab(paste("Expression of", df_y_axis)) +
        ggplot2::geom_vline(color="grey", xintercept=(first_median - first_mad), size=line_size) +
        ggplot2::geom_vline(color="grey", xintercept=(first_median + first_mad), size=line_size) +
        ggplot2::geom_vline(color="darkgrey", xintercept=first_median, size=line_size) +
        ggplot2::geom_hline(color="grey", yintercept=(second_median - second_mad), size=line_size) +
        ggplot2::geom_hline(color="grey", yintercept=(second_median + second_mad), size=line_size) +
        ggplot2::geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
        ggplot2::geom_point(colour=grDevices::hsv(mydist$dist, 1, mydist$dist), alpha=0.6, size=size) +
        ggplot2::theme(legend.position="none")
    if (!is.null(gvis_filename)) {
        hpgl_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename)
    }
    return(first_vs_second)
}

#' Make a pretty scatter plot between two sets of numbers with a
#' linear model superimposed and some supporting statistics.
#'
#' @param df  a dataframe likely containing two columns
#' @param tooltip_data   a df of tooltip information for gvis
#' graphs.
#' @param gvis_filename   a filename to write a fancy html graph.
#' @param cormethod   what type of correlation to check?
#' @param size   size of the dots on the plot.
#' @param verbose   be verbose?
#' @param identity   add the identity line?
#' @param loess   add a loess estimation?
#' @param gvis_trendline  add a trendline to the gvis plot?  There are a couple possible types, I think linear is the most common.
#' @param first  first column to plot
#' @param second  second column to plot
#' @param base_url  a base url to add to the plot
#' @param pretty_colors  colors!
#' @return a list including a ggplot2 scatter plot and some
#' histograms.  This plot provides a "bird's eye"
#' view of two data sets.  This plot assumes a (potential) linear
#' correlation between the data, so it calculates the correlation
#' between them.  It then calculates and plots a robust linear model
#' of the data using an 'SMDM' estimator (which I don't remember how
#' to describe, just that the document I was reading said it is good).
#' The median/mad of each axis is calculated and plotted as well.  The
#' distance from the linear model is finally used to color the dots on
#' the plot.  Histograms of each axis are plotted separately and then
#' together under a single cdf to allow tests of distribution
#' similarity.  This will make a fun clicky googleVis graph if
#' requested.
#' @seealso \link[robust]{lmRob} \link[stats]{weights} \link{hpgl_histogram}
#' @examples
#' \dontrun{
#'  hpgl_linear_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe,
#'                      gvis_filename="html/fun_scatterplot.html")
#' }
#' @export
hpgl_linear_scatter <- function(df, tooltip_data=NULL, gvis_filename=NULL, cormethod="pearson",
                                size=2, verbose=FALSE, loess=FALSE, identity=FALSE, gvis_trendline=NULL,
                                first=NULL, second=NULL, base_url=NULL, pretty_colors=TRUE) {
    hpgl_env <- environment()
    df <- data.frame(df[,c(1,2)])
    df <- df[complete.cases(df),]
    correlation <- cor.test(df[,1], df[,2], method=cormethod, exact=FALSE)
    df_columns <- colnames(df)
    df_x_axis <- df_columns[1]
    df_y_axis <- df_columns[2]
    colnames(df) <- c("first","second")
    model_test <- try(robustbase::lmrob(formula=second ~ first, data=df, method="SMDM"), silent=TRUE)
    if (class(model_test) == "try-error") {
        model_test <- try(lm(formula=second ~ first, data=df), silent=TRUE)
    }
    if (class(model_test) == "try-error") {
        model_test <- try(glm(formula=second ~ first, data=df), silent=TRUE)
    }
    if (class(model_test) == "try-error") {
        message("Could not perform a linear modelling of the data.")
        message("Going to perform a scatter plot without linear model.")
        plot <- hpgl_scatter(df)
        ret <- list(data=df, scatter=plot)
        return(ret)
    } else {
        linear_model <- model_test
    }
    linear_model <- try(robustbase::lmrob(formula=second ~ first, data=df, method="SMDM"))
    linear_model_summary <- summary(linear_model)
    linear_model_rsq <- linear_model_summary$r.squared
    linear_model_weights <- stats::weights(linear_model, type="robustness", na.action=NULL)
    linear_model_intercept <- stats::coef(linear_model_summary)[1]
    linear_model_slope <- stats::coef(linear_model_summary)[2]
    first_median <- summary(df$first)["Median"]
    second_median <- summary(df$second)["Median"]
    first_mad <- stats::mad(df$first, na.rm=TRUE)
    second_mad <- stats::mad(df$second, na.rm=TRUE)
    line_size <- size / 2
    first_vs_second <- ggplot2::ggplot(df, ggplot2::aes_string(x="first", y="second"), environment=hpgl_env) +
        ggplot2::xlab(paste("Expression of", df_x_axis)) +
        ggplot2::ylab(paste("Expression of", df_y_axis)) +
        ggplot2::geom_vline(color="grey", xintercept=(first_median - first_mad), size=line_size) +
        ggplot2::geom_vline(color="grey", xintercept=(first_median + first_mad), size=line_size) +
        ggplot2::geom_hline(color="grey", yintercept=(second_median - second_mad), size=line_size) +
        ggplot2::geom_hline(color="grey", yintercept=(second_median + second_mad), size=line_size) +
        ggplot2::geom_hline(color="darkgrey", yintercept=second_median, size=line_size) +
        ggplot2::geom_vline(color="darkgrey", xintercept=first_median, size=line_size) +
        ggplot2::geom_abline(colour="grey", slope=linear_model_slope, intercept=linear_model_intercept, size=line_size)
    if (isTRUE(pretty_colors)) {
        first_vs_second <- first_vs_second +
            ggplot2::geom_point(size=size, alpha=0.4,
                                colour=grDevices::hsv(linear_model_weights * 9/20,
                                                      linear_model_weights/20 + 19/20,
                                                      (1.0 - linear_model_weights)))
    } else {
        first_vs_second <- first_vs_second +
            ggplot2::geom_point(colour="black", size=size, alpha=0.4)
    }
    if (loess == TRUE) {
        first_vs_second <- first_vs_second +
            ggplot2::geom_smooth(method="loess")
    }
    if (identity == TRUE) {
        first_vs_second <- first_vs_second +
            ggplot2::geom_abline(colour="darkgreen", slope=1, intercept=0, size=1)
    }
    first_vs_second <- first_vs_second +
        ggplot2::theme(legend.position="none") +
        ggplot2::theme_bw()

    if (!is.null(gvis_filename)) {
        if (verbose) {
            message("Generating an interactive graph.")
        }
        hpgl_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename,
                          trendline=gvis_trendline, base_url=base_url)
    }
    if (!is.null(first) & !is.null(second)) {
        colnames(df) <- c(first, second)
    } else if (!is.null(first)) {
        colnames(df) <- c(first, "second")
    } else if (!is.null(second)) {
        colnames(df) <- c("first", second)
    }
    x_histogram <- hpgl_histogram(data.frame(df[, 1]), verbose=verbose, fillcolor="lightblue", color="blue")
    y_histogram <- hpgl_histogram(data.frame(df[, 2]), verbose=verbose, fillcolor="pink", color="red")
    both_histogram <- hpgl_multihistogram(df, verbose=verbose)
    plots <- list(data=df, scatter=first_vs_second, x_histogram=x_histogram,
                  y_histogram=y_histogram, both_histogram=both_histogram,
                  correlation=correlation, lm_model=linear_model, lm_summary=linear_model_summary,
                  lm_weights=linear_model_weights, lm_rsq=linear_model_rsq,
                  first_median=first_median, first_mad=first_mad,
                  second_median=second_median, second_mad=second_mad)
    if (verbose) {
        message(sprintf("Calculating correlation between the axes using:", cormethod))
        message(correlation)
        message("Calculating linear model between the axes")
        message(linear_model_summary)
        message("Generating histogram of the x axis.")
        message("Generating histogram of the y axis.")
        message("Generating a histogram comparing the axes.")
    }
    return(plots)
}

#' Make a pretty MA plot from the output of voom/limma/eBayes/toptable.
#'
#' @param counts  df of linear-modelling, normalized counts by sample-type,
#'        which is to say the output from voom/voomMod/hpgl_voom().
#' @param de_genes  df from toptable or its friends containing p-values.
#' @param pval_cutoff   a cutoff defining significant from not.
#' @param alpha   how transparent to make the dots.
#' @param logfc_cutoff a fold change cutoff
#' @param pval the name of the pvalue column to use for cutoffs
#' @param size   how big are the dots?
#' @param tooltip_data   a df of tooltip information for gvis
#' @param gvis_filename   a filename to write a fancy html graph.
#'        graphs.
#' @param ... more options for you
#' @return a ggplot2 MA scatter plot.  This is defined as the rowmeans
#' of the normalized counts by type across all sample types on the
#' x-axis, and the log fold change between conditions on the y-axis.
#' Dots are colored depending on if they are 'significant.'  This will
#' make a fun clicky googleVis graph if requested.
#' @seealso \link{hpgl_gvis_ma_plot} \link[limma]{toptable}
#' \link[limma]{voom} \link{hpgl_voom}
#' \link[limma]{lmFit} \link[limma]{makeContrasts}
#' \link[limma]{contrasts.fit}
#' @examples
#' \dontrun{
#' hpgl_ma_plot(voomed_data, toptable_data, gvis_filename="html/fun_ma_plot.html")
#' ## Currently this assumes that a variant of toptable was used which
#' ## gives adjusted p-values.  This is not always the case and I should
#' ## check for that, but I have not yet.
#' }
#' @export
hpgl_ma_plot <- function(counts, de_genes, pval_cutoff=0.05, alpha=0.5, logfc_cutoff=1, pval="adjpval",
                         size=2, tooltip_data=NULL, gvis_filename=NULL, ...) {
    hpgl_env <- environment()
    if (pval == "adjpval") {
        pval_column <- "adj.P.Val"
        aes_color <- "(adjpval <= pval_cutoff)"
    } else {
        pval_column <- "P.Value"
        aes_color <- "(pval <= pval_cutoff)"
    }
    df <- data.frame("avg" = rowMeans(counts[rownames(de_genes),]),
                     "logfc" = de_genes[["logFC"]],
                     "pval" = de_genes[["P.Value"]],
                     "adjpval" = de_genes[[pval_column]])
    df[["adjpval"]] <- as.numeric(format(df[["adjpval"]], scientific=FALSE))
    df[["pval"]] <- as.numeric(format(df[["pval"]], scientific=FALSE))
    df$state <- ifelse(df[["adjpval"]] > pval_cutoff, "pinsig",
                ifelse(df[["adjpval"]] <= pval_cutoff & df[["logfc"]] >= logfc_cutoff, "upsig",
                ifelse(df[["adjpval"]] <= pval_cutoff & df[["logfc"]] <= (-1 * logfc_cutoff), "downsig", "fcinsig")))
    num_pinsig <- sum(df[["state"]] == "pinsig")
    num_upsig <- sum(df[["state"]] == "upsig")
    num_downsig <- sum(df[["state"]] == "downsig")
    num_fcinsig <- sum(df[["state"]] == "fcinsig")
    plt <- ggplot2::ggplot(df, ggplot2::aes_string(x="avg", y="logfc", color=aes_color),
                           environment=hpgl_env) +
        ggplot2::geom_hline(yintercept=c((logfc_cutoff * -1), logfc_cutoff), color="red", size=(size / 2)) +
        ggplot2::geom_point(stat="identity", size=size, alpha=alpha, ggplot2::aes_string(shape="as.factor(state)", fill=aes_color)) +
        ggplot2::scale_shape_manual(name="state", values=c(21,22,23,24),
                                    labels=c(
                                        paste0("Down Sig.: ", num_downsig),
                                        paste0("FC Insig.: ", num_fcinsig),
                                        paste0("P Insig.: ", num_pinsig),
                                        paste0("Up Sig.: ", num_upsig)),
                                    guide=ggplot2::guide_legend(override.aes=ggplot2::aes(size=3, fill="grey"))) +
        ggplot2::scale_color_manual(values=c("FALSE"="darkred","TRUE"="darkblue")) +
        ggplot2::scale_fill_manual(values=c("FALSE"="darkred","TRUE"="darkblue")) +
        ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(size=3))) +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=-90)) +
        ggplot2::xlab("Average Count (Millions of Reads)") +
        ggplot2::ylab("log fold change") +
        ggplot2::theme_bw()
    if (!is.null(gvis_filename)) {
        hpgl_gvis_ma_plot(counts, de_genes, tooltip_data=tooltip_data, filename=gvis_filename, ...)
    }
    return(plt)
}
## Consider using these options for the kind of pretty graph Eva likes.
##ggplot(mydata) + aes(x=x, y=y) + scale_x_log10() + scale_y_log10() +
##+   stat_density2d(geom="tile", aes(fill=..density..^0.25), contour=FALSE) +
##+   scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256))

#'   Make a ggplot graph of the number of non-zero genes by sample.
#' Made by Ramzi Temanni <temanni at umd dot edu>
#'
#' @param data an expt, expressionset, or dataframe.
#' @param design   a design matrix.
#' @param colors   a color scheme.
#' @param labels   how do you want to label the graph?
#'   'fancy' will use directlabels() to try to match the labels with the positions without overlapping
#'   anything else will just stick them on a 45' offset next to the graphed point
#' @param title   add a title?
#' @param ... rawr
#' @return a ggplot2 plot of the number of non-zero genes with respect to each library's CPM
#' @seealso \link[ggplot2]{geom_point} \link[directlabels]{geom_dl}
#' @examples
#' \dontrun{
#'  nonzero_plot = hpgl_nonzero(expt=expt)
#'  nonzero_plot  ## ooo pretty
#' }
#' @export
hpgl_nonzero <- function(data, design=NULL, colors=NULL, labels=NULL, title=NULL, ...) {
    hpgl_env <- environment()
    names <- NULL
    data_class <- class(data)[1]
    if (data_class == "expt") {
        design <- data$design
        colors <- data$colors
        names <- data$samplenames
        data <- Biobase::exprs(data$expressionset)
    } else if (data_class == "ExpressionSet") {
        data <- Biobase::exprs(data)
    } else if (data_class == "matrix" | data_class == "data.frame") {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }

    if (is.null(labels)) {
        if (is.null(names)) {
            labels <- colnames(data)
        } else {
            labels <- names
        }
    } else if (labels[1] == "boring") {
        if (is.null(names)) {
            labels <- colnames(data)
        } else {
            labels <- names
        }
    }

    shapes <- as.integer(as.factor(design$batch))
    non_zero <- data.frame("id"=colnames(data),
                           "nonzero_genes"=colSums(data >= 1),
                           "cpm"=colSums(data) * 1e-6,
                           "condition"=design$condition,
                           "batch"=design$batch)
    non_zero_plot <- ggplot2::ggplot(data=non_zero, ggplot2::aes_string(x="cpm", y="nonzero_genes"), environment=hpgl_env, fill=colors, shape=shapes) +
        ## geom_point(stat="identity", size=3, colour=hpgl_colors, pch=21) +
        ggplot2::geom_point(ggplot2::aes_string(fill="colors"), colour="black", pch=21, stat="identity", size=3) +
        ggplot2::scale_fill_manual(name="Condition", values=levels(as.factor(colors)), labels=levels(as.factor(design$condition))) +
        ggplot2::ylab("Number of non-zero genes observed.") +
        ggplot2::xlab("Observed CPM") +
        ggplot2::theme_bw()
    if (!is.null(labels)) {
        if (labels[[1]] == "fancy") {
            non_zero_plot <- non_zero_plot +
                directlabels::geom_dl(ggplot2::aes_string(label="labels"),
                                      method="smart.grid", colour=colors)
        } else {
            non_zero_plot <- non_zero_plot +
                ggplot2::geom_text(ggplot2::aes_string(x="cpm", y="nonzero_genes", label="labels"),
                                   angle=45, size=4, vjust=2)
        }
    }
    if (!is.null(title)) {
        non_zero_plot <- non_zero_plot + ggplot2::ggtitle(title)
    }
    non_zero_plot <- non_zero_plot +
        ggplot2::theme(axis.ticks=ggplot2::element_blank(), axis.text.x=ggplot2::element_text(angle=90))
    return(non_zero_plot)
}


#'   Plot all pairwise MA plots in an experiment.
#'
#' Use affy's ma.plot() on every pair of columns in a data set to help
#' diagnose problematic samples.
#'
#' @param data an expt expressionset or data frame
#' @param log   is the data in log format?
#' @param ... more options are good
#' @return a list of affy::maplots
#' @seealso \link[affy]{ma.plot}
#' @examples
#' \dontrun{
#'  ma_plots = hpgl_pairwise_ma(expt=some_expt)
#' }
#' @export
hpgl_pairwise_ma <- function(data, log=NULL, ...) {
    ##require.auto('affy')
    data_class <- class(data)[1]
    if (data_class == "expt") {
        design <- data[["design"]]
        colors <- data[["colors"]]
        data <- Biobase::exprs(data[["expressionset"]])
    } else if (data_class == 'ExpressionSet') {
        data <- Biobase::exprs(data)
    } else if (data_class == 'matrix' | data_class == 'data.frame') {
        data <- as.data.frame(data)  ## some functions prefer matrix, so I am keeping this explicit for the moment
    } else {
        stop("This function currently only understands classes of type: expt, ExpressionSet, data.frame, and matrix.")
    }
    plot_list <- list()
    for (c in 1:(length(colnames(data)) - 1)) {
        nextc <- c + 1
        for (d in nextc:length(colnames(data))) {
            first <- as.numeric(data[, c])
            second <- as.numeric(data[, d])
            if (max(first) > 1000) {
                if (is.null(log)) {
                    message("I suspect you want to set log=TRUE for this.")
                    message("In fact, I am so sure, I am doing it now.")
                    message("If I am wrong, set log=FALSE, but I'm not.")
                    log <- TRUE
                }
            } else if (max(first) < 80) {
                if (!is.null(log)) {
                    message("I suspect you want to set log=FALSE for this.")
                    message("In fact, I am so  sure, I am doing it now.")
                    message("If I am wrong, set log=TRUE.")
                    log <- FALSE
                }
            }
            firstname <- colnames(data)[c]
            secondname <- colnames(data)[d]
            name <- paste0(firstname, "_", secondname)
            if (isTRUE(log)) {
                first <- log2(first + 1.0)
                second <- log2(second + 1.0)
            }
            m <- first - second
            a <- (first + second) / 2
            affy::ma.plot(A=a, M=m, plot.method="smoothScatter", show.statistics=TRUE, add.loess=TRUE)
            title(paste0("MA of ", firstname, " vs ", secondname))
            plot_list[[name]] = grDevices::recordPlot()
        }
    }
    return(plot_list)
}

#'   Make a pretty scatter plot between two sets of numbers.
#'
#' @param df  a dataframe likely containing two columns
#' @param gvis_filename   a filename to write a fancy html graph.
#' @param tooltip_data   a df of tooltip information for gvis
#' @param size   the size of the dots on the graph.
#' @param color   color of the dots on the graph.
#' @return a ggplot2 scatter plot.
#' @seealso \link{hpgl_gvis_scatter} \link[ggplot2]{geom_point}
#' \link{hpgl_linear_scatter}
#' @examples
#' \dontrun{
#' hpgl_scatter(lotsofnumbers_intwo_columns, tooltip_data=tooltip_dataframe,
#'              gvis_filename="html/fun_scatterplot.html")
#' }
#' @export
hpgl_scatter <- function(df, tooltip_data=NULL, color="black", gvis_filename=NULL, size=2) {
    hpgl_env <- environment()
    df <- data.frame(df[,c(1,2)])
    df <- df[complete.cases(df),]
    df_columns <- colnames(df)
    df_x_axis <- df_columns[1]
    df_y_axis <- df_columns[2]
    colnames(df) <- c("first","second")
    first_vs_second <- ggplot2::ggplot(df, ggplot2::aes_string(x="first", y="second"), environment=hpgl_env) +
        ggplot2::xlab(paste("Expression of", df_x_axis)) +
        ggplot2::ylab(paste("Expression of", df_y_axis)) +
        ggplot2::geom_point(colour=color, alpha=0.6, size=size) +
        ggplot2::theme(legend.position="none")
    if (!is.null(gvis_filename)) {
        hpgl_gvis_scatter(df, tooltip_data=tooltip_data, filename=gvis_filename)
    }
    return(first_vs_second)
}

#'   Make a pretty Volcano plot!
#'
#' @param toptable_data  a dataframe from limma's toptable which
#' includes log(fold change) and an adjusted p-value.
#' @param p_cutoff   a cutoff defining significant from not.
#' @param fc_cutoff   a cutoff defining the minimum/maximum fold change
#' for interesting.  This is log, so I went with +/- 0.8 mostly
#' arbitrarily as the default.
#' @param alpha   how transparent to make the dots.
#' @param size   how big are the dots?
#' @param gvis_filename  a filename to write a fancy html graph.
#' @param tooltip_data   a df of tooltip information for gvis.
#' @param ...  I love parameters!
#' @return a ggplot2 MA scatter plot.  This is defined as the
#' -log10(p-value) with respect to log(fold change).  The cutoff
#' values are delineated with lines and mark the boundaries between
#' 'significant' and not.  This will make a fun clicky googleVis graph
#' if requested.
#' @seealso \link{hpgl_gvis_ma_plot} \link[limma]{toptable}
#' \link[limma]{voom} \link{hpgl_voom} \link[limma]{lmFit}
#' \link[limma]{makeContrasts} \link[limma]{contrasts.fit}
#' @examples
#' \dontrun{
#'  hpgl_volcano_plot(toptable_data, gvis_filename="html/fun_ma_plot.html")
#' ## Currently this assumes that a variant of toptable was used which
#' ## gives adjusted p-values.  This is not always the case and I should
#' ## check for that, but I have not yet.
#' }
#' @export
hpgl_volcano_plot <- function(toptable_data, tooltip_data=NULL, gvis_filename=NULL,
                              fc_cutoff=0.8, p_cutoff=0.05, size=2, alpha=0.6, ...) {
    hpgl_env <- environment()
    low_vert_line <- 0.0 - fc_cutoff
    horiz_line <- -1 * log10(p_cutoff)
    toptable_data$modified_p <- -1 * log10(toptable_data$P.Value)
    plt <- ggplot2::ggplot(toptable_data,
                           ggplot2::aes_string(x="logFC", y="modified_p", color="(P.Value <= p_cutoff)"),
                           environment=hpgl_env) +
        ggplot2::geom_hline(yintercept=horiz_line, color="black", size=size) +
        ggplot2::geom_vline(xintercept=fc_cutoff, color="black", size=size) +
        ggplot2::geom_vline(xintercept=low_vert_line, color="black", size=size) +
        ggplot2::geom_point(stat="identity", size=size, alpha=alpha) +
        ## theme(axis.text.x=element_text(angle=-90)) +
        ggplot2::xlab("log fold change") +
        ggplot2::ylab("-log10(adjusted p value)") +
        ggplot2::theme(legend.position="none")
    if (!is.null(gvis_filename)) {
        hpgl_gvis_volcano_plot(toptable_data, fc_cutoff=fc_cutoff, p_cutoff=p_cutoff, tooltip_data=tooltip_data, filename=gvis_filename)
    }
    return(plt)
}