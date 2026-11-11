## grange.R: Various functions to assist when working with ranged datatypes.

#' @include 01_hpgltools.R
NULL

#' I was thinking that once I start using this binned GR, I will start wanting to add
#' stuff to it, so let us start with GC content.
#'
#' Ok, this works but is annoyingly slow.  Perhaps making it an apply() will help?
add_gc_content <- function(bsg, gr, pattern = "GC") {
  elements <- length(gr)
  gc_vector <- rep(0, elements)
  for (b in seq_along(gr)) {
    bin <- gr[b]
    seq <- Biostrings::getSeq(bsg, seqnames(bin), start = BiocGenerics::start(bin), end = BiocGenerics::end(bin))
    gc_vector[b] <- Biostrings::letterFrequency(seq, pattern, as.prob = TRUE)
    mesg(b, " of ", elements, ".")
  }
  mcols(gr)[["gc"]] <- as.numeric(gc_vector)
  return(gr)
}

#' Next I want every bin to include the following mcols:
#'   1. CDS which may be NA when the bin does not overlap with any gene
#'   2. fivep which may be NA when the bin is entirely in a CDS, but otherwise points to the
#' downstream gene.
#'   3. threep which may be NA when the bin is entirely in a CDS, but otherwise points to the
#' upstream gene.
crossref_bins_genes <- function(bins, annotation_gff, annotation_type = "mRNA",
                                exclude_pattern = "SCAF", window = 1000, output = NULL) {
  bins = binned_gr
  annotation_gff = "~/libraries/genome/gff/lpanamensis_v68.gff"
  annotation_gr <- gff2gr(annotation_gff)
  if (!is.null(annotation_type)) {
    annotation_idx <- mcols(annotation_gr)[["type"]] == annotation_type
    annotation_gr <- annotation_gr[annotation_idx]
  }

  if (!is.null(exclude_pattern)) {
    uninteresting <- grepl(x = seqnames(annotation_gr), pattern = exclude_pattern)
    annotation_gr <- annotation_gr[!uninteresting]
  }

  annotation_df <- mcols(annotation_gr)
  ## We cannot use the various join functions if the mcols have start/end in them; so take that out
  ## after we are done using it.
  original_annotation <- annotation_df
  annotation_df[["gene_start"]] <- IRanges::start(annotation_gr)
  annotation_df[["gene_end"]] <- IRanges::end(annotation_gr)
  annotation_df[["gene_strand"]] <- GenomicRanges::strand(annotation_gr)
  mcols(annotation_gr) <- annotation_df
  
  mcols(bins)[["fivep"]] <- NA
  mcols(bins)[["CDS"]] <- NA
  mcols(bins)[["threep"]] <- NA
  final_mcols <- mcols(bins)
    
  joined <- plyranges::join_overlap_intersect(bins, annotation_gr)
  mcols(joined) <- mcols(joined)[, c("ID", "gene_start", "gene_end", "gene_strand", "description")]
  bin_df <- mcols(bins)
  gene_df <- mcols(joined)
  merged_df <- merge(bin_df, gene_df, all.x = TRUE, by = "row.names")
  rownames(merged_df) <- merged_df[["Row.names"]]
  merged_df[["Row.names"]] <- NULL
  final_mcols <- merged_df[rownames(bin_df), ]
  final_mcols[["CDS"]] <- final_mcols[["ID"]]
  final_mcols[["ID"]] <- NULL

  ## Now get rid of the extra start/end before we get yelled at.
  mcols(annotation_gr) <- original_annotation

  upstream <- plyranges::flank_upstream(annotation_gr, width = window)
  joined_upstream <- plyranges::join_overlap_intersect(bins, upstream)
  mcols(joined_upstream) <- mcols(joined_upstream)[, c("ID")]
  colnames(mcols(joined_upstream)) <- "ID"
  bin_df <- mcols(bins)
  gene_df <- mcols(joined_upstream)
  merged_df <- merge(bin_df, gene_df, all.x = TRUE, by = "row.names")
  rownames(merged_df) <- merged_df[["Row.names"]]
  merged_df[["Row.names"]] <- NULL
  merged_df <- merged_df[rownames(bin_df), ]
  final_mcols[["fivep"]] <- merged_df[["ID"]]

  downstream <- flank_downstream(annotation_gr, width = window)
  joined_downstream <- plyranges::join_overlap_intersect(bins, downstream)
  mcols(joined_downstream) <- mcols(joined_downstream)[, c("ID")]
  colnames(mcols(joined_downstream)) <- "ID"
  bin_df <- mcols(bins)
  gene_df <- mcols(joined_downstream)
  merged_df <- merge(bin_df, gene_df, all.x = TRUE, by = "row.names")
  rownames(merged_df) <- merged_df[["Row.names"]]
  merged_df[["Row.names"]] <- NULL
  merged_df <- merged_df[rownames(bin_df), ]
  final_mcols[["threep"]] <- merged_df[["ID"]]

  mcols(bins) <- final_mcols
  if (!is.null(output)) {
    plyranges::write_gff3(x = bins, file = output)
  }
  return(bins)
}

#' Create a windowed GRange from a genome.  Assume a BSGenome, but have methods
#' to handle fasta etc.
#'
#' We would like to have a binned summarized experiment of expression and variants, so I figured a
#' quick and dirty way to handle it was to generate a fasta and recount our reads/variants using it.
#'
#' I want to have a little fun with this and add a mcol with GC content per bin.  I want this as
#' practice for adding more difficult things later.
#'
#' @param genome Input genome, by default bsgenome.
#' @param window_size what it says on the tin.
#' @param add_tails Add an extra bin at the end of each chromosome with the remaining portion.
#' @param output Output gff or other file(s) (right now just gff).
#' @param exclude_pattern Regex used to exclude regions which are not interesting. (scaffolds)
#' @return List with a new GRange and perhaps some useful information, maybe # bins/chromsoome?
#' @export
genome_windows <- function(input_genome, window_size = 1000, add_tails = TRUE, output = NULL,
                           exclude_pattern = "SCAF") {
  chr_info <- seqinfo(input_genome)
  all_chromosomes <- as.data.frame(chr_info)
  all_gr <- GRanges(chr_info)
  if (!is.null(exclude_pattern)) {
    uninteresting <- grepl(x = rownames(all_chromosomes), pattern = exclude_pattern)
    wanted_names <- rownames(all_chromosomes)[!uninteresting]
    all_chromosomes <- all_chromosomes[wanted_names, ]
    chr_info <- chr_info[wanted_names]
    all_gr <- GRanges(chr_info)
  }
  binned_gr <- GRanges()
  for (g in seq_along(all_gr)) {
    this_gr <- all_gr[g]
    tile_n <- lengths(this_gr) / window_size
    fragments <- GenomicRanges::tile(this_gr, n = tile_n)
    ## I am going to do a round-trip through data.table in order to ensure that I can get the tails.
    fragment_dt <- as.data.table(fragments[[1]])
    fragment_df <- as.data.frame(fragment_dt)
    if (isTRUE(add_tails)) {
      ## I need to add one more for the end of the chromosome, tile returns a list, so grab the first element.
      this_region <- fragments[[1]]
      final_fragment <- length(this_region)
      start_of_final_endpoint <- BiocGenerics::end(this_region[final_fragment]) + 1
      end_of_final_endpoint <- BiocGenerics::end(this_gr)
      final_length <- end_of_final_endpoint - start_of_final_endpoint
      name <- as.character(seqnames(this_region))[1]
      last_row <- c(name, start_of_final_endpoint, end_of_final_endpoint, final_length, "*")
      fragment_df <- rbind(as.data.frame(fragment_dt), last_row)
      fragment_dt <- as.data.table(fragment_df)
    }
    new_names <- paste0(fragment_dt[["seqnames"]], "_", fragment_dt[["start"]], "_", fragment_dt[["end"]])
    rownames(fragment_df) <- new_names
    final_gr <- GRanges(fragment_df)
    binned_gr <- c(binned_gr, final_gr)
  } ## End iterating over every chromosome.
  
  if (!is.null(output)) {
    written <- plyranges::write_gff3(x = binned_gr, file = output)
  }
  return(binned_gr)
}

#' Extract a grange from gene x to gene y
#'
#' @param gr Input grange.
#' @param from Starting point.
#' @param to Ending point.
#' @param column What mcol column to use to extract IDs.
#' @param type_column Use this feature type.
#' @param type Type of entries from which to create the granges.
#' @param padding_start Add this many genes as padding to the grange before the start.
#' @param padding_end Add this many genes as padding to the grange after the end.
#' @export
gr_from_to <- function(gr, from, to, column = "gene_id", type_column = "type",
                       type = NULL, padding_start = 1, padding_end = 1) {
  meta <- as.data.frame(mcols(gr))
  meta_na <- is.na(meta[[column]])
  if (sum(meta_na) > 0) {
    message("Dropping ", sum(meta_na), " rows with NA from column ", column, ".")
    gr <- gr[!meta_na, ]
    meta <- mcols(gr)
  }

  if (!is.null(type)) {
    mesg("Extracting elements of type ", type, " from the granges.")
    wanted_idx <- meta[[type_column]] == type
    gr <- gr[wanted_idx, ]
    meta <- mcols(gr)
  }

  beginning_id <- meta[[column]] == from
  if (sum(beginning_id, na.rm = TRUE) != 1) {
    stop("This only works if we get one gene for the beginning.")
  }
  beginning_num <- which(beginning_id) - padding_start
  beginning_gr <- gr[beginning_num, ]
  end_id <- meta[[column]] == to
  end_num <- which(end_id) + padding_end
  if (sum(end_id) != 1) {
    stop("This only works if we get one gene for the end.")
  }
  end_gr <- gr[end_num, ]
  wanted_begin <- GenomicRanges::start(beginning_gr)
  wanted_end <- GenomicRanges::end(end_gr)
  annotation_starts <- GenomicRanges::start(gr)
  annotation_ends <- GenomicRanges::end(gr)
  ## Take the set of genes where our wanted beginning position is
  ## larger than the end of the previous gene
  ## and the wanted end is smaller than the start of the next
  ## Or stick with start vs start and end vs end?
  wanted_idx <- annotation_starts >= wanted_begin &
    annotation_ends <= wanted_end
  mesg("This range results in ", sum(wanted_idx), " genes in the subset.")
  region <- gr[wanted_idx, ]
  mcols(region)[, "orientation"] <- TRUE
  minus_idx <- GenomicRanges::strand(region) == "-" | GenomicRanges::strand(region) == -1
  mcols(region)[minus_idx, "orientation"] <- FALSE
  mcols(region)[, "gene_name"] <- mcols(region)[, "gene_id"]
  mcols(region)[, "gene_biotype"] <- "gene"
  return(region)
}

#' Stealing the LoadTrackFile function from ggcoverage and making it more helpful.
#'
#' @param se Input summarized experiment.
#' @param track_column Column name with track files.
#' @param region_string String describing the chromosome(s).
#' @param name_column colData column containing the desired row names.
#' @param gene_name Specific gene name(s)
#' @param padding Specific number of nucleotides to pad the coverage region.
#' @param bin_size What it says on the tin.
#' @param transform Perform a log transformation on the coverage?
#' @param convert Perform a cpmish conversion on the coverage?
#' @param norm Normalize the coverage?
#' @param type_column colData column describing the data types.
#' @param group_column colData column describing the conditions of interest.
#' @param cores Parallelize this?
#' @param extend Expand the resulting grange?
#' @export
load_se_tracks <- function(se, track_column = "deeptools_coverage", region_string = NULL,
                           name_column = "Parent", gene_name = NULL, padding = 100,
                           bin_size = 10, transform = NULL, convert = "cpm", norm = NULL,
                           type_column = "media", group_column = "media_type",
                           cores = NULL, extend = NULL) {
  norm_factor <- libsize_factor(se)
  metadata <- colData(se)
  coverage_gr <- NULL
  track_files <- colData(se)[[track_column]]
  names(track_files) <- rownames(metadata)
  format <- tools::file_ext(track_files[1])
  cores <- 1

  if (is.null(region_string)) {
    mesg("No 'region' specified; extracting coverage for an example range\n(<=100,000 bases, first annotated sequence)")
    if (format == "bam") {
      seqnames <- Rsamtools::scanBamHeader(track_files[1]) |>
        lapply(function(x) x$targets) |>
        unname() |>
        unlist()
      coverage_gr <- GenomicRanges::GRanges(seqnames = names(seqnames[1]),
                                            IRanges::IRanges(start = 1, end = min(1e+05, seqnames[1])))
    } else if (format %in% c("wig", "bw", "bedgraph")) {
      coverage_gr <- range(rtracklayer::import(track_files[1]))
      seqnames <- as.character(Seqinfo::seqnames(coverage_gr))
      if (GenomicRanges::width(coverage_gr) <= 1e+05) {
        coverage_gr <- GenomicRanges::resize(coverage_gr, width = 1e+05)
      }
    }
    mesg("Coverage extracted from sequence/chromosome: ", names(seqnames[1]))
  } else {
    mesg("Extracting coverage for the region ", region_string)
    coverage_gr <- ggcoverage_prepare(coverage_gr, region_string = region_string,
                                      gene_name = gene_name, name_column = name_column,
                                      extend = extend)
  }
  if (is.null(cores)) {
    cores <- parallel::detectCores()
  }
  track_list <- list()
  if (format %in% c("wig", "bw", "bedgraph")) {
    if (bin_size == 1) {
      stop("To visualize single nucleotide resolution, please use bam file!")
    } else {
      if (cores == 1) {
        for (t in seq_along(track_files)) {
          sample_name <- names(track_files)[t]
          track_file <- as.character(track_files[t])
          track_gr <- rtracklayer::import(track_file, which = coverage_gr)
          track_df <- as.data.frame(track_gr)
          track_df[["TrackFile"]] <- track_file
          track_df[["sampleid"]] <- sample_name
          track_list[[sample_name]] <- track_df
        }
      } else {
        BiocParallel::register(BiocParallel::MulticoreParam(workers = cores), default = TRUE)
        track_list <- BiocParallel::bplapply(
          track_files, BPPARAM = BiocParallel::MulticoreParam(),
          FUN = ggcoverage_import_bw, coverage_gr, metadata)
      }
    }
  } else if (format == "bam") {
    if (cores == 1) {
      lapply(track_files, "ggcoverage" %:::% "index_bam")
    } else {
      BiocParallel::register(BiocParallel::MulticoreParam(workers = cores),
                             default = TRUE)
      BiocParallel::bplapply(track_files, BPPARAM = BiocParallel::MulticoreParam(),
                             FUN = "ggcoverage" %:::% "index_bam")
    }
    if (bin_size == 1) {
      if (cores == 1) {
        track_list <- lapply(track_files, "ggcoverage" %:::% "single_nuc_cov",
                             bin_size)
      } else {
        track_list <- BiocParallel::bplapply(
          track_files, BPPARAM = BiocParallel::MulticoreParam(),
          FUN = "ggcoverage" %:::% "single_nuc_cov",
          bin_size)
      }
    } else {
      if (norm == "None") {
        message("Calculating coverage with GenomicAlignments when 'norm = None'")
        if (cores == 1) {
          track_list <- lapply(track_files,
                               FUN = "ggcoverage" %:::% "import_bam_ga", coverage_gr, bin_size)
        } else {
          track_list <- BiocParallel::bplapply(
            track_files, BPPARAM = BiocParallel::MulticoreParam(),
            FUN = "ggcoverage" %:::% "import_bam_ga", coverage_gr, bin_size)
        }
      } else {
        message("Calculate coverage with bamCoverage when 'norm != None'")
        if (is.null(bamcoverage.path)) {
          bamcoverage.path <- Sys.which("bamCoverage")
          if (bamcoverage.path == "") {
            stop("Can not find bamCoverage automatically, please specify 'bamcoverage.path'")
          }
        }
        if (cores == 1) {
          bc_extra_parameters <- NULL
          track_list <- lapply(track_files, "ggcoverage" %:::% "bam_coverage",
                               bamcoverage.path, bin_size, norm,
                               bc_extra_parameters, coverage_gr)
        } else {
          track_list <- BiocParallel::bplapply(
            track_files, BPPARAM = BiocParallel::MulticoreParam(), FUN = "ggcoverage" %:::% "bam_coverage",
            bamcoverage.path, bin_size, norm, bc_extra_parameters, coverage_gr)
        }
      }
    }
  } else if (format == "txt") {
    if (bin_size == 1) {
      stop("To visualize single nucleotide, please use bam file!")
    } else {
      if (cores == 1) {
        track_list <- lapply(track_files, "ggcoverage" %:::% "import_txt")
      } else {
        BiocParallel::register(BiocParallel::MulticoreParam(workers = cores),
                               default = TRUE)
        track_list <- BiocParallel::bplapply(
          track_files, BPPARAM = BiocParallel::MulticoreParam(), FUN = "ggcoverage" %:::% "import_txt")
      }
    }
  }

  convert <- "libsize"
  ## do a CPMish normalization before rbind
  if (is.null(convert) || isFALSE(convert)) {
    convert <- "none"
  }
  if (convert == "libsize") {
    mesg("Performing a libsize normalization of the coverage.")
    for (num in seq_along(track_list)) {
      sample_name <- names(track_list)[num]
      track_list[[sample_name]][["score"]] <- track_list[[sample_name]][["score"]] / norm_factor[sample_name]
    }
  }
  if (is.null(transform) || isFALSE(transform)) {
    transform <- "none"
  }
  if (transform == "log2") {
    mesg("log2 transforming the coverage.")
    for (num in seq_along(track_list)) {
      track_list[[num]][["score"]] <- log2(track_list[[num]][["score"]] + 1)
    }
  }

  all_tracks_df <- do.call(rbind, track_list)
  if ("character" %in% class(metadata)) {
    used_metadata <- extract_metadata(metadata)
  } else if (!is.null(metadata)) {
    used_metadata <- metadata
  } else {
    message("No metadata provided, returning coverage as is.")
    used_metadata <- NULL
  }
  if (is.null(used_metadata)) {
    all_tracks_df[["Type"]] <- all_tracks_df[["TrackFile"]]
    all_tracks_df[["Group"]] <- all_tracks_df[["TrackFile"]]
    all_tracks_df[["TrackFile"]] <- NULL
  } else {
    message("Merging the track information to the metadata, this may take a while.")
    all_tracks_df <- merge(all_tracks_df, used_metadata, by.x = "sampleid", by.y = "row.names")
    if (is.null(type_column)) {
      warning("The type column is NULL, setting it to the track file.")
      all_tracks_df[["Type"]] <- all_tracks_df[["TrackFile"]]
    } else {
      all_tracks_df[["Type"]] <- all_tracks_df[[type_column]]
    }
    if (is.null(group_column)) {
      warning("The group column is NULL, setting it to the track file.")
      all_tracks_df[["Group"]] <- all_tracks_df[["TrackFile"]]
    } else {
      all_tracks_df[["Group"]] <- all_tracks_df[[group_column]]
    }
  }
  all_tracks_df[["Group"]] <- as.character(all_tracks_df[["Group"]])
  all_tracks_df[["Type"]] <- as.character(all_tracks_df[["Type"]])
  return(all_tracks_df)
}

#' Copied the prepare function from ggcoverage and made some small changes.
#' Note, that extend is in nucleotides while padding is in genes.
#'
#' @param gr Input grange
#' @param region_string String describing the chromosome(s) of interest.
#' @param gene_name Specific gene of interest.
#' @param wanted_type mcols type of interest.
#' @param name_column mcols column containing the gene names.
#' @param extend Number of nucleotides to expand the resulting region.
ggcoverage_prepare <- function(gr, region_string = NULL, gene_name = "HNRNPC",
                               wanted_type = "gene",
                               name_column = "gene_name", extend = 2000) {
  if (!is.null(region_string)) {
    region_split <- unlist(strsplit(x = region_string, split = ":"))
    region_chr <- region_split[1]
    region_start_end <- unlist(strsplit(x = region_split[2], split = "-"))
    if (length(region_start_end) == 1) {
      region_start <- as.numeric(gsub(pattern = ",", replacement = "",
                                      x = region_start_end[1]))
      region_end <- region_start
    } else if (length(region_start_end) == 2) {
      region_start <- as.numeric(gsub(pattern = ",", replacement = "",
                                      x = region_start_end[1]))
      region_end <- as.numeric(gsub(pattern = ",", replacement = "",
                                    x = region_start_end[2]))
    }
  } else {
    type <- NULL
    gr_info <- gr |>
      as.data.frame() |>
      dplyr::filter(type == wanted_type)
    wanted_rows <- gr_info[, name_column] == gene_name
    used_info <- gr_info[wanted_rows, ]
    region_chr <- as.character(used_info[["seqnames"]])
    region_start <- used_info[["start"]]
    region_end <- used_info[["end"]]
  }

  if (!is.null(extend)) {
    region_start <- region_start - extend
    if (region_start[1] < 1) {
      region_start <- 1
    }
    region_end <- region_end[length(region_end)] + extend
  }
  if (region_start == region_end) {
    message("The start position is same as end position, and the extension is 0. Automatically increase by 1!")
    region_end <- region_start + 1
  }
  gr <- GenomicRanges::GRanges(seqnames = region_chr,
                               ranges = IRanges::IRanges(region_start, region_end))
  return(gr)
}

# track_list <- lapply(track_files, ggcoverage_import_bw, coverage_gr, metadata)

ggcoverage_import_bw <- function(x, gr, sample_id) {
  sample_name <- names(x)
  message("TESTME: ", sample_name, ".")
  track_df <- as.data.frame(rtracklayer::import(x, which = gr))
  track_df[["TrackFile"]] <- as.character(x)
  track_df[["sampleid"]] <- names(x)
  return(track_df)
}
