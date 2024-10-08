% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_edger.R
\name{edger_pairwise}
\alias{edger_pairwise}
\title{Set up a model matrix and set of contrasts to do pairwise comparisons using EdgeR.}
\usage{
edger_pairwise(
  input = NULL,
  conditions = NULL,
  batches = NULL,
  model_cond = TRUE,
  model_batch = TRUE,
  model_intercept = FALSE,
  alt_model = NULL,
  extra_contrasts = NULL,
  annot_df = NULL,
  force = FALSE,
  keepers = NULL,
  edger_method = "long",
  keep_underscore = FALSE,
  ...
)
}
\arguments{
\item{input}{Dataframe/vector or expt class containing data, normalization
state, etc.}

\item{conditions}{Factor of conditions in the experiment.}

\item{batches}{Factor of batches in the experiment.}

\item{model_cond}{Include condition in the experimental model?}

\item{model_batch}{Include batch in the model?  In most cases this is a good
thing(tm).}

\item{model_intercept}{Use an intercept containing model?}

\item{alt_model}{Alternate experimental model to use?}

\item{extra_contrasts}{Add some extra contrasts to add to the list of
pairwise contrasts. This can be pretty neat, lets say one has conditions
A,B,C,D,E and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
de_vs_cb = (E-D)-(C-B),"}

\item{annot_df}{Annotation information to the data tables?}

\item{force}{Force edgeR to accept inputs which it should not have to deal with.}

\item{keepers}{Ask for a specific set of contrasts instead of all.}

\item{edger_method}{I found a couple/few ways of doing edger in the manual,
choose with this.}

\item{...}{The elipsis parameter is fed to write_edger() at the end.}
}
\value{
List including the following information:
 contrasts = The string representation of the contrasts performed.
 lrt = A list of the results from calling glmLRT(), one for each contrast.
 contrast_list = The list of each call to makeContrasts()
 I do this to avoid running into the limit on # of contrasts addressable by topTags()
 all_tables = a list of tables for the contrasts performed.
}
\description{
This function performs the set of possible pairwise comparisons using EdgeR.
}
\details{
Like the other _pairwise() functions, this attempts to perform all pairwise
contrasts in the provided data set.  The details are of course slightly
different when using EdgeR.  Thus, this uses the function
choose_binom_dataset() to try to ensure that the incoming data is appropriate
for EdgeR (if one normalized the data, it will attempt to revert to raw
counts, for example). It continues on to extract the conditions and batches
in the data, choose an appropriate experimental model, and run the EdgeR
analyses as described in the manual.  It defaults to using an experimental
batch factor, but will accept a string like 'sva' instead, in which case it
will use sva to estimate the surrogates, and append them to the experimental
design.  The edger_method parameter may be used to apply different EdgeR code
paths as outlined in the manual.  If you want to play with non-standard data,
the force argument will round the data and shoe-horn it into EdgeR.
}
\examples{
\dontrun{
 expt <- create_expt(metadata = "metadata.xlsx", gene_info = annotations)
 pretend <- edger_pairwise(expt, model_batch = "sva")
}
}
\seealso{
[edgeR] [deseq_pairwise()] [ebseq_pairwise()] [limma_pairwise()] [basic_pairwise()]
 DOI:10.12688/f1000research.8987.2
}
