% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_fetch.R
\name{stsFetchTsneInput}
\alias{stsFetchTsneInput}
\title{stsFetchTsneInput}
\usage{
stsFetchTsneInput(
  qdt,
  qgr,
  region_size = NULL,
  file_format = NULL,
  qwin = 50,
  qmet = "summary",
  cap_value = 20,
  fetch_FUN = NULL,
  agg_FUN = NULL,
  high_on_right = TRUE,
  bfc = BiocFileCache::BiocFileCache(),
  n_cores = getOption("mc.cores", 1),
  rname = digest::digest(list(qgr, qdt, qwin, qmet, fetch_FUN)),
  force_overwrite = FALSE,
  verbose = TRUE,
  skip_checks = FALSE
)
}
\arguments{
\item{qdt}{data.table containing, file, tall_var, wide_var, and optionally
norm_factor}

\item{qgr}{GRanges of regions to fetch}

\item{region_size}{numeric, if provided will be used to set all widths on
qgr.}

\item{file_format}{character, describing file format. must be one of bam or
bw.  Will be automatically determined from file extension if NULL.}

\item{qwin}{number of datapoints to use to view each region}

\item{qmet}{strategy to use in each region, summary or sample}

\item{cap_value}{maximum allowed value.  useful for outlier handling.}

\item{fetch_FUN}{function used to load and process data from files.}

\item{agg_FUN}{function used to aggregate any duplicate tall_var+wide_var mappings in
qdt.  Should accept a numeric vector and return a single value.}

\item{high_on_right}{if TRUE, profiles where highest signal is on the left
are flipped.}

\item{bfc}{BiocFileCache object to use to cache data}

\item{n_cores}{number of cores to use. Defaults to value of mc.cores or 1 if
mc.cores is not set.}

\item{rname}{rname to use with cache.  Default is a digest of arguments used
to fetch data.}

\item{force_overwrite}{if TRUE, any contents of cache are overwritten. Default is FALSE.}

\item{verbose}{if TRUE outputs progress and log messages}

\item{skip_checks}{if TRUE, tall_var and wide_var completeness checks are skipped.
tall_var/wide_var combinations must be complete to run StsRunTsne(). Default is FALSE.}
}
\value{
a tidy data.table of profile data.
}
\description{
retrieve data from series of either bigwig or bam files to a tidy data.table
for later input to stsRunTsne
}
\details{
validates inputs prepares query_gr fetches either bam or bigwig based on
extension of first file transform and normalizes profiles after fetch
}
\examples{
data("query_gr")
bw_files = dir(system.file('extdata', package = "chiptsne"), pattern = ".bw$", full.names = TRUE)
cfg_dt = data.table(file = bw_files)
cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
cfg_dt = cfg_dt[tall_var \%in\% c("ESH1", "HUES48", "HUES64")]
cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
profile_dt = stsFetchTsneInput(cfg_dt, query_gr)
profile_dt
}
