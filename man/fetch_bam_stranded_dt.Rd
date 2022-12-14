% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_fetch.R
\name{fetch_bam_stranded_dt}
\alias{fetch_bam_stranded_dt}
\title{fetch_bam_stranded_dt}
\usage{
fetch_bam_stranded_dt(
  qdt,
  qgr,
  qwin = 50,
  qmet = "summary",
  agg_FUN = sum,
  bfc = NULL,
  n_cores = getOption("mc.cores", 1),
  rname = NULL,
  force_overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{qdt}{data.table containing, file, tall_var, wide_var, and optionally
norm_factor}

\item{qgr}{GRanges of regions to fetch}

\item{qwin}{number of datapoints to use to view each region}

\item{qmet}{strategy to use in each region, summary or sample}

\item{agg_FUN}{function used to aggregate any duplicate tall_var+wide_var mappings in
qdt. Should accept a numeric vector and return a single value.}

\item{bfc}{BiocFileCache object to use to cache data}

\item{n_cores}{number of cores to use. Defaults to value of mc.cores or 1 if
mc.cores is not set.}

\item{rname}{rname to use with cache.  Default is a digest of arguments used
to fetch data.}

\item{force_overwrite}{if TRUE, any contents of cache are overwritten.}
}
\value{
a tidy data.table of profile data.
}
\description{
fetches a tidy data.table of scores from bigwigs, automatically uses cache.
In contrast to fetch_bam_dt, fragment size is not calculated and fragments
are not extended and pileup per strand is stored.
}
\examples{
data("query_gr")
bam_files = dir(system.file('extdata', package = "chiptsne"),
    pattern = ".bam$", full.names = TRUE)
cfg_dt = data.table(file = bam_files)
cfg_dt[, c("tall_var", "wide_var") := tstrsplit(basename(file), "_", keep = 1:2)]
cfg_dt = cfg_dt[tall_var \%in\% c("ESH1", "HUES48", "HUES64")]
cfg_dt[, norm_factor := ifelse(wide_var == "H3K4me3", .3, 1)]
profile_dt = fetch_bam_stranded_dt(cfg_dt, query_gr)
profile_dt
}
