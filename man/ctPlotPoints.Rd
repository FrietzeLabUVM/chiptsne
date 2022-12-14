% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ctPlots.R
\name{ctPlotPoints}
\alias{ctPlotPoints}
\title{ctPlotPoints}
\usage{
ctPlotPoints(
  sts,
  feature_name = NULL,
  signal_name = NULL,
  xmin = -Inf,
  xmax = Inf,
  profile_value = ssvQC:::val2var[sts$signal_config$plot_value],
  profile_value_label = ssvQC:::val2lab[sts$signal_config$plot_value],
  bg_color = "gray60",
  agg_FUN = max,
  point_size = 1
)
}
\arguments{
\item{sts}{A ChIPtSNE object with ssvQC.prepSignal already called.}

\item{feature_name}{Feature name present in sts.  With default of NULL, first
feature name will be used.}

\item{signal_name}{Signal name present in sts.  With default of NULL, first
signal name will be used.}

\item{xmin}{The min range of x-values to apply agg_FUN to.}

\item{xmax}{The max range of x-values to apply agg_FUN to.}

\item{profile_value}{Value to use for profiles. Default is the plot_value
defined in the signal config of sts.}

\item{profile_value_label}{Label to use for profile scale. Default is the
plot_value label defined in the signal config of sts.}

\item{bg_color}{Color to use for plot background. Default is "gray60".}

\item{agg_FUN}{A function appled to all y_ values in xmin to xmax range. Must
accept single numeric vector and return 1 value.}

\item{point_size}{Size to plot points at. Default is 1.}
}
\value{
ggplot where each profile is summarized by a single point in t-SNE space.
}
\description{
ctPlotPoints
}
\examples{
data(ex_sts)
ctPlotPoints(sts)

}
