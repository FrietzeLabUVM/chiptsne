% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_velocity.R
\name{prep_velocity}
\alias{prep_velocity}
\title{prep_velocity}
\usage{
prep_velocity(
  tsne_dt,
  tall_var_a,
  tall_var_b,
  id_to_plot = NULL,
  max_plotted = 500,
  delta.min = 0,
  delta.max = Inf,
  angle.min = 0,
  angle.max = 360,
  drop_backgroud = FALSE
)
}
\arguments{
\item{tsne_dt}{tidy data.table of tsne results.  variables names must
include: tx, ty, id, and tall_var}

\item{tall_var_a}{an item in tsne_dt$tall_var. the origin for calculating
angle/distance.}

\item{tall_var_b}{an item in tsne_dt$tall_var. the destination for calculating
angle/distance.}

\item{id_to_plot}{character.  ids in tsne_dt$id. Default of NULL causes all
ids to be used.}

\item{max_plotted}{the maximum number of ids to plot.}

\item{delta.min}{numeric. ids below this delta marked as background.}

\item{delta.max}{numeric. ids above this delta marked as background.}

\item{angle.min}{numeric. ids with angle below this marked as background.}

\item{angle.max}{numeric. ids with angle above this marked as background.}

\item{drop_backgroud}{if TRUE, all background velocities are dropped. Default
is FALSE.}
}
\value{
a tidy data.table of velocity measurements from tall_var_a to tall_var_b
}
\description{
prep_velocity
}
\examples{
data("tsne_dt")
vel_dt = prep_velocity(tsne_dt,
    unique(tsne_dt$tall_var)[2], unique(tsne_dt$tall_var)[3])
vel_dt
}
