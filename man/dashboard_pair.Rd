\name{dashboard_pair}
\alias{dashboard_pair}
\title{Compute scores and display a three-panel dashboard for a variable pair}
\description{
  Runs \code{\link{calculate_scores}} and renders a three-panel
  \code{ggplot2} figure composed with \code{patchwork}. A text summary is
  also printed to the console.

  \describe{
    \item{Panel A (top-left)}{Horizontal bar chart of point estimates with
      bootstrap CI error bars.}
    \item{Panel B (top-right)}{Overlaid histograms of the bootstrap
      distributions for each score key, with dashed vertical lines at the
      point estimates.}
    \item{Panel C (bottom, full-width)}{Heatmap of the element-wise
      \code{W x P} matrix. The main diagonal (polarization) or
      anti-diagonal (consensus) is drawn as a dashed separator. Cell values
      are printed inside each tile.}
  }

  For \code{score_type = "polarization"} only the per-person scores
  (\code{<x>_per_person}, \code{<y>_per_person}) are visualized. They are
  recomputed inside every bootstrap resample, so panels A and B show the
  bootstrap distribution of the per-person scores themselves. The aggregate
  triangle scores and \code{overall} are still returned in the result but
  are not plotted.
}
\usage{
dashboard_pair(
  df,
  x,
  y,
  score_type   = "polarization",
  start_value  = 0L,
  end_value    = 10L,
  p            = 1.0,
  q            = 1.0,
  kernel       = "power",
  dropna       = FALSE,
  B_boot       = 2000L,
  ci           = 0.95,
  random_state = 42L
)
}
\arguments{
  \item{df}{A data frame. Survey data with one row per respondent.}
  \item{x}{Character. Column name of the first variable.}
  \item{y}{Character. Column name of the second variable.}
  \item{score_type}{Character. \code{"polarization"} or \code{"consensus"}.
    Default \code{"polarization"}.}
  \item{start_value}{Integer. Minimum value on the response scale. Default \code{0}.}
  \item{end_value}{Integer. Maximum value on the response scale. Default \code{10}.}
  \item{p}{Numeric. Distance kernel parameter. Default \code{1.0}.}
  \item{q}{Numeric. Agreement exponent. Default \code{1.0}.}
  \item{kernel}{Character. \code{"power"} or \code{"gaussian"}. Default \code{"power"}.}
  \item{dropna}{Logical. Passed to \code{\link{calculate_scores}}. Default \code{FALSE}.}
  \item{B_boot}{Integer. Number of bootstrap replications. Default \code{2000}.}
  \item{ci}{Numeric. Confidence level. Default \code{0.95}.}
  \item{random_state}{Integer. Random seed. Default \code{42}.}
}
\value{
  A named list with one element:
  \describe{
    \item{\code{scores}}{The full result of \code{\link{calculate_scores}}
      called with \code{keep_dists = TRUE}, including bootstrap distribution
      vectors.}
  }

  The three-panel figure is printed as a side effect.
}
\examples{
\dontrun{
# Polarization dashboard
dash <- dashboard_pair(df, "idemus", "idekemalist",
                       score_type = "polarization", B_boot = 2000)
dash$scores$point$overall

# Quick exploration with fewer bootstrap samples
dash2 <- dashboard_pair(df, "welsatsoc", "welsatfri",
                        score_type = "polarization", B_boot = 200)

# Consensus mode
dash3 <- dashboard_pair(df, "ideleft", "ideright",
                        score_type = "consensus", B_boot = 500)
}
}
\seealso{
  \code{\link{calculate_scores}},
  \code{\link{weight_matrix_visualization}}
}
