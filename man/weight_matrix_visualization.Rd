\name{weight_matrix_visualization}
\alias{weight_matrix_visualization}
\title{Display a heatmap of the weight matrix}
\description{
  Generates and prints an annotated \code{ggplot2} heatmap of the weight
  matrix produced by \code{\link{weight_matrix_generator}}. Each cell
  displays its numeric value; cells with absolute value above 0.5 use white
  text for contrast. Useful for inspecting the geometric effect of
  \code{p}, \code{q}, \code{type}, and \code{kernel} before running any
  data through the scorer.
}
\usage{
weight_matrix_visualization(
  start_value = 0,
  end_value   = 10,
  p           = 1.0,
  q           = 1.0,
  type        = "polarization",
  kernel      = "power"
)
}
\arguments{
  \item{start_value}{Integer. Lower bound of the response scale. Default \code{0}.}
  \item{end_value}{Integer. Upper bound of the response scale. Default \code{10}.}
  \item{p}{Numeric. Distance exponent or sigma bandwidth. Default \code{1.0}.}
  \item{q}{Numeric. Agreement exponent. Default \code{1.0}.}
  \item{type}{Character. \code{"polarization"} or \code{"consensus"}.
    Default \code{"polarization"}.}
  \item{kernel}{Character. \code{"power"} or \code{"gaussian"}.
    Default \code{"power"}.}
}
\value{
  The \code{ggplot} object (returned invisibly). The plot is printed as a
  side effect.
}
\examples{
# Default polarization weight matrix
weight_matrix_visualization()

# Steep distance exponent: only extreme disagreements count
weight_matrix_visualization(p = 3, q = 1, type = "polarization")

# Gaussian kernel comparison
weight_matrix_visualization(p = 0.3, q = 1,
                            type = "polarization", kernel = "gaussian")

# Consensus geometry on a 1-5 Likert scale
weight_matrix_visualization(start_value = 1, end_value = 5,
                            type = "consensus")

# High q: suppresses cells far from midpoint
weight_matrix_visualization(p = 1, q = 3, type = "polarization")
}
\seealso{
  \code{\link{weight_matrix_generator}},
  \code{\link{dashboard_pair}}
}
