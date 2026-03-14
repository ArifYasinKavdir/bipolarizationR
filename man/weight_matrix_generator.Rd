\name{weight_matrix_generator}
\alias{weight_matrix_generator}
\title{Generate a normalised weight matrix for polarization or consensus scoring}
\description{
  Builds an \code{(N+1) x (N+1)} weight matrix where each cell \code{W[i, j]}
  encodes the geometric contribution of a respondent pair whose answers are
  \code{i} and \code{j} on a common integer scale. The matrix is normalised
  so its maximum absolute value equals 1.

  Each cell weight is the product of a distance term \code{d} and an
  agreement term \code{a}, combined via the chosen kernel.
}
\usage{
weight_matrix_generator(
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
  \item{p}{Numeric. Distance exponent (\code{"power"} kernel) or sigma
    bandwidth (\code{"gaussian"} kernel). Default \code{1.0}.}
  \item{q}{Numeric. Agreement exponent. \code{q = 0} disables the agreement
    term; \code{q > 1} suppresses cells far from the diagonal/anti-diagonal.
    Default \code{1.0}.}
  \item{type}{Character. Scoring geometry:
    \describe{
      \item{\code{"polarization"}}{Distance \code{d} = normalised absolute
        difference; agreement \code{a} = proximity to the scale midpoint
        (highest when both answers are near the centre).}
      \item{\code{"consensus"}}{Distance \code{d} = how far the pair sum
        lies above the scale midpoint; agreement \code{a} = inverse of
        absolute difference (highest when answers are equal).}
    }
    Default \code{"polarization"}.}
  \item{kernel}{Character. How \code{d} maps to a weight:
    \describe{
      \item{\code{"power"}}{\code{M = (d^p) * (a^q)}. Smooth and monotone.}
      \item{\code{"gaussian"}}{\code{M = (1 - exp(-d^2 / 2p^2)) * sign(d) * (a^q)}.
        S-shaped; weight saturates at large distances.}
    }
    Default \code{"power"}.}
}
\value{
  A numeric matrix of shape \code{(N+1) x (N+1)} where
  \code{N = end_value - start_value}. All values lie in \code{[-1, 1]}
  and the maximum absolute value equals 1.
}
\examples{
# Default polarization weight matrix (0-10 scale)
W <- weight_matrix_generator()
max(abs(W))   # 1

# Steep distance curve: only extreme disagreements count
W2 <- weight_matrix_generator(p = 3, q = 1, type = "polarization")

# Gaussian kernel with narrow bandwidth
W3 <- weight_matrix_generator(p = 0.3, q = 1,
                               type = "polarization", kernel = "gaussian")

# Consensus geometry on a 1-5 scale
W4 <- weight_matrix_generator(start_value = 1, end_value = 5,
                               type = "consensus")
}
\seealso{
  \code{\link{weight_matrix_visualization}},
  \code{\link{calculate_scores}}
}
