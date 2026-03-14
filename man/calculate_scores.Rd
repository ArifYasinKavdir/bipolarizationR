\name{calculate_scores}
\alias{calculate_scores}
\title{Calculate polarization or consensus scores for a pair of survey variables}
\description{
  Builds a joint frequency table for (\code{first_variable},
  \code{second_variable}), multiplies it element-wise by a geometrically
  derived weight matrix, and summarises the result as a scalar score per
  triangle (polarization) or per half of the anti-diagonal (consensus).
  Non-parametric bootstrap resampling provides standard errors and
  confidence intervals.
}
\usage{
calculate_scores(
  df,
  first_variable,
  second_variable,
  start_value  = 0L,
  end_value    = 10L,
  p            = 1.0,
  q            = 1.0,
  score_type   = "polarization",
  kernel       = "power",
  dropna       = FALSE,
  B            = 2000L,
  ci           = 0.95,
  random_state = 42L,
  keep_dists   = FALSE
)
}
\arguments{
  \item{df}{A data frame. Survey data with one row per respondent.}
  \item{first_variable}{Character. Column name of the first variable —
    becomes the row axis of the crosstab.}
  \item{second_variable}{Character. Column name of the second variable —
    becomes the column axis.}
  \item{start_value}{Integer. Minimum value on the response scale. Default \code{0}.}
  \item{end_value}{Integer. Maximum value on the response scale. Default \code{10}.}
  \item{p}{Numeric. Distance kernel parameter. For \code{kernel = "power"},
    this is the distance exponent; for \code{kernel = "gaussian"}, this is
    the sigma bandwidth. Default \code{1.0}.}
  \item{q}{Numeric. Agreement exponent — controls steepness of the agreement
    term. \code{q = 0} disables it; higher values suppress extreme cells.
    Default \code{1.0}.}
  \item{score_type}{Character. \code{"polarization"} (opposing camps) or
    \code{"consensus"} (shared conviction). Default \code{"polarization"}.}
  \item{kernel}{Character. \code{"power"} (monotone, \code{M = d^p * a^q})
    or \code{"gaussian"} (S-shaped, saturating). Default \code{"power"}.}
  \item{dropna}{Logical. Informational; both columns are always inner-aligned
    (NAs removed from either). Default \code{FALSE}.}
  \item{B}{Integer. Number of bootstrap replications. Default \code{2000}.}
  \item{ci}{Numeric. Confidence level, e.g. \code{0.95} for a 95\% interval.
    Default \code{0.95}.}
  \item{random_state}{Integer or \code{NULL}. Random seed for reproducible
    bootstrap samples. Pass \code{NULL} for a random seed. Default \code{42}.}
  \item{keep_dists}{Logical. If \code{TRUE}, include the raw bootstrap
    distribution vectors (\code{dist}) in the returned bootstrap list.
    Default \code{FALSE}.}
}
\value{
  A named list with four elements:

  \describe{
    \item{\code{type}}{Character. Echo of \code{score_type}.}

    \item{\code{point}}{Named list of point estimates. For
      \code{score_type = "polarization"}: \code{first_variable} (lower
      triangle score), \code{second_variable} (upper triangle score), and
      \code{overall} (sum of both). For \code{score_type = "consensus"}:
      \code{negatives}, \code{positives}, and \code{overall}. All types
      also include diagnostic elements prefixed with \code{.}:
      \code{.sparsity}, \code{.wp_matrix}, \code{.x_labels},
      \code{.n_above}/\code{.n_below} (polarization) or
      \code{.n_anti_above}/\code{.n_anti_below} (consensus).}

    \item{\code{bootstrap}}{Named list with one entry per public score key.
      Each entry contains \code{mean}, \code{se}, \code{ci} (length-2
      numeric vector \code{c(lower, upper)}), and optionally \code{dist}
      (numeric vector of length \code{B}) if \code{keep_dists = TRUE}.}

    \item{\code{meta}}{List with \code{B}, \code{ci}, and \code{kernel}.}
  }
}
\examples{
\dontrun{
# Basic polarization score
result <- calculate_scores(df, "idemus", "idekemalist",
                           score_type = "polarization", B = 500)
result$point$overall
result$bootstrap$overall$ci
result$bootstrap$overall$se

# Consensus with gaussian kernel
result2 <- calculate_scores(df, "ideleft", "ideright",
                            score_type = "consensus",
                            kernel = "gaussian", p = 0.5, B = 2000)

# Keep bootstrap distributions for custom plots
result3 <- calculate_scores(df, "idemus", "idekemalist",
                            keep_dists = TRUE)
hist(result3$bootstrap$overall$dist)
}
}
\seealso{
  \code{\link{weight_matrix_generator}},
  \code{\link{dashboard_pair}},
  \code{\link{weight_matrix_visualization}}
}
