#' Generate a normalised weight matrix for polarization or consensus scoring.
#'
#' Each cell weight W[i, j] is the product of a distance term \code{d} and an
#' agreement term \code{a}. The matrix is always normalised so its maximum
#' absolute value equals 1.
#'
#' @param start_value Integer. Lower bound of the response scale. Default 0.
#' @param end_value   Integer. Upper bound of the response scale. Default 10.
#' @param p           Numeric. Distance exponent (power kernel) or sigma
#'   bandwidth (gaussian kernel). Default 1.0.
#' @param q           Numeric. Agreement exponent — controls steepness of the
#'   agreement term. Default 1.0.
#' @param type        Character. \code{"polarization"} or \code{"consensus"}.
#'   \itemize{
#'     \item \strong{polarization}: \code{d} = normalised absolute distance;
#'       \code{a} = how far the pair is from both extremes (highest at midpoints).
#'     \item \strong{consensus}: \code{d} = how far above the scale midpoint the
#'       pair sum lies; \code{a} = normalised agreement (inverse distance).
#'   }
#' @param kernel      Character. \code{"power"} or \code{"gaussian"}.
#'   \itemize{
#'     \item \strong{power}:    \code{M = (d^p) * (a^q)}
#'     \item \strong{gaussian}: \code{M = (1 - exp(-d^2 / 2p^2)) * sign(d) * (a^q)}
#'   }
#'
#' @return A numeric matrix of shape \code{(N+1) x (N+1)} where
#'   \code{N = end_value - start_value}, with values in \code{[-1, 1]} and
#'   maximum absolute value equal to 1.
#'
#' @export
#'
#' @examples
#' W <- weight_matrix_generator(0, 10, p = 1, q = 1, type = "polarization")
#' max(abs(W))  # 1
weight_matrix_generator <- function(start_value = 0,
                                    end_value   = 10,
                                    p           = 1.0,
                                    q           = 1.0,
                                    type        = "polarization",
                                    kernel      = "power") {
  N    <- end_value - start_value
  vals <- seq(start_value, end_value)
  n    <- length(vals)

  # Replicate numpy meshgrid(i, j):
  #   I_mat[r, c] = vals[c]   (column values, equivalent to Python's I)
  #   J_mat[r, c] = vals[r]   (row    values, equivalent to Python's J)
  I_mat <- outer(rep(1L, n), vals)
  J_mat <- outer(vals, rep(1L, n))

  if (type == "polarization") {
    d <- abs(I_mat - J_mat) / N
    a <- 1 - abs((I_mat + J_mat) - N) / N
  } else if (type == "consensus") {
    d <- (I_mat + J_mat) / N - 1
    a <- (N - abs(I_mat - J_mat)) / N
  } else {
    stop("type must be 'polarization' or 'consensus'")
  }

  if (kernel == "power") {
    M <- (d ^ p) * (a ^ q)
  } else if (kernel == "gaussian") {
    g <- (1 - exp(-(d ^ 2) / (2 * p ^ 2))) * sign(d)
    M <- g * (a ^ q)
  } else {
    stop("kernel must be 'power' or 'gaussian'")
  }

  max_val <- max(abs(M))
  if (max_val > 0) M <- M / max_val
  M
}
