# --------------------------------------------------------------------------- #
# Internal helpers                                                             #
# --------------------------------------------------------------------------- #

# Normalise each mask region of A independently by its own sum.
.normalize_by_masks <- function(A, masks) {
  out <- matrix(0.0, nrow = nrow(A), ncol = ncol(A))
  for (m in masks) {
    s <- sum(A[m])
    if (s > 0) out[m] <- A[m] / s
  }
  out
}


# Normalise all cells of A by the grand total so the whole matrix sums to 1.
.normalize_globally <- function(A) {
  total <- sum(A)
  if (total > 0) A / total else matrix(0.0, nrow = nrow(A), ncol = ncol(A))
}


# Compute a single point-estimate score for the pair (x, y).
#
# Returns a named list whose public keys depend on score_type:
#   polarization : x, y, overall  (+ diagnostics prefixed with ".")
#   consensus    : negatives, positives, overall  (+ diagnostics)
.single_run_score <- function(df,
                               x,
                               y,
                               start_value,
                               end_value,
                               p,
                               q,
                               score_type,
                               dropna,
                               kernel  = "power",
                               formula = "global") {
  full_scale <- seq(start_value, end_value)
  n          <- length(full_scale)

  # Extract columns and inner-align (drop NAs in either)
  x_vals <- suppressWarnings(as.integer(df[[x]]))
  y_vals <- suppressWarnings(as.integer(df[[y]]))
  valid  <- !is.na(x_vals) & !is.na(y_vals)
  x_vals <- x_vals[valid]
  y_vals <- y_vals[valid]

  # Build joint count matrix, forcing all scale values to appear
  ct <- table(
    factor(x_vals, levels = full_scale),
    factor(y_vals, levels = full_scale)
  )
  A <- matrix(as.double(ct), nrow = n, ncol = n)

  W <- weight_matrix_generator(
    start_value = start_value, end_value = end_value,
    p = p, q = q, type = score_type, kernel = kernel
  )

  if (!identical(dim(W), dim(A))) {
    stop(sprintf(
      "weight matrix shape (%d x %d) must match count table shape (%d x %d)",
      nrow(W), ncol(W), nrow(A), ncol(A)
    ))
  }

  # Index matrices (1-based), equivalent to numpy np.indices(A.shape) 0-based
  r_mat <- outer(seq_len(n), rep(1L, n))   # r_mat[r,c] = r
  c_mat <- outer(rep(1L, n), seq_len(n))   # c_mat[r,c] = c

  if (score_type == "polarization") {
    masks <- list(
      above = r_mat < c_mat,   # x < y  (upper triangle)
      below = r_mat > c_mat,   # x > y  (lower triangle)
      diag  = r_mat == c_mat
    )
    n_above <- sum(A[masks$above])
    n_below <- sum(A[masks$below])

    if (formula == "global") {
      total       <- sum(A)
      P           <- .normalize_globally(A)
      WP          <- P * W
      score_above <- sum(WP[masks$above])
      score_below <- sum(WP[masks$below])
      pct_above   <- if (total > 0) n_above / total * 100 else 0
      pct_below   <- if (total > 0) n_below / total * 100 else 0

      result <- list(
        type       = "polarization",
        overall    = sum(WP),
        .sparsity  = mean(A == 0),
        .wp_matrix = WP,
        .x_labels  = full_scale,
        .n_above   = n_above,
        .n_below   = n_below
      )
      result[[x]]                       <- score_below
      result[[y]]                       <- score_above
      result[[paste0(x, "_per_person")]] <- if (pct_below > 0) score_below / pct_below else 0
      result[[paste0(y, "_per_person")]] <- if (pct_above > 0) score_above / pct_above else 0
      return(result)

    } else {  # "per_triangle"
      P  <- .normalize_by_masks(A, masks)
      WP <- P * W

      result <- list(
        type       = "polarization",
        overall    = sum(WP),
        .sparsity  = mean(A == 0),
        .wp_matrix = WP,
        .x_labels  = full_scale,
        .n_above   = n_above,
        .n_below   = n_below
      )
      result[[x]] <- sum(WP[masks$below])
      result[[y]] <- sum(WP[masks$above])
      return(result)
    }

  } else if (score_type == "consensus") {
    # s_mat = r_mat + c_mat  (1-based sums, range 2 to 2n)
    # Python (0-based) anti-diag: I + J == n - 1
    # Converting: (r-1) + (c-1) = n - 1  =>  r + c = n + 1
    s_mat <- r_mat + c_mat
    masks <- list(
      anti_above = s_mat <= n,         # Python: s < n-1  (both low)
      anti_below = s_mat >= (n + 2L),  # Python: s > n-1  (both high)
      anti_diag  = s_mat == (n + 1L)   # Python: s == n-1
    )
    P  <- .normalize_by_masks(A, masks)
    WP <- P * W

    list(
      type          = "consensus",
      negatives     = sum(WP[masks$anti_above]),
      positives     = sum(WP[masks$anti_below]),
      overall       = sum(WP),
      .sparsity     = mean(A == 0),
      .wp_matrix    = WP,
      .x_labels     = full_scale,
      .n_anti_above = sum(A[masks$anti_above]),
      .n_anti_below = sum(A[masks$anti_below])
    )

  } else {
    stop("score_type must be 'polarization' or 'consensus'")
  }
}


# Run a non-parametric bootstrap around the point estimate.
# Returns list(point = ..., summary = ...) where summary has per-key
# mean, se, ci, dist (always present internally; stripped later if needed).
.bootstrap_scores <- function(df,
                               x,
                               y,
                               start_value,
                               end_value,
                               p,
                               q,
                               score_type,
                               dropna,
                               B,
                               ci,
                               seed,
                               kernel  = "power",
                               formula = "global") {
  if (!is.null(seed)) set.seed(seed)

  point       <- .single_run_score(df, x, y, start_value, end_value, p, q,
                                    score_type, dropna, kernel, formula)
  public_keys <- names(point)[names(point) != "type" &
                                !startsWith(names(point), ".")]

  nr   <- nrow(df)
  dist <- setNames(
    lapply(public_keys, function(k) numeric(B)),
    public_keys
  )

  for (b in seq_len(B)) {
    idx <- sample.int(nr, nr, replace = TRUE)
    s   <- .single_run_score(df[idx, , drop = FALSE], x, y,
                              start_value, end_value, p, q,
                              score_type, dropna, kernel, formula)
    for (k in public_keys) dist[[k]][b] <- s[[k]]
  }

  alpha  <- (1 - ci) / 2
  probs  <- c(alpha, 1 - alpha)

  summary <- setNames(
    lapply(public_keys, function(k) {
      a <- dist[[k]]
      list(
        mean = mean(a),
        se   = sd(a),
        ci   = as.numeric(quantile(a, probs = probs, names = FALSE)),
        dist = a
      )
    }),
    public_keys
  )

  list(point = point, summary = summary)
}


# --------------------------------------------------------------------------- #
# Public API                                                                   #
# --------------------------------------------------------------------------- #

#' Calculate polarization or consensus scores for a pair of survey variables.
#'
#' Builds a joint frequency table for (\code{first_variable},
#' \code{second_variable}), multiplies it by a weighted kernel matrix, and
#' summarises the result as a scalar score per triangle (polarization) or per
#' half of the anti-diagonal (consensus). Bootstrap resampling provides
#' standard errors and confidence intervals.
#'
#' @param df             Data frame. Survey data; one row per respondent.
#' @param first_variable  Character. Column name — row axis of the crosstab.
#' @param second_variable Character. Column name — column axis.
#' @param start_value    Integer. Minimum value on the response scale. Default 0.
#' @param end_value      Integer. Maximum value on the response scale. Default 10.
#' @param p              Numeric. Distance exponent (power) or sigma (gaussian). Default 1.0.
#' @param q              Numeric. Agreement exponent. Default 1.0.
#' @param score_type     Character. \code{"polarization"} or \code{"consensus"}. Default \code{"polarization"}.
#' @param kernel         Character. \code{"power"} or \code{"gaussian"}. Default \code{"power"}.
#' @param formula        Character. \code{"global"} (default) or \code{"per_triangle"}.
#'   \code{"global"} divides all cells by the grand total so the whole matrix
#'   sums to 1; polarization scores reflect each triangle's share of the full
#'   sample and per-person keys (\code{<var>_per_person}) are added.
#'   \code{"per_triangle"} is the legacy behaviour: each triangle is normalised
#'   independently.
#' @param dropna         Logical. Informational flag; both columns are always
#'   inner-aligned (NAs dropped). Default \code{FALSE}.
#' @param B              Integer. Number of bootstrap replications. Default 2000.
#' @param ci             Numeric. Confidence level (e.g. 0.95 for 95\%). Default 0.95.
#' @param random_state   Integer or \code{NULL}. Seed for reproducibility. Default 42.
#' @param keep_dists     Logical. If \code{TRUE}, include raw bootstrap arrays
#'   (\code{dist}) in the returned bootstrap list. Default \code{FALSE}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{type}}{Echo of \code{score_type}.}
#'   \item{\code{point}}{Named list of point estimates. For polarization with
#'     \code{formula="global"}: \code{first_variable}, \code{second_variable},
#'     \code{overall}, \code{<first_variable>_per_person},
#'     \code{<second_variable>_per_person}, plus diagnostics.
#'     With \code{formula="per_triangle"}: \code{first_variable},
#'     \code{second_variable}, \code{overall}, plus diagnostics.
#'     For consensus: \code{negatives}, \code{positives}, \code{overall}, plus
#'     diagnostics.}
#'   \item{\code{bootstrap}}{Named list per public score key, each containing
#'     \code{mean}, \code{se}, \code{ci} (length-2 vector), and optionally
#'     \code{dist}.}
#'   \item{\code{meta}}{List with \code{B}, \code{ci}, \code{kernel}, \code{formula}.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- calculate_scores(df, "idemus", "idekemalist",
#'                            score_type = "polarization", B = 500)
#' result$point$overall
#' result$bootstrap$overall$ci
#' }
calculate_scores <- function(df,
                              first_variable,
                              second_variable,
                              start_value  = 0L,
                              end_value    = 10L,
                              p            = 1.0,
                              q            = 1.0,
                              score_type   = "polarization",
                              kernel       = "power",
                              formula      = "global",
                              dropna       = FALSE,
                              B            = 2000L,
                              ci           = 0.95,
                              random_state = 42L,
                              keep_dists   = FALSE) {
  res <- .bootstrap_scores(
    df          = df,
    x           = first_variable,
    y           = second_variable,
    start_value = start_value,
    end_value   = end_value,
    p           = p,
    q           = q,
    score_type  = score_type,
    dropna      = dropna,
    B           = B,
    ci          = ci,
    seed        = random_state,
    kernel      = kernel,
    formula     = formula
  )

  boot <- res$summary
  if (!keep_dists) {
    for (k in names(boot)) boot[[k]]$dist <- NULL
  }

  list(
    type      = res$point$type,
    point     = res$point,
    bootstrap = boot,
    meta      = list(B = B, ci = ci, kernel = kernel, formula = formula)
  )
}
