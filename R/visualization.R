# --------------------------------------------------------------------------- #
# Weight matrix heatmap                                                        #
# --------------------------------------------------------------------------- #

#' Display a heatmap of the weight matrix.
#'
#' Useful for understanding the geometric effect of \code{p}, \code{q},
#' \code{type}, and \code{kernel} before running any data through the scorer.
#'
#' @inheritParams weight_matrix_generator
#'
#' @return The \code{ggplot} object (invisibly), and prints the plot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' weight_matrix_visualization(p = 2, q = 1, type = "polarization")
#' weight_matrix_visualization(p = 0.3, q = 1, type = "polarization", kernel = "gaussian")
#' }
weight_matrix_visualization <- function(start_value = 0,
                                        end_value   = 10,
                                        p           = 1.0,
                                        q           = 1.0,
                                        type        = "polarization",
                                        kernel      = "power") {
  W          <- weight_matrix_generator(start_value, end_value, p, q, type, kernel)
  full_scale <- seq(start_value, end_value)
  n          <- length(full_scale)

  # Long format: x_lab = column value (y-axis of crosstab), y_lab = row value (x-axis)
  # expand.grid: first arg varies fastest (fills down each column in column-major)
  wm_df        <- expand.grid(y_lab = full_scale, x_lab = full_scale)
  wm_df$value  <- as.vector(W)   # column-major = matches expand.grid order

  vmax <- max(abs(W))
  if (vmax == 0) vmax <- 1

  p_out <- ggplot(wm_df, aes(x = x_lab, y = y_lab, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", value),
                  colour = abs(value) > 0.5),
              size = 3, show.legend = FALSE) +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#D6604D",
      midpoint = 0,
      limits   = c(-vmax, vmax),
      name     = "Weight"
    ) +
    scale_x_continuous(breaks = full_scale) +
    scale_y_reverse(breaks = full_scale) +
    labs(
      title = sprintf(
        "%s Weight Matrix — kernel=%s (p=%g, q=%g)",
        tools::toTitleCase(type), kernel, p, q
      ),
      x = "Column value",
      y = "Row value"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(size = 9)
    )

  print(p_out)
  invisible(p_out)
}


# --------------------------------------------------------------------------- #
# Full dashboard                                                                #
# --------------------------------------------------------------------------- #

#' Compute scores and display a three-panel dashboard for a variable pair.
#'
#' The dashboard has three panels:
#' \itemize{
#'   \item \strong{A} — Bootstrap confidence intervals (horizontal bar chart).
#'   \item \strong{B} — Bootstrap distributions (overlaid histograms).
#'   \item \strong{C} — Per-cell W×P heatmap with diagonal/anti-diagonal separator.
#' }
#' A text summary is also printed to the console.
#'
#' For \code{score_type = "polarization"} only the per-person scores
#' (\code{<x>_per_person}, \code{<y>_per_person}) are visualized. They are
#' recomputed inside every bootstrap resample, so panels A and B show the
#' bootstrap distribution of the per-person scores themselves. The aggregate
#' triangle scores and \code{overall} are still returned in the result but
#' are not plotted.
#'
#' @param df         Data frame. Survey data; one row per respondent.
#' @param x          Character. Column name of the first variable.
#' @param y          Character. Column name of the second variable.
#' @param score_type Character. \code{"polarization"} or \code{"consensus"}.
#' @param start_value,end_value Integer. Scale bounds. Defaults 0 and 10.
#' @param p,q        Numeric. Kernel and agreement parameters.
#' @param kernel     Character. \code{"power"} or \code{"gaussian"}.
#' @param dropna     Logical. Passed to \code{calculate_scores}.
#' @param B_boot     Integer. Bootstrap replications. Default 2000.
#' @param ci         Numeric. Confidence level. Default 0.95.
#' @param random_state Integer. Seed for reproducibility. Default 42.
#'
#' @return A list with element \code{scores} containing the full
#'   \code{calculate_scores} result (including bootstrap distributions).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dash <- dashboard_pair(df, "welsatsoc", "welsatfri",
#'                        score_type = "polarization", B_boot = 200)
#' dash$scores$point$overall
#' }
dashboard_pair <- function(df,
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
                            random_state = 42L) {
  res <- calculate_scores(
    df             = df,
    first_variable  = x,
    second_variable = y,
    start_value    = start_value,
    end_value      = end_value,
    p              = p,
    q              = q,
    score_type     = score_type,
    kernel         = kernel,
    dropna         = dropna,
    B              = B_boot,
    ci             = ci,
    random_state   = random_state,
    keep_dists     = TRUE
  )

  point <- res$point
  boot  <- res$bootstrap

  # Keys shown in panels A and B: per-person scores only for polarization;
  # aggregate scores stay available in the returned result.
  if (score_type == "polarization") {
    keys <- names(boot)[endsWith(names(boot), "_per_person")]
    if (length(keys) == 0) {  # e.g. formula = "per_triangle" has no per-person keys
      keys <- names(boot)[names(boot) != "overall"]
    }
  } else {
    keys <- names(boot)
  }

  # ── Print text summary ──────────────────────────────────────────────────
  sep <- strrep("=", 54)
  cat("\n", sep, "\n", sep = "")
  cat(sprintf("  %s DASHBOARD  —  %s  x  %s\n",
              toupper(score_type), x, y))
  cat(sep, "\n", sep = "")
  cat(sprintf("  kernel=%s  |  B=%d  |  CI=%d%%\n",
              kernel, B_boot, round(ci * 100)))
  cat(sprintf("  Diagnostics: sparsity=%.2f%%\n", point$.sparsity * 100))
  cat(strrep("-", 54), "\n", sep = "")
  cat("  Bootstrap results:\n")
  for (k in keys) {
    ci_k <- boot[[k]]$ci
    cat(sprintf("    %10s: %.4f  [%.4f, %.4f]  SE=%.4f\n",
                k, point[[k]], ci_k[1], ci_k[2], boot[[k]]$se))
  }
  cat(sep, "\n\n", sep = "")

  # ── Shared data ──────────────────────────────────────────────────────────
  WP         <- point$.wp_matrix
  full_scale <- point$.x_labels
  n          <- length(full_scale)
  vmax       <- max(abs(WP))
  if (vmax == 0) vmax <- 1

  point_vals <- vapply(keys, function(k) point[[k]], numeric(1))
  ci_lo      <- vapply(keys, function(k) boot[[k]]$ci[1], numeric(1))
  ci_hi      <- vapply(keys, function(k) boot[[k]]$ci[2], numeric(1))

  # ── Panel A: Bootstrap CIs (horizontal bar chart) ────────────────────────
  bar_colours <- ifelse(point_vals >= 0, "#4878CF", "#D65F5F")

  df_bars <- data.frame(
    key   = factor(keys, levels = rev(keys)),
    value = point_vals,
    lo    = ci_lo,
    hi    = ci_hi,
    fill  = bar_colours,
    stringsAsFactors = FALSE
  )

  p_a <- ggplot(df_bars, aes(x = value, y = key)) +
    geom_col(aes(fill = fill), alpha = 0.75, width = 0.5, show.legend = FALSE) +
    scale_fill_identity() +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.9) +
    geom_text(
      aes(
        label = sprintf("%.4f", value),
        hjust = ifelse(value >= 0, -0.15, 1.15)
      ),
      size = 3
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
    labs(
      title = sprintf("Bootstrap Results\n(B=%d, %d%% CI)", B_boot, round(ci * 100)),
      x     = "Score",
      y     = NULL
    ) +
    theme_minimal(base_size = 10)

  # ── Panel B: Bootstrap distributions (histograms) ────────────────────────
  palette <- c("#4878CF", "#D65F5F", "#6ACC65", "#B47CC7", "#C4AD66")

  boot_list <- lapply(seq_along(keys), function(i) {
    k <- keys[i]
    data.frame(key = k, value = boot[[k]]$dist, stringsAsFactors = FALSE)
  })
  boot_long <- do.call(rbind, boot_list)
  boot_long$key <- factor(boot_long$key, levels = keys)

  vlines_df <- data.frame(
    key   = factor(keys, levels = keys),
    value = point_vals,
    stringsAsFactors = FALSE
  )

  key_colours <- setNames(
    palette[(seq_along(keys) - 1L) %% length(palette) + 1L],
    keys
  )

  p_b <- ggplot(boot_long, aes(x = value, fill = key, colour = key)) +
    geom_histogram(bins = 35, alpha = 0.55, position = "identity",
                   linewidth = 0.2) +
    geom_vline(data = vlines_df,
               aes(xintercept = value, colour = key),
               linetype = "dashed", linewidth = 1.2, show.legend = FALSE) +
    scale_fill_manual(values = key_colours, name = NULL) +
    scale_colour_manual(values = key_colours, name = NULL) +
    labs(
      title = "Bootstrap Distributions\n(all keys)",
      x     = "Score (bootstrap)",
      y     = "Count"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.text = element_text(size = 7))

  # ── Panel C: W×P heatmap ─────────────────────────────────────────────────
  # WP[r, c]: r = x-variable index (row), c = y-variable index (column)
  # In the plot: x-axis = y-variable values, y-axis = x-variable values
  # expand.grid: x_lab (x-axis) varies fastest; y_lab (y-axis) varies slowest
  # as.vector(t(WP)) gives row-major order matching expand.grid(x_lab, y_lab)
  wp_df        <- expand.grid(x_lab = full_scale, y_lab = full_scale)
  wp_df$value  <- as.vector(t(WP))   # transpose for row-major / expand.grid order

  thresh <- vmax * 0.4

  p_c <- ggplot(wp_df, aes(x = x_lab, y = y_lab, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    geom_text(
      data = subset(wp_df, value != 0),
      aes(label  = sprintf("%.3f", value),
          colour = abs(value) > thresh),
      size = 2.2, show.legend = FALSE
    ) +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#D6604D",
      midpoint = 0,
      limits   = c(-vmax, vmax),
      name     = "Score"
    ) +
    scale_x_continuous(breaks = full_scale) +
    scale_y_reverse(breaks = full_scale) +
    labs(
      title = sprintf("Per-cell %s Scores  (W \u00d7 P matrix)",
                      tools::toTitleCase(score_type)),
      x = y,
      y = x
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_text(size = 9)
    )

  # Add diagonal separator
  if (score_type == "polarization") {
    # Main diagonal: x_val == y_val  (where x_lab = y_var values, y_lab = x_var values)
    p_c <- p_c + geom_abline(slope = 1, intercept = 0,
                              colour = "black", linewidth = 1.5,
                              linetype = "dashed", alpha = 0.6)
  } else {
    # Anti-diagonal: x_val + y_val = end_value + start_value
    intercept_val <- end_value + start_value
    p_c <- p_c + geom_abline(slope = -1, intercept = intercept_val,
                              colour = "black", linewidth = 1.5,
                              linetype = "dashed", alpha = 0.6)
  }

  # ── Compose with patchwork ───────────────────────────────────────────────
  dashboard <- (p_a | p_b) / p_c +
    plot_layout(heights = c(1, 2)) +
    plot_annotation(
      title = sprintf(
        "%s Dashboard  [%s]  —  %s  \u00d7  %s",
        tools::toTitleCase(score_type), kernel, x, y
      ),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )

  print(dashboard)
  invisible(list(scores = res))
}
