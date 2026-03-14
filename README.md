# bipolarizationR

An R package for detecting and quantifying **polarization** and **consensus** in survey data using weighted kernel scoring and non-parametric bootstrap inference.

This library is written based on polarization research conducted by **Zübeyir Nişancı**, **Belkıs Yüce** and **Arif Yasin Kavdır**. It is a faithful R port of the [`bipolarization`](https://github.com/ArifYasinKavdir/bipolarization) Python library.

---

## Installation

### From GitHub (recommended)

```r
# install.packages("devtools")
devtools::install_github("ArifYasinKavdir/bipolarizationR")
```

### Local installation (development)

```r
devtools::install_local("path/to/bipolarizationR")
# or, for interactive development:
devtools::load_all("path/to/bipolarizationR")
```

**Dependencies:** `ggplot2`, `patchwork` (installed automatically).

---

## Quick Start

```r
library(bipolarizationR)

# df is your survey data frame
result <- calculate_scores(df, "idemus", "idekemalist", score_type = "polarization")

result$point$overall          # point estimate
result$bootstrap$overall$ci   # 95% confidence interval
result$bootstrap$overall$se   # bootstrap standard error
```

---

## Concepts

### The Scoring Idea

For a pair of variables `(x, y)` measured on a common integer scale `[start_value, end_value]`:

1. Build an `(N+1) × (N+1)` joint count matrix `A` where `N = end_value − start_value`.
   `A[i, j]` = number of respondents who answered `i` on `x` and `j` on `y`.

2. Normalise each region of `A` independently (upper triangle, lower triangle, diagonal) to get a probability matrix `P`.

3. Build a weight matrix `W` encoding how much each cell `(i, j)` contributes geometrically.

4. The score is `sum(W * P)` over the relevant region.

### Weight Matrix

Each weight `W[i, j]` is the product of two terms:

| Term      | Symbol | Meaning                                          |
| --------- | ------ | ------------------------------------------------ |
| Distance  | `d`    | How far apart `i` and `j` are on the scale       |
| Agreement | `a`    | How much `i` and `j` tend toward the same region |

The matrix is always normalised so its maximum absolute value equals 1.

### Kernels

| Kernel       | Formula                                    | Effect                                                                                       |
| ------------ | ------------------------------------------ | -------------------------------------------------------------------------------------------- |
| `"power"`    | `M = (d^p) * (a^q)`                        | Smooth, monotone. `p < 1` emphasises moderate distances; `p > 1` emphasises extremes.        |
| `"gaussian"` | `M = (1 − exp(−d²/2p²)) * sign(d) * (a^q)` | S-shaped. `p` is the bandwidth sigma; small `p` makes weight rise steeply near the diagonal. |

### Polarization vs Consensus

#### Polarization

The main diagonal (`i == j`) is the "agreement axis". The two off-diagonal triangles represent two groups:

- **Lower triangle** (`i > j`): respondents where `x > y` — the "x-leaning" group.
- **Upper triangle** (`i < j`): respondents where `x < y` — the "y-leaning" group.

Weight is highest for cells that are far from the diagonal (high disagreement) **and** close to the middle of the scale. High score → the two variables split respondents into opposing camps.

#### Consensus

The anti-diagonal (`i + j == N`) is the "midpoint axis". The two halves represent:

- **Negative half** (`i + j < N`): both ratings tend low.
- **Positive half** (`i + j > N`): both ratings tend high.

Weight is highest where the pair sum is far from `N` **and** `i ≈ j` (strong internal consistency). High score → respondents rate both variables similarly (either both high or both low).

### Bootstrap Inference

All scores are accompanied by confidence intervals from non-parametric bootstrap resampling:

1. Compute the point estimate on the full dataset.
2. Draw `B` bootstrap samples (resampling rows with replacement).
3. Compute the score on each sample.
4. Report `mean`, `se` (standard deviation of bootstrap distribution), and the `ci`-level percentile interval.

---

## API Reference

### `calculate_scores()`

Main entry point. Computes point estimates and bootstrap CIs for a variable pair.

```r
calculate_scores(
  df,
  first_variable,
  second_variable,
  start_value  = 0,
  end_value    = 10,
  p            = 1.0,
  q            = 1.0,
  score_type   = "polarization",
  kernel       = "power",
  dropna       = FALSE,
  B            = 2000,
  ci           = 0.95,
  random_state = 42,
  keep_dists   = FALSE
)
```

**Required arguments**

| Argument          | Type       | Description                            |
| ----------------- | ---------- | -------------------------------------- |
| `df`              | data.frame | Survey data; one row per respondent    |
| `first_variable`  | character  | Column name — row axis of the crosstab |
| `second_variable` | character  | Column name — column axis              |

**Keyword arguments**

| Argument       | Default          | Description                                                 |
| -------------- | ---------------- | ----------------------------------------------------------- |
| `start_value`  | `0`              | Minimum value on the response scale                         |
| `end_value`    | `10`             | Maximum value on the response scale                         |
| `p`            | `1.0`            | Distance exponent (`power`) or sigma bandwidth (`gaussian`) |
| `q`            | `1.0`            | Agreement exponent                                          |
| `score_type`   | `"polarization"` | `"polarization"` or `"consensus"`                           |
| `kernel`       | `"power"`        | `"power"` or `"gaussian"`                                   |
| `dropna`       | `FALSE`          | NaN handling flag (both columns are always inner-aligned)   |
| `B`            | `2000`           | Number of bootstrap replications                            |
| `ci`           | `0.95`           | Confidence level (e.g. `0.95` → 95% CI)                     |
| `random_state` | `42`             | Integer seed for reproducibility; pass `NULL` for random    |
| `keep_dists`   | `FALSE`          | If `TRUE`, include raw bootstrap arrays in the result       |

**Returns:** named list — see [Return Value Schema](#return-value-schema).

---

### `weight_matrix_generator()`

Low-level function. Builds and returns the raw weight matrix.

```r
weight_matrix_generator(
  start_value = 0,
  end_value   = 10,
  p           = 1.0,
  q           = 1.0,
  type        = "polarization",
  kernel      = "power"
)
```

**Returns:** numeric matrix of shape `(N+1) × (N+1)`, values in `[−1, 1]`, max absolute value = 1.

---

### `dashboard_pair()`

Runs `calculate_scores()` and displays a three-panel figure. Also prints a text summary to the console.

```r
dashboard_pair(
  df,
  x,
  y,
  score_type   = "polarization",
  start_value  = 0,
  end_value    = 10,
  p            = 1.0,
  q            = 1.0,
  kernel       = "power",
  dropna       = FALSE,
  B_boot       = 2000,
  ci           = 0.95,
  random_state = 42
)
```

| Panel             | Content                                                              |
| ----------------- | -------------------------------------------------------------------- |
| **A** (top-left)  | Horizontal bar chart of point estimates with bootstrap CI error bars |
| **B** (top-right) | Overlaid histograms of the bootstrap distributions                   |
| **C** (bottom)    | Heatmap of the W×P matrix with diagonal/anti-diagonal separator      |

**Returns:** `list(scores = <calculate_scores result>)`

---

### `weight_matrix_visualization()`

Displays an annotated heatmap of the weight matrix. Useful for inspecting the effect of parameters before running data through the scorer.

```r
weight_matrix_visualization(
  start_value = 0,
  end_value   = 10,
  p           = 1.0,
  q           = 1.0,
  type        = "polarization",
  kernel      = "power"
)
```

**Returns:** the `ggplot` object (invisibly).

---

## Return Value Schema

```r
list(
  type = "polarization",    # or "consensus"

  point = list(
    # --- polarization ---
    first_variable  = <float>,  # score from lower triangle (x > y)
    second_variable = <float>,  # score from upper triangle (x < y)
    overall         = <float>,  # total score (sum of both triangles)

    # --- consensus ---
    negatives = <float>,        # score from anti-upper half (both low)
    positives = <float>,        # score from anti-lower half (both high)
    overall   = <float>,

    # --- diagnostics (both types, prefixed with ".") ---
    .sparsity  = <float>,       # fraction of empty cells in the count matrix
    .wp_matrix = <matrix>,      # full W×P matrix (N+1 × N+1)
    .x_labels  = <integer>,     # scale tick labels
    # polarization only:
    .n_above = <float>,         # respondent count in upper triangle
    .n_below = <float>,         # respondent count in lower triangle
    # consensus only:
    .n_anti_above = <float>,
    .n_anti_below = <float>
  ),

  bootstrap = list(
    overall = list(
      mean = <float>,           # bootstrap mean
      se   = <float>,           # bootstrap standard error
      ci   = c(lower, upper),   # percentile interval
      dist = <numeric vector>   # present only if keep_dists = TRUE
    )
    # one entry per public score key
  ),

  meta = list(
    B      = <integer>,
    ci     = <float>,
    kernel = <character>
  )
)
```

> **Sparsity warning:** when `.sparsity` is high (> 0.5), many cells are empty. This can bias scores and widen CIs. Consider collapsing the scale or increasing sample size.

---

## Parameter Tuning Guide

### Choosing `score_type`

| Goal                                              | Recommended      |
| ------------------------------------------------- | ---------------- |
| Detect opposing camps (high on one, low on other) | `"polarization"` |
| Detect shared conviction (both high or both low)  | `"consensus"`    |

### Choosing `kernel`

| Kernel       | When to use                                                                                             |
| ------------ | ------------------------------------------------------------------------------------------------------- |
| `"power"`    | Default; interpretable; monotone distance effect                                                        |
| `"gaussian"` | Sharp threshold: cells close to diagonal get near-zero weight, weight rises steeply beyond distance `p` |

### Choosing `p`

| Value (power kernel) | Effect                                            |
| -------------------- | ------------------------------------------------- |
| `p < 1`              | Concave — moderate disagreements weighted heavily |
| `p = 1`              | Linear — proportional to distance                 |
| `p > 1`              | Convex — only large disagreements count           |

For the gaussian kernel, `p` = sigma: smaller `p` means faster saturation.

### Choosing `q`

| Value   | Effect                                                                         |
| ------- | ------------------------------------------------------------------------------ |
| `q = 0` | Agreement term disabled — pure distance scoring                                |
| `q = 1` | Linear agreement — default, balanced                                           |
| `q > 1` | Only cells near scale midpoint/anti-diagonal are weighted; extremes suppressed |

### Choosing `B`

| `B`       | When appropriate               |
| --------- | ------------------------------ |
| `200–500` | Quick exploration, prototyping |
| `2000`    | Standard reporting             |
| `5000+`   | Publication-quality CIs        |

---

## Examples

### Basic polarization score

```r
library(bipolarizationR)

result <- calculate_scores(df, "idemus", "idekemalist", score_type = "polarization")

cat(sprintf("Overall polarization: %.4f\n", result$point$overall))
cat(sprintf("95%% CI: [%.4f, %.4f]\n",
            result$bootstrap$overall$ci[1],
            result$bootstrap$overall$ci[2]))
cat(sprintf("SE: %.4f\n", result$bootstrap$overall$se))
```

### Consensus score with gaussian kernel

```r
result <- calculate_scores(
  df, "ideleft", "ideright",
  score_type = "consensus",
  kernel     = "gaussian",
  p          = 0.5,
  q          = 1.0,
  B          = 2000
)
result$point$overall
```

### Inspect the weight matrix before scoring

```r
# See how weights look for polarization with a steep distance curve
weight_matrix_visualization(p = 2.0, q = 1.0, type = "polarization", kernel = "power")

# Compare with gaussian kernel
weight_matrix_visualization(p = 0.3, q = 1.0, type = "polarization", kernel = "gaussian")
```

### Full dashboard

```r
dash <- dashboard_pair(df, "idemus", "idekemalist",
                       score_type = "polarization", B_boot = 2000)

# Access the underlying scores
dash$scores$point$overall
```

### Compare multiple pairs

```r
pairs <- list(
  c("idemus",    "ideathe"),
  c("idemus",    "idekemalist"),
  c("ideleft",   "ideright"),
  c("idesec",    "ideislamist")
)

results <- lapply(pairs, function(p) {
  r  <- calculate_scores(df, p[1], p[2], score_type = "polarization", B = 500)
  ci <- r$bootstrap$overall$ci
  cat(sprintf("%s x %s: %.4f  [%.4f, %.4f]\n",
              p[1], p[2], r$point$overall, ci[1], ci[2]))
  r
})
```

### Access raw bootstrap distributions

```r
result <- calculate_scores(df, "idemus", "idekemalist", keep_dists = TRUE)

dist <- result$bootstrap$overall$dist
hist(dist, breaks = 50, main = "Bootstrap distribution of overall polarization")
abline(v = result$point$overall, col = "red", lty = 2)
```

### Access the W×P matrix directly

```r
result <- calculate_scores(df, "idemus", "idekemalist")
wp     <- result$point$.wp_matrix
labels <- result$point$.x_labels

# Find the most polarizing cell
idx <- which(wp == max(wp), arr.ind = TRUE)
cat(sprintf("Highest cell: (%d, %d) = %.4f\n",
            labels[idx[1]], labels[idx[2]], wp[idx]))
```

---

<!-- ## Citation

If you use this package in academic work, please cite the underlying research:

> Nişancı, Z., Yüce, B., & Kavdır, A. Y. (2024). _bipolarizationR: Weighted kernel scoring for polarization and consensus in survey data._ R package version 0.1.0. -->
