# CLAUDE.md

Guidance for Claude Code (and future agents) working in the **LatencyA** R package.

## What this package does

Latency analysis in epidemiological studies via **weighted cumulative exposure**
Cox models. Splines (or polynomials) weight time-varying exposures by their timing
relative to the outcome; that weight function is what estimates the "latency."

## Source layout (`R/`)

- `CoxNCSpline.R`  — Cox model with natural cubic splines + `extract_CoxNCSpline`
- `CoxPoly.R`      — Cox model with polynomial terms + `extract_CoxPoly`
- `CoxKnotSearch.R`— knot-number search + `extract_CoxKnotsearch` + `CoxKnotsearch_boot`
- `CoxTermsearch.R`— term/model selection + `extract_CoxTermsearch` + `CoxTermsearch_boot`
- `resample_clusters.R` — cluster (id) resampling used by the bootstraps

## The four `extract_*` functions — conventions

All four accept `lag` as **a scalar or a vector** and return **one row per lag**.

- The basis row per lag is selected with `match(lag, 0:(latency-1))` (or
  `outer(lag, 0:K, "^")` for the polynomial), always with `drop = FALSE` so the
  single-lag case is just the 1-row case — no `if (length(lag) > 1)` branching.
- Per-lag log HR: `as.vector(B_lag %*% coef)`.
- Per-lag variance (only the standalone pair): `rowSums((B_lag %*% vcov) * B_lag)`
  — the diagonal of `B_lag %*% vcov %*% t(B_lag)` without forming the full matrix.

Return shapes differ **by design**:

| function                | columns                      | has analytic variance? |
|-------------------------|------------------------------|------------------------|
| `extract_CoxNCSpline`   | `lag, log_HR, log_HR_var`    | yes                    |
| `extract_CoxPoly`       | `lag, log_HR, log_HR_var`    | yes                    |
| `extract_CoxKnotsearch` | `lag, log_HR`                | no — variance is bootstrapped |
| `extract_CoxTermsearch` | `lag, log_HR`                | no — variance is bootstrapped |

> History: these used to collapse a vector of lags into a single **cumulative**
> estimate via `colSums(B_lag)`. That behavior was removed in favor of per-lag
> output. If cumulative is ever needed again, it is `colSums(B_lag)` on the same
> basis — add it as an explicit `cumulative = TRUE` arg rather than overloading
> vector input.

## Bootstrap contract (`*_boot`)

`CoxKnotsearch_boot` / `CoxTermsearch_boot` are always called with a **single**
`lag`. They consume the extractor via `extract_*(...)$log_HR` (length-1) inside
`log_HR[i] <-` and `furrr::future_map_dbl`, both of which require a scalar per
call. If you change an extractor's return shape, keep `$log_HR` length-1 for the
single-lag path or you will break these.

Variance/percentiles come from the bootstrap distribution, not analytically.

### Parallelism

- `parallel = FALSE` → serial `for` loop.
- `parallel = TRUE`  → `furrr::future_map_dbl` under `plan(parallel_plan)`.
- `parallel_plan = c("multicore", "multisession")`, validated with `match.arg`;
  **default `"multicore"`**. Note `multicore` does not work on Windows / inside
  RStudio — `multisession` is the portable fallback.

## Testing gotcha

`devtools::load_all()` + `multisession` will fail: multisession spawns fresh R
processes that load the **installed** package, not the `load_all` dev code. To
test multisession changes, `devtools::install()` first, then `library(LatencyA)`.
`multicore` forks the current process, so it sees `load_all` code fine.

## Routine after editing `R/`

Run `devtools::document()` to regenerate `man/*.Rd` and `NAMESPACE`.

## Packaging

`.claude/` and `CLAUDE.md` are excluded from the build via `.Rbuildignore`
(this repo targets CRAN). `.claude/settings.local.json` is gitignored (local-only).
