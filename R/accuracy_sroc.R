draw_accuracy_reitsma <- function(
  diagnostic_path = "data/diagnostic.csv",
  n_draws = 20000,
  min_spec = 0.90,
  oversample_mult = 3,
  Sigma_source = c("vcov", "between"),
  seed = NULL
) {
  Sigma_source <- match.arg(Sigma_source)
  if (!is.null(seed)) set.seed(seed)

  # read data
  diagnostic <- readr::read_csv(diagnostic_path, show_col_types = FALSE)

  # basic checks
  req <- c("TP", "FP", "FN", "TN")
  if (!all(req %in% names(diagnostic))) {
    stop("diagnostic.csv must contain columns: TP, FP, FN, TN")
  }

  # coerce and drop invalid rows
  diagnostic <- diagnostic |>
    dplyr::mutate(
      TP = as.numeric(TP),
      FP = as.numeric(FP),
      FN = as.numeric(FN),
      TN = as.numeric(TN)
    ) |>
    dplyr::filter(
      !is.na(TP), !is.na(FP), !is.na(FN), !is.na(TN),
      TP >= 0, FP >= 0, FN >= 0, TN >= 0,
      (TP + FN) > 0,
      (FP + TN) > 0
    )

  # fit reitsma
  fit <- mada::reitsma(diagnostic, TP = "TP", FN = "FN", FP = "FP", TN = "TN")

  coeffs <- as.numeric(fit$coefficients)
  names(coeffs) <- c("logit_sens", "logit_fpr")

  # choose covariance for sampling
  Sigma <- NULL
  if (Sigma_source == "vcov") {
    # uncertainty in pooled mean (asymptotic)
    Sigma <- tryCatch(stats::vcov(fit), error = function(e) fit$vcov)
  } else {
    # between-study heterogeneity (if available in object)
    Sigma <- if (!is.null(fit$Psi)) {
      fit$Psi
    } else if (!is.null(fit$Sigma.b)) {
      fit$Sigma.b
    } else {
      stop("Between-study covariance not found in fit object (Psi/Sigma.b).")
    }
  }

  # oversample then reject on spec
  k <- as.integer(oversample_mult * n_draws)

  make_draws <- function(n) {
    z <- MASS::mvrnorm(n = n, mu = coeffs, Sigma = Sigma)
    out <- data.frame(
      sens = plogis(z[, "logit_sens"]),
      spec = 1 - plogis(z[, "logit_fpr"])
    )
    out[out$spec >= min_spec, , drop = FALSE]
  }

  acc <- make_draws(k)

  while (nrow(acc) < n_draws) {
    need <- n_draws - nrow(acc)
    extra <- make_draws(max(2L * need, 1000L))
    acc <- dplyr::bind_rows(acc, extra)
  }

  acc <- acc[seq_len(n_draws), , drop = FALSE]

  # return both draws and the fitted object for traceability
  list(
    acc_draws = acc,
    fit = fit,
    Sigma_source = Sigma_source,
    min_spec = min_spec
  )
}