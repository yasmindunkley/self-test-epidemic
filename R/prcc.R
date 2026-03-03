# ------------------ PARTIAL CORRELATION RANKING ------------------

predictor_vars <- c("R0","CFR","epsilon","D","sigma","sens","spec","c","r")
outcome_vars   <- c("delta_peak_I_rel")

# ------------------ PARTIAL CORRELATION RANKING FUNCTION ------------------
partial_cor_rank <- function(df, y_name, x_names) {
  # pick only the needed columns, drop rows with any NA in those cols
  subdf_all <- df %>%
    dplyr::select(dplyr::all_of(c(x_names, y_name))) %>%
    dplyr::filter(stats::complete.cases(.))
  
  # if too few rows/cols to run partial corr, return NA row per predictor
  if (nrow(subdf_all) < 3L) {
    return(
      tibble::tibble(
        var = x_names,
        pcor = NA_real_,
        p.value = NA_real_,
        n = nrow(subdf_all),
        rank_abs = NA_integer_,
        outcome = y_name,
        sign = NA_character_
      )
    )
  }
  
  # ensure everything is numeric for ppcor
  subdf_all <- dplyr::mutate(subdf_all, dplyr::across(dplyr::everything(), as.numeric))
  
  res <- lapply(x_names, function(xi) {
    others <- setdiff(x_names, xi)
    subdf  <- subdf_all[, c(xi, others, y_name), drop = FALSE]
    
    pc <- try(ppcor::pcor(subdf, method = "spearman"), silent = TRUE)
    
    if (inherits(pc, "try-error")) {
      tibble::tibble(var = xi, pcor = NA_real_, p.value = NA_real_, n = nrow(subdf))
    } else {
      tibble::tibble(
        var = xi,
        pcor = unname(pc$estimate[xi, y_name]),
        p.value = unname(pc$p.value[xi, y_name]),
        n = nrow(subdf)
      )
    }
  })
  
  dplyr::bind_rows(res) %>%
    dplyr::arrange(dplyr::desc(abs(.data$pcor))) %>%
    dplyr::mutate(
      rank_abs = dplyr::row_number(),
      outcome  = y_name,
      sign     = dplyr::case_when(
        is.na(.data$pcor) ~ NA_character_,
        .data$pcor >= 0   ~ "+",
        TRUE              ~ "−"
      )
    )
}

# --- run for each outcome ---
pcor_results <- purrr::map(
  outcome_vars,
  ~partial_cor_rank(all_sims_df, .x, predictor_vars)
)
names(pcor_results) <- outcome_vars

# one combined table
pcor_results_tbl <- purrr::list_rbind(pcor_results)

pcor_results_tbl