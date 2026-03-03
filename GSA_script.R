rm(list = ls())

list.files("R", full.names = TRUE) |>
  lapply(source)

set.seed(1234)

# ========== Master simulation design (50 for trial - increase to 20,000 per manuscript methods) ==========
samples <- 50  
n_sims <- 1       
r_values <- NULL  

# ---------- global thresholds ----------
YEAR_CAP_DAYS <- 1825
TIME_FACTOR <- 150
MINDAY_FACTOR <- 10
THR_DAILY <- 1e-6
CI_SMALL_THR <- 0.01
K_SMOOTH <- 7

param_ranges <- list(
  R0 = c(1.1, 10),
  CFR = c(0.01, 0.5),
  epsilon = c(1/12, 1),
  D = c(3, 21),
  sigma = c(0, 0.15),
  sens = c(0.1, 1),  # placeholders – will be replaced by SRoC posterior
  spec = c(0.1, 1),  # placeholders – will be replaced by SRoC posterior
  c = c(1e-6, 0.8),  
  r = c(0, 0.15)
)

# ---------- sampling ----------
lhs_raw <- randomLHS(samples, length(param_ranges)) %>% as.data.frame()
names(lhs_raw) <- names(param_ranges)

lhs_parameters <- as.data.frame(matrix(nrow = samples, ncol = length(param_ranges)))
colnames(lhs_parameters) <- names(param_ranges)

# scale except sens/spec (will overwrite below)
for (param in setdiff(names(param_ranges), c("sens","spec"))) {
  lhs_parameters[[param]] <-
    param_ranges[[param]][1] +
    (lhs_raw[[param]] *
       (param_ranges[[param]][2] - param_ranges[[param]][1]))
}

# ---------- Build posterior_draws_all with samples = 20000 directly (via function) ----------
acc_obj <- draw_accuracy_reitsma(
  diagnostic_path = "data/diagnostic.csv",
  n_draws = samples,
  min_spec = 0.90,
  oversample_mult = 3,
  Sigma_source = "between",  
  seed = NULL
)

posterior_draws_all <- acc_obj$acc_draws

# overwrite sens/spec with posterior draws
lhs_parameters <- lhs_parameters %>%
  mutate(
    sens = posterior_draws_all$sens,
    spec = posterior_draws_all$spec
  )

# ------------------ MASTER SIM WITH HARD-STOP + NAMED OUTPUTS ------------------

simulation_results <- vector("list", n_sims)

for (sim in seq_len(n_sims)) {
  
  message("Running simulation ", sim, " of ", n_sims)
  
  df_sim <- lhs_parameters
  
  rows_out <- vector("list", nrow(df_sim))
  
  for (i in seq_len(nrow(df_sim))) {
    
    p <- df_sim[i, ]
    
    L_val <- 1 / p$epsilon
    min_day_val <- ceiling(MINDAY_FACTOR * (L_val + p$D))
    max_day_val <- min(YEAR_CAP_DAYS,
                       ceiling(TIME_FACTOR * (L_val + p$D)))
    
    base <- get_ST_final(
      R0=p$R0, CFR=p$CFR, epsilon=p$epsilon, D=p$D,
      sigma=0, sens=p$sens, spec=p$spec,
      c=p$c, r=p$r,
      time_end=max_day_val,
      return_traj=TRUE
    )
    
    df_b <- base$traj
    hs_b <- find_hard_stop(df_b,
                           thr_daily=THR_DAILY,
                           min_day=min_day_val,
                           k_smooth=K_SMOOTH)
    b_row <- df_b[df_b$time == hs_b, ]
    
    inter <- get_ST_final(
      R0=p$R0, CFR=p$CFR, epsilon=p$epsilon, D=p$D,
      sigma=p$sigma, sens=p$sens, spec=p$spec,
      c=p$c, r=p$r,
      time_end=max_day_val,
      return_traj=TRUE
    )
    
    df_i <- inter$traj
    hs_i <- find_hard_stop(df_i,
                           thr_daily=THR_DAILY,
                           min_day=min_day_val,
                           k_smooth=K_SMOOTH)
    i_row <- df_i[df_i$time == hs_i, ]
    
    deaths_averted <- b_row$CM - i_row$CM
    
    peak_prev_baseline <- max(df_b$I + df_b$Ii, na.rm = TRUE)
    peak_prev_interv   <- max(df_i$I + df_i$Ii, na.rm = TRUE)
    
    delta_peak_I_abs <- peak_prev_baseline - peak_prev_interv
    
    delta_peak_I_rel <- ifelse(peak_prev_baseline > 0,
                               delta_peak_I_abs / peak_prev_baseline,
                               NA_real_)
    
    t_peak_b <- df_b$time[which.max(df_b$I + df_b$Ii)]
    t_peak_i <- df_i$time[which.max(df_i$I + df_i$Ii)]
    
    df_i <- df_i %>%
      dplyr::mutate(
        daily_tests = pmax(p$sigma * (S + E + I + R), 0)
      )
    
    tests_total_hs <- sum(
      df_i$daily_tests[df_i$time <= hs_i],
      na.rm = TRUE
    )
    
    tests_per_death_averted_hs <- ifelse(deaths_averted > 0,
                                         tests_total_hs / deaths_averted,
                                         NA_real_)
    
    tests_per_case_averted_hs <- ifelse((b_row$CI - i_row$CI) > 0,
                                        tests_total_hs /
                                          (b_row$CI - i_row$CI),
                                        NA_real_)
    
    rows_out[[i]] <- tibble(
      R0 = p$R0,
      CFR = p$CFR,
      epsilon = p$epsilon,
      D = p$D,
      sigma = p$sigma,
      sens = p$sens,
      spec = p$spec,
      c = p$c,
      r = p$r,
      time_end = max_day_val,
      min_day = min_day_val,
      baseline_end_day = hs_b,
      intervention_end_day = hs_i,
      time_to_end_diff = hs_i - hs_b,
      baseline_cases_end = b_row$CI,
      baseline_deaths_end = b_row$CM,
      baseline_isolations_end = b_row$CIso,
      intervention_cases_end = i_row$CI,
      intervention_deaths_end = i_row$CM,
      intervention_isolations_end = i_row$CIso,
      cases_averted = b_row$CI - i_row$CI,
      deaths_averted = deaths_averted,
      peak_I_baseline = peak_prev_baseline,
      peak_I_intervention = peak_prev_interv,
      delta_peak_I_abs = delta_peak_I_abs,
      delta_peak_I_rel = delta_peak_I_rel,
      t_peak_b = t_peak_b,
      t_peak_i = t_peak_i,
      tests_total_hs = tests_total_hs,
      tests_per_death_averted = tests_per_death_averted_hs,
      tests_per_case_averted = tests_per_case_averted_hs,
      Simulation = sim
    )
  }
  
  df_out <- bind_rows(rows_out)
  
  simulation_results[[sim]] <- df_out
}

all_sims_df <- bind_rows(simulation_results)