# --- hard stop using trailing smoothing on daily CI ---
# hard stop = last day when smoothed daily incidence has been ≤ threshold for the past 100 days
## after explicitly the lastmin window

find_hard_stop <- function(df, thr_daily = 1e-6, min_day = 0, k_smooth = 7) {
  daily <- c(NA_real_, diff(df$CI))
  sm <- stats::filter(daily, rep(1 / k_smooth, k_smooth), sides = 1, method = "convolution")
  sm <- as.numeric(sm)
  
  ok <- sm <= thr_daily
  ok[is.na(ok)] <- FALSE
  
  tvec <- df$time
  win_days <- 100L
  i_min <- which(tvec >= min_day)[1]; if (is.na(i_min)) i_min <- length(tvec)
  
  start <- max(win_days + 1L, 2L, i_min)
  hs <- NA_integer_
  for (k in start:length(tvec)) {
    idx <- (k - win_days + 1L):k
    if (min(idx) < 1L) next
    # enforce window STARTS after min_day:
    if ((k - win_days + 1L) < i_min) next
    if (!all(ok[idx])) next
    hs <- tvec[k]                  # end of first quiet 100-day streak
    # if you prefer the FIRST day of the streak: hs <- tvec[k - win_days + 1L]
    break
  }
  if (is.na(hs)) hs <- tail(tvec, 1)
  hs
}

