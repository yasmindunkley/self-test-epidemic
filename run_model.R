get_ST_final <- function(R0, CFR, epsilon, D, sigma, sens, spec, c, r, time_end,
                             return_traj = FALSE, return_I = TRUE) {
  c_eff <- max(c, 1e-6)
  parameters <- c(R0=R0, CFR=CFR, epsilon=epsilon, D=D,
                  sigma=sigma, sens=sens, spec=spec, c=c_eff, r=r)
  times <- seq(0, time_end, by = 1)
  
  init <- c(S = 1 - 1e-6, E = 0, I = 1e-6, R = 0,
            Si = 0, Ei = 0, Ii = 0, Ri = 0,
            CI = 0, CM = 0, CIso = 0)
  
  result <- as.data.frame(deSolve::ode(
    y = init, times = times, func = self_test, parms = parameters,
    method = "lsoda", nonnegative = seq_along(init)
  ))
  
  final_row <- result[nrow(result), ]
  out <- list(
    final_CI   = final_row$CI,
    final_CM   = final_row$CM,
    final_S    = final_row$S,
    final_Si   = final_row$Si,
    final_Ii   = final_row$Ii,
    final_Ei   = final_row$Ei,
    final_R    = final_row$R,
    final_Ri   = final_row$Ri,
    final_CIso = final_row$CIso
  )
  if (return_I)    out$final_I <- final_row$I
  if (return_traj) out$traj    <- result
  out
}