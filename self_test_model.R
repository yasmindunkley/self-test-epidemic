
self_test <- function(time, state, parameters, ...) {
  with(as.list(c(parameters, state)), {
    
    # Derived parameters
    mu <- CFR / D
    gamma <- mu * (1 - CFR) / CFR
    beta <- R0 * (gamma + mu)
    N <- S + E + I + R + Si + Ei + Ii + Ri
    d <- 1.5 * D
    delta <- 1 / (c * d)  
    beta_i <- beta * (1 - c)  
    foi <- ((beta * I) / N) + ((beta_i * Ii) / N)
    
    ## differentials
    dS <- (delta * Si) -(foi * S) -(sigma * ((1 - spec) * S)) 
    dE <- (foi * S) -(sigma * ((1 - spec) * E)) -(epsilon * E) +(delta * Ei)
    dI <- (epsilon * E) -(sigma * (sens * I)) -(mu * I) -(gamma * I) -(r * I) +(delta * Ii)
    dR <- (gamma * I) -(sigma * ((1 - spec) * R)) + (delta * Ri)
    dSi <- (sigma * ((1 - spec) * S)) - (delta * Si) - ((1 - c) * foi * Si)
    dEi <- (sigma * ((1 - spec) * E)) - (epsilon * Ei) + ((1 - c) * foi * Si) - (delta * Ei)
    dIi <- (epsilon * Ei) + (sigma * (sens * I)) - (mu * Ii) - (gamma * Ii) - (delta * Ii) +(r * I)
    dRi <- (gamma * Ii) + (sigma * ((1 - spec) * R)) - (delta * Ri)
    dCI <- (epsilon * E) + (epsilon * Ei)
    dCM <- (mu * I) + (mu * Ii)
    
    dCIso <- sigma * ((1 - spec) * (S + E + R) + sens * I) + (r * I)
    
    # Return the list of derivatives
    list(c(dS, dE, dI, dR, dSi, dEi, dIi, dRi, dCI, dCM, dCIso))
  })
}
