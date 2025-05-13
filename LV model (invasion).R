library(deSolve)
library(tidyverse)

### M2 Invasion model
### Model specification
M2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
    dS_dt = r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S 
    return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
  })
}

### Event function
### P2 invade P1 monoculture
P2_invade = function(times, state, parms){
  with(as.list(c(state, parms)), {
    P1 = P1
    P2 = P2 + 1
    return(c(P1, P2))
  })
}

### Model parameters
times = seq(0, 250000, by = 0.1)
state = c(H = 0.1, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.2)
parms = c(r = 1, K = 10,
          a1 = 0.2, a2 = 0.1, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
          b1 = 0.35, b2 = 0.05, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
          o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)


### Model application w/ event function
### P2 invade P1-monoculture
pop_size_1 <- ode(func = LV, times = times,
                  y = state_1, parms = parms,
                  events = list(func = eventfun_2invade, time = invasion))
  
### Data manipulation
Data <- as.data.frame(rbind(pop_size_1, pop_size_2))
Data$Scenario <-rep(c("N2 invade N1-monoculture", "N1 invade N2-monoculture"),
                    each = length(times))

