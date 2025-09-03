library(tidyverse)
library(deSolve)

#Parameter setting
parms <- c(
           r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

times <- seq(0, 100000, by = 1)

state <- c(H = 0.1, P1H = 0, P2H = 0, P1 = 0.1, P2 = 0.1, S = 0.2)

#Function setting
dP1_dt = function(times, state, parms){
  M2 <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
      dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
      dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
      dP1_dt = 0
      dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
      dS_dt = (r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S)
      return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
    })
    }
  pop_size = ode(func = M2, times = times, y = state, parms = parms)
  pop_size %>% as.data.frame()
  parms[c('H', 'P1H', 'P2H', 'P1', 'P2', 'S')] = 
    pop_size[length(times), c("H", "P1H", "P2H", "P1", "P2", "S")]
  #return(quasiequ)
  with(as.list(parms), {
    (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
  })
}
#Need to solve cycling problem first!
dP1_dt(times, state, parms)


P1 = seq(0.01, 0.02, by = 0.01)

dP1_dt_P1 = as.data.frame(cbind(P1,
                                matrix(0, 
                                       nrow = length(P1),
                                       ncol = 1)))
names(dP1_dt_P1) = c("P1", "dP1_dt")

###For-loop
for (i in 1:length(P1)) {
  temp_state = state
  temp_state["P1"] = P1[i]
  growthRate_P1 = dP1_dt(times, temp_state, parms)
  dP1_dt_P1[i, "dP1_dt"] = growthRate_P1
}
dP1_dt_P1
