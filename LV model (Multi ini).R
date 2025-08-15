library(deSolve)
library(tidyverse)

pop_sizes = function(times, states, parms){
  ### Build up Mechanism 2 ----
  M2 <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
      dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
      dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
      dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
      dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
      dS_dt = (r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S)# * epsilon
      return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
    })
  }
  ### Model parameters ----
  times <- seq(0, 8000, by = 0.1)
  parms <- c(epsilon = 1,
             r = 1, K = 10,
             a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
             b1 = 0.2, b2 = 0.45, m1 = 0.027, m2 = 0.001, e1H = 0.5, e2H = 0.5,
             o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)
  
  ### Model application (Run simulations)----
  ode(func = M2, times = times, y = states, parms = parms)
}

states_table = tibble(
  H = rep(0.7258128, 4),
  P1H = rep(0.3909225, 4),
  P2H = rep(0, 4),
  P1 = rep(2.2271095, 4),
  P2 = c(0.01, 0.1, 1, 10),
  S = rep(0.8368879, 4)
)

states_list = pmap(states_table, function(H, P1H, P2H, P1, P2, S){
  c(H = H, P1H = P1H, P2H = P2H, P1 = P1, P2 = P2, S = S)
}) #Put each row in the function which put each H, P1H, P2H, P1, P2 and S in the row into the vector that names with H, P1H, P2H, P1, P2 and S.

pop_sizes_list = 
  map(states_list, ~ pop_sizes(times, .x, parms))
  #using map to put states in to the function pop_sizes .x -> A list or atomic vector.
