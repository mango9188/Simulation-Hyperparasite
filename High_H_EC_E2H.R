# For H to change the equilibrium + Invasion
# In order to know how many H is needed to invade the system.
library(deSolve)
library(tidyverse)


times <- seq(0, 20000, by = 0.01)
ini_High_H <- c(H = 0.5, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.5) #ini_1
parms_EC_E2H =
  c(r = 1, K = 10,
    a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
    b1 = 0.2, b2 = 0.42, m1 = 0.06, m2 = 0.01, e1H = 0.5, e2H = 0.5,
    o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)

TimeSeries = function(times, state, parms){
  M2 <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
      dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
      dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
      dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
      dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
      dS_dt = (r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S)
      return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
    })
  }
  pop_size = ode(func = M2, times = times, y = state, parms = parms)
  pop_size %>%
    as.data.frame() #%>%
  #filter(time %% 1 == 0)
}
High_H_EC_E2H = TimeSeries(times, ini_High_H, parms_EC_E2H) #E_C



# Run a for-loop with bisection method ----
#Define the function
times <- seq(0, 20000, by = 10)
M2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
    dS_dt = (r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S)
    return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
  })
}

### Read the time series from previous simulation to further subset different initial states

#function that will output the competition outcome
f_outcome = function(state, parms){
  pop_size = ode(func = M2,
                 times = times,
                 y = state,
                 parms = parms,
                 maxsteps = 20000)
  Outcome =
    paste(ifelse(pop_size[nrow(pop_size), c("H", "P1H", "P2H", "P1", "P2")] > 1e-7, "T", "F"), collapse = "")
  return(Outcome)
}

#function of bisection method
Bisection_H = function(min_H, max_H, state, parms){
  ###The goal is to find a maximum value of H that make competition outcome == E2H
  ###Min
  Min_state = state
  Min_state["H"] = min_H
  min_E = f_outcome(state = Min_state, parms = parms)
  
  ###Max
  Max_state = state
  Max_state["H"] = max_H
  max_E = f_outcome(state = Max_state, parms = parms)
  
  if(min_E == max_E){
    print(paste0("The equilibrium of minimum H is equal to maximum H (", min_H, ")"))
    return(NA)
  }
  
  loop_count = 0
  while(max_H - min_H > 1e-6 && loop_count < 80){
    ###Mid
    mid_H = (min_H + max_H)/2
    Mid_state = state
    Mid_state["H"] = mid_H
    mid_E = f_outcome(state = Mid_state, parms = parms)
    
    ###determine the value of mid_H
    if(mid_E == min_E){
      min_H = mid_H
    }else if(mid_E == max_E){
      max_H = mid_H
    }else{
      return(NA)
    }
    loop_count = loop_count + 1
  }
  return((min_H + max_H)/2)
}

#The timing of invasion
x = seq(0, 5, by = 0.01)
Invasion_timing = 10^x
#dim(High_H_EC_E2H)[1]


#Set up
High_H_EC_E2H$H_press = 0

#Run a for-loop
Start_time = Sys.time()
for (i in Invasion_timing) {
  #setting the initial condition
  #temp_state = round(unlist(High_H_EC_E2H[i,  c("H", "P1H", "P2H", "P1", "P2", "S")]), 8) #use unlist() to make values become vectors
  temp_state = round(unlist(High_H_EC_E2H[round(High_H_EC_E2H$time, 2) == round(i, 2),  c("H", "P1H", "P2H", "P1", "P2", "S")]), 8)
  
  #run with the bisection method
  High_H_EC_E2H[i,]$H_press = 
    Bisection_H(min_H = High_H_EC_E2H[i, "H"],
                max_H = 25,
                state = temp_state,
                parms = parms_EC_E2H)
  print(paste0("Total: ", length(Invasion_timing), ", Now: ", which(Invasion_timing == i)))
}
End_time = Sys.time()
End_time - Start_time

#saveRDS(High_H_EC_E2H, "High_H_EC_E2H_invadeTest1")
