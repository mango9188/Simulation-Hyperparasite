# For H to change the equilibirum
library(deSolve)
library(tidyverse)
library(pracma)

1
# Perform time scale separation, see how big H is to cross basin of attraction----
# Here, I set H as a constant
### Build up model with constant H
Hcon <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
    dS_dt = r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S 
    return(list(c(dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
  })
}

#Model parameters
times <- seq(0, 10000, by = 0.1)
state <- c(P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.2)
parms <- c(H = 9999,
           r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.01, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

H_press = seq(0.1, 1, by = 0.1)
Data = data.frame(H_press,
                  matrix(0, nrow = length(H_press), ncol = length(state)))
names(Data) = c("H_press", "P1H", "P2H", "P1", "P2", "S")

start = Sys.time()
for (i in 1:dim(Data)[1]) {
  temp_parms = parms
  temp_parms["H"] = Data$H_press[i]
  quasi_pop_size = ode(func = Hcon,
                       times = times,
                       y = state,
                       parms = temp_parms)
  
  if (sd(quasi_pop_size[(nrow(quasi_pop_size)-round(length(times)*0.3)):nrow(quasi_pop_size),"S"]) > 1e-8){
    paste0("Cycles happen when H = ", temp_parms["H"])
    k = 1
    avg = c(rep(NA, 5)) #create a vector to store the long term average of each species 
    for (j in names(Data)[2:6]) { #run a loop that can input each species
      peaks = findpeaks(quasi_pop_size[, j])[ ,2]
      periods = peaks[length(peaks)] - peaks[length(peaks) - 1]
      avg[k] = mean(quasi_pop_size[(length(times) - periods + 1):length(times), j]) #store the long term average in to the vector that just created
      k = k + 1 #increase k so that we can put the value of next species into the vector
    }
    k = 0 #reset k (it's unnecessary actually)
    Data[i, 2:6] = avg
  }else{
    Data[i, 2:6] = quasi_pop_size[nrow(quasi_pop_size), -1]
  }
}
end = Sys.time()
end-start


Data = mutate(Data, dH_dt = with(as.list(parms), {
  (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H_press - (d * H_press)}))

ggplot(d = Data, mapping = aes(x = H_press, y = dH_dt))+
  geom_line()

View(Data)

Data = readRDS(file = "m1 = 005, m2 = 001 from Ec to E2H")

# Even though we found the basin of attraction, we can not make use of it. More H is needed to make E_C become E_2H.


# How large an H is required to shift the system's equilibrium point?-----------
### Build up Mechanism 2
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

### Model parameters
times <- seq(0, 10000, by = 0.1)
state <- c(H = 888888, P1H = 0.2333761, P2H = 0.1687427, P1 = 1.4767783, P2 = 0.4522385, S = 0.9095539)

parms <- c(r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.01, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

H_press = seq(0.5, 5, by = 0.5)
comp_out = data.frame(H_press,
                  matrix(0, nrow = length(H_press), ncol = length(state) + 1))
names(comp_out) = c("H_press", "H", "P1H", "P2H", "P1", "P2", "S", "Equilibrium")

start = Sys.time()
for (i in 1:dim(comp_out)[1]) {
  temp_state = state
  temp_state["H"] = comp_out$H_press[i]
  pop_size = ode(func = M2,
                       times = times,
                       y = temp_state,
                       parms = parms)
  
  if (sd(pop_size[(nrow(pop_size)-round(length(times)*0.3)):nrow(pop_size),"S"]) > 1e-8){
    paste0("Cycles happen when H = ", temp_state["H"])
    k = 1
    avg = c(rep(NA, 6)) #create a vector to store the long term average of each species 
    for (j in names(Data)[2:7]) { #run a loop that can input each species
      peaks = findpeaks(pop_size[, j])[ ,2]
      periods = peaks[length(peaks)] - peaks[length(peaks) - 1]
      avg[k] = mean(pop_size[(length(times) - periods + 1):length(times), j]) #store the long term average in to the vector that just created
      k = k + 1 #increase k so that we can put the value of next species into the vector
    }
    k = 0 #reset k (it's unnecessary actually)
    comp_out[i, 2:7] = avg
  }else{
    comp_out[i, 2:7] = pop_size[nrow(pop_size), -1]
  }
}
end = Sys.time()
end-start

### Data analysis
extinct_thres = 1e-7

comp_out$Outcome =
  apply(comp_out[, c("H", "P1H", "P2H", "P1", "P2")], 1, function(row){
    paste(ifelse(row > extinct_thres, "T", "F"), collapse = "")
  })


# Run a for-loop to apply it to fit parameter sets in each grids---------
### Build up Mechanism 2
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

### Read the data to further subset different initial states and parameters
P.space = 
  readRDS("Pre4A1_forPSpace") %>%
  filter(Stable_E == "P2H,C")

###Create a function to determine the range of H_press

### Model parameters
times <- seq(0, 10000, by = 1)
state <- c(H = 0, P1H = 0, P2H = 0, P1 = 0, P2 = 0, S = 0)

parms <- c(r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

H_press = seq(1, 6, by = 1)
parameter = expand_grid(P.space[, c("m1", "m2", "H", "P1H", "P2H", "P1", "P2", "S")], H_press)
parameter = mutate(parameter, H_press = H_press + H)
H_series = data.frame(parameter,
                      matrix(0, nrow = length(H_press), ncol = length(state)))
names(H_series) = c("m1", "m2", "H", "P1H", "P2H", "P1", "P2", "S", "H_press", "H_t", "P1H_t", "P2H_t", "P1_t", "P2_t", "S_t")


start = Sys.time()
for (i in 1:dim(H_series)[1]) {
  temp_state = state
  temp_state[c("H", "P1H", "P2H", "P1", "P2", "S")] = 
    H_series[i, c("H_press", "P1H", "P2H", "P1", "P2", "S")] %>%
    unlist()
  temp_parms = parms
  temp_parms["m1"] = H_series$m1[i]
  temp_parms["m2"] = H_series$m2[i]
  
  pop_size = ode(func = M2,
                 times = times,
                 y = temp_state,
                 parms = temp_parms)
  H_series[i, 10:15] = pop_size[nrow(pop_size), -1]
  
}

end = Sys.time()
end-start

### Data analysis
extinct_thres = 1e-7

H_series$Outcome =
  apply(H_series[, c("H_t", "P1H_t", "P2H_t", "P1_t", "P2_t")], 1, function(row){
    paste(ifelse(row > extinct_thres, "T", "F"), collapse = "")
  })
H_series =
  H_series %>%
  mutate(delta_H = H_press - H)

saveRDS(H_series, file = "H to change the equilibirum from Pre4A1 Press1")

# Run a for-loop with bisection method ----
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

### Read the data to further subset different initial states and parameters
P.space = 
  readRDS("Pre4A1_forPSpace") %>%
  filter(Stable_E == "P2H,C")
P.space = P.space[,-c(9:11)]
P.space$H_press = 0
head(P.space)

### Model parameters
times <- c(0, 10000)
state <- c(H = 0, P1H = 0, P2H = 0, P1 = 0, P2 = 0, S = 0)

parms <- c(r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

# parameter = expand_grid(P.space[, c("m1", "m2", "H", "P1H", "P2H", "P1", "P2", "S")])

#function that will output the competition outcome
f_outcome = function(state, parms){
  pop_size = ode(func = M2,
                 times = times,
                 y = state,
                 parms = parms,
                 maxsteps = 50000)
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
  }
  
  loop_count = 0
  while(max_H - min_H > 1e-9 && loop_count < 100){
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

Start_time = Sys.time()
for (i in 1:dim(P.space)[1]) {
  #setting the initial condition
  temp_state = unlist(P.space[i,  c("H", "P1H", "P2H", "P1", "P2", "S")])
  
  #setting the parameters
  temp_parms = parms
  temp_parms[c("m1", "m2")] = P.space[i, c("m1", "m2")]
  
  #run with the bisection method
  P.space[i,]$H_press = Bisection_H(min_H = min(P.space[i, "H"]),
                                    max_H = 25,
                                    state = temp_state,
                                    parms = temp_parms)
  if(i %% 10 == 0) print(paste0("Total: ", nrow(P.space), ", Now: ", i))
}
Ending_time = Sys.time()
Ending_time - Start_time
