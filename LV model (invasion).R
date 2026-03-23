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

### Event function----
#### P2 invade P1 monoculture----
P2_invade = function(times, state, parms){
  with(as.list(c(state, parms)), {
    P1 = P1
    P2 = P2 + 1
    return(c(P1, P2))
  })
}
### H invasion-----
H_invasion = function(times, state, parms){
  with(as.list(c(state, parms)), {
    H = H + 8
    P1H = P1H
    P2H = P2H
    P1 = P1
    P2 = P2
    S = S
    return(c(H, P1H, P2H, P1, P2, S))
  })
}

.
#time-----
times = seq(0, 7000, by = 0.1)

#state----
ini_High_H <- c(H = 1, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.5) #ini_1 (P1win)
ini_Low_H <- c(H = 0.01, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.5) #ini_2 (P2win)

#parms----
# parms_E1H_E2H = 
#   c(epsilon = 1,
#     r = 1, K = 10,
#     a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
#     b1 = 0.2, b2 = 0.45, m1 = 0.03, m2 = 0.001, e1H = 0.5, e2H = 0.5,
#     o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)
# parms_EC_E2H = 
#   c(epsilon = 1,
#     r = 1, K = 10,
#     a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
#     b1 = 0.2, b2 = 0.45, m1 = 0.045, m2 = 0.001, e1H = 0.5, e2H = 0.5,
#     o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

parms_EC_E2H <- c(r = 1, K = 10,
                  a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
                  b1 = 0.2, b2 = 0.42, m1 = 0.06, m2 = 0.005, e1H = 0.5, e2H = 0.5,
                  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)
parms_E1H_E2H <- c(r = 1, K = 10,
                  a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
                  b1 = 0.2, b2 = 0.42, m1 = 0.04, m2 = 0.005, e1H = 0.5, e2H = 0.5,
                  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)



#Plotting----
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
    as.data.frame() %>%
    filter(time %% 1 == 0) %>%
    #filter(time <= 500 | time >= 4500) %>% 
    #mutate(Time_Window = ifelse(time <= 500, "Transient", "Equilibrium")) %>%
    #mutate(Time_Window = factor(Time_Window, levels = c("Transient", "Equilibrium"))) %>%
    pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"),
                 names_to = "species", values_to = "biomass") %>%
    ggplot(mapping = aes(x = time, y = biomass, color = species)) +
    labs(x = "Time", y = "Abundance") +
    geom_line(lwd = 1) +
    #facet_wrap(~ Time_Window, scales = "free_x") +
    scale_y_continuous(limits = c(0, 9.5))+
    scale_colour_manual(labels = State_labels,
                        values = State_values)+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.spacing = unit(0.2, "lines"))
}
TimeSeriesInvasion = function(times, state, parms, InvasionTimes){
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
  pop_size = ode(func = M2, times = times, y = state, parms = parms,
                 events = list(func = H_invasion, time = InvasionTimes))
  pop_size %>%
    as.data.frame() %>%
    # filter(time %% 1 == 0) %>%
    # filter(time <= 500 | time >= 4500) %>% 
    # mutate(Time_Window = ifelse(time <= 500, "Transient", "Equilibrium")) %>%
    # mutate(Time_Window = factor(Time_Window, levels = c("Transient", "Equilibrium"))) %>%
    pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"),
                 names_to = "species", values_to = "biomass") %>%
    ggplot(mapping = aes(x = time, y = biomass, color = species)) +
    labs(x = "Time", y = "Abundance") +
    geom_line(lwd = 1) +
    #facet_wrap(~ Time_Window, scales = "free_x") +
    scale_y_continuous(limits = c(0, 9.5))+
    scale_colour_manual(labels = State_labels,
                        values = State_values)+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.spacing = unit(0.2, "lines"))
}

TimeSeries(times, ini_High_H, parms_EC_E2H) #Coexist
TimeSeries(times, ini_Low_H, parms_EC_E2H)  #E2H

TimeSeries(times, ini_High_H, parms_E1H_E2H) #E1H
TimeSeries(times, ini_Low_H, parms_E1H_E2H) #E2H

TimeSeriesInvasion(times, ini_High_H, parms_EC_E2H, InvasionTimes = c(4000)) #Coexist -> EP2H
TimeSeriesInvasion(times, ini_Low_H, parms_EC_E2H, InvasionTimes = c(10)) #EP2H -> Coexist

TimeSeriesInvasion(times, ini_High_H, parms_E1H_E2H, InvasionTimes = c(150)) #EP1H -> EP2H (if H = H+5)
TimeSeriesInvasion(times, ini_Low_H, parms_E1H_E2H, InvasionTimes = c(1)) #EP2H

ggsave("ASS time series ini_Low_H parms_E1H_E2H Int2000.png", width = 16, height = 11, units = "cm", dpi = 800)
ggsave("ASS time series ini_High_H parms_EC_E2H.png", width = 25, height = 11, units = "cm", dpi = 800)
ggsave("ASS time series ini_Low_H parms_EC_E2H.png", width = 25, height = 11, units = "cm", dpi = 800)
ggsave("ASS time series ini_High_H parms_E1H_E2H.png", width = 25, height = 11, units = "cm", dpi = 800)
ggsave("ASS time series ini_Low_H parms_E1H_E2H.png", width = 25, height = 11, units = "cm", dpi = 800)

ggsave("ASS time series ini_High_H parms_E1H_E2H EP1H to EP2H.png", width = 25, height = 11, units = "cm", dpi = 800)
ggsave("ASS time series ini_High_H parms_EC_E2H EC to EP2H.png", width = 25, height = 11, units = "cm", dpi = 800)

ggsave("Invasion_NEW.png", width = 25, height = 11, units = "cm", dpi = 800)
