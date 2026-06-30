library(ggplot2)
library(tidyverse)
library(deSolve)
library(gganimate)
library(gifski)
#Theme setting======
C = 
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.y.right = element_text(size = 18),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.text.align = 0,
    #panel.border = element_rect(color = "grey", fill = NA, size = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
theme_set(C)

# Result 4_1: the time series invasion/ parameter switching ----
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


## Pesticide application (i.e. changing parameters m in a specific time)----
times1 = seq(0, 500, by = 1)
times2 = seq(500, 3000, by = 0.1)

Initial = c(H = 0.7383519, P1H = 0.2883162, P2H = 0.1024059, P1 = 1.67909, P2 = 0.2674836, S = 0.9609926)
parms_Low_m1 <- c(r = 1, K = 10,
                  a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
                  b1 = 0.2, b2 = 0.42, m1 = 0.06, m2 = 0.01, e1H = 0.5, e2H = 0.5,
                  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)

parms_High_m1 <- c(r = 1, K = 10,
               a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
               b1 = 0.2, b2 = 0.42, m1 = 0.07, m2 = 0.01, e1H = 0.5, e2H = 0.5,
               o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)

TimeSeriesdisturbance = function(times1, times2, state, parms1, parms2){
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
  pop_size1 = 
    ode(func = M2, times = times1, y = state, parms = parms1)
  
  new_state = 
    pop_size1 %>%
    tail(1) %>%
    unlist()
  new_state = new_state[2:7]
  names(new_state) = c("H", "P1H", "P2H", "P1", "P2" , "S")
  
  pop_size2 = 
    ode(func = M2, times = times2, y = new_state, parms = parms2)
  
  pop_size1 = 
    pop_size1 %>%
    as.data.frame()
  
  pop_size2 = 
    pop_size2 %>%
    as.data.frame()
  
  rbind(pop_size1, pop_size2) %>%
    # filter(time %% 1 == 0) %>%
    # filter(time <= 500 | time >= 4500) %>% 
    # mutate(Time_Window = ifelse(time <= 500, "Transient", "Equilibrium")) %>%
    # mutate(Time_Window = factor(Time_Window, levels = c("Transient", "Equilibrium"))) %>%
    pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"),
                 names_to = "species", values_to = "biomass") %>%
    ggplot(mapping = aes(x = time, y = biomass, color = species)) +
    labs(x = "Time", y = "Biomass") +
    geom_line(lwd = 1) +
    #facet_wrap(~ Time_Window, scales = "free_x") +
    scale_y_continuous(limits = c(0, 5))+
    scale_colour_manual(name = "Species",
                        labels = State_labels,
                        values = State_values)+
    #theme_minimal()+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.spacing = unit(0.2, "lines"))+
    transition_reveal(along = time)
}

Disturbance = TimeSeriesdisturbance(times1, times2, Initial, parms_Low_m1, parms_High_m1)

GIF <- animate(Disturbance, 
               bg = 'transparent',
               renderer = gifski_renderer(loop = FALSE), 
               width = 1600, height = 900)

#### Save GIF
anim_save(filename = ".gif", GIF)

#AFTER disturbance (m1 = 0.07): 1.083019 0.3119674 0.2303794 1.253033 0.4102449 1.257898


## H invasion----
H_invasion = function(times, state, parms){
  with(as.list(c(state, parms)), {
    H = H + 5.3
    P1H = P1H
    P2H = P2H
    P1 = P1
    P2 = P2
    S = S
    return(c(H, P1H, P2H, P1, P2, S))
  })
}

times = seq(0, 5000, by = 0.1)
#Initial = c(H = 1.083019, P1H = 0.3119674, P2H = 0.2303794, P1 = 1.253033, P2 = 0.4102449, S = 1.257898) #m1 = 0.07
Initial = c(H = 0.7383519, P1H = 0.2883162, P2H = 0.1024059, P1 = 1.67909, P2 = 0.2674836, S = 0.9609926) #m1 = 0.06

parms_EC_E2H <- c(r = 1, K = 10,
                  a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
                  b1 = 0.2, b2 = 0.42, m1 = 0.06, m2 = 0.01, e1H = 0.5, e2H = 0.5,
                  o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)

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
    #filter(time >= 2400 & time <= 2700) %>% 
    # mutate(Time_Window = ifelse(time <= 500, "Transient", "Equilibrium")) %>%
    # mutate(Time_Window = factor(Time_Window, levels = c("Transient", "Equilibrium"))) %>%
    pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"),
                 names_to = "species", values_to = "biomass") %>%
    ggplot(mapping = aes(x = time, y = biomass, color = species)) +
    labs(x = "Time", y = "Biomass") +
    geom_line(lwd = 1) +
    #facet_wrap(~ Time_Window, scales = "free_x") +
    scale_y_continuous(limits = c(0, 4.2))+
    scale_colour_manual(name = "Species",
                        labels = State_labels,
                        values = State_values)+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      panel.spacing = unit(0.2, "lines"))+
      #legend.position = "none")#+
    transition_reveal(along = time)
  #axis.title.y = element_blank(), 
  #axis.text.y = element_blank())
}

Invasion = 
  TimeSeriesInvasion(times, Initial, parms_EC_E2H, InvasionTimes = c(500))

GIF <- animate(Invasion, 
               bg = 'transparent',
               renderer = gifski_renderer(loop = FALSE), 
               width = 400, height = 300)

#### Save GIF
anim_save(filename = "m1 = 006 H = 5p3 invasion success.gif", GIF)

## Patch work-----
Disturbance + Invasion + 
  plot_layout(widths = c(0.2)) & 
  theme(plot.tag = element_text(size = 10, face = "bold"))

Disturbance + labs(tag = "(A)") + Invasion + labs(tag = "(B)") + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", plot.tag = element_text(size = 15, face = "bold"))

