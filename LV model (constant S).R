library(deSolve)
library(tidyverse)

# Here, I set S as a constant.

### Build up Mechanism 2 ----
Scon_M2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
    return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt)))
  })
}

#Model parameters
times <- seq(0, 10000, by = 0.01)
state <- c(H = 0.1, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01)
parms <- c(S = 0.328,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

### Model application ----
pop_size = ode(func = Scon_M2, times = times, y = state, parms = parms)

tail(pop_size)

### Plot the result ----
## Plotting
pop_size %>%
  as.data.frame() %>%
  filter(time %% 10 == 0) %>%
  #filter(time < 1500) %>%
  #View()
  pivot_longer(cols = c("P1H", "P2H", "P1", "P2", "H"), #"H", "P1H", "P2H", "P1", "P2", "S" 
               names_to = "species", values_to = "biomass") %>%
  #filter(species == c("H","S")) %>%
  ggplot(mapping = aes(x = time, y = biomass, color = species)) +
  labs(x = "Time", y = "Biomass", title = "Small H") +  #paste0("r =", parms["r"])) expression(α[1] == 0.35)
  geom_line(lwd = 1) +
  scale_colour_manual(labels = c("H" = "Hyper", "P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "S" = "Host"),
                      values = c("H" = "#C03728", "P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", "S" = "#00AF66"))


View(pop_size)
sd(pop_size[(nrow(pop_size)-round(length(times)*0.35)):nrow(pop_size),"S"]) > 1e-8

### Theme setting ----
A = 
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 10),
    legend.text.align = 0)
theme_set(A)


### Analytic way ----
parms <- c(H = 0.002,
           r = 1, K = 10,
           a1 = 0.45, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.15, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

S1 = with((as.list(parms)), {
  S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
  return(S1)
})

S2 = with((as.list(parms)), {
  S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
  return(S2)
})

if(S1 > S2){
  print("P2 will win the competition.")
  P2 = with((as.list(parms)), {
    P2 = r * (1 - (S2/K)) * (m2 + o2) / (a2 * (m2 + o2 + psi2 * b2 * H))
    return(P2)
  })
  P2H = with((as.list(parms)), {
    P2H = P2 * (b2 * H) / (m2 + o2)
    return(P2H)
  })
}else if (S1 < S2){
  print("P1 will win the competition.")
  P1 = with((as.list(parms)), {
    P1 = r * (1 - (S1/K)) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
    return(P1)
  })
  P1H = with((as.list(parms)), {
    P1H = P1 * (b1 * H) / (m1 + o1)
    return(P1H)
  })
}else{
  print("The system might coexist.")
}

