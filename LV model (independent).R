library(deSolve)
library(tidyverse)

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
times <- seq(0, 40000, by = 0.1)
#state <- c(H = 1, P1H = 0, P2H = 0, P1 = 2, P2 = 2, S = 5)
state <- c(H = 0.1, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 2)
#state <- c(H = 0.1638255, P1H = 0.02361139, P2H =  0.1168105, P1 =  0.6125324, P2 = 1.346811, S = 0.4553895)
#Change alpha and beta--
parms <- c(epsilon = 1,
           r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

#parms <- c(r = 1, K = 10,a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.0536, e1H = 0.5, e2H = 0.5,o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

### Model application ----
pop_size = ode(func = M2, times = times, y = state, parms = parms)

tail(pop_size)

### Plot the result ----
## Plotting
pop_size %>%
  as.data.frame() %>%
  filter(time %% 10 == 0) %>%
  #filter(time > 1000) %>%
  #View()
  pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"), #"H", "P1H", "P2H", "P1", "P2", "S" 
             names_to = "species", values_to = "biomass") %>%
  #filter(species == c("H","S")) %>%
  ggplot(mapping = aes(x = time, y = biomass, color = species)) +
  labs(x = "Time", y = "Biomass") +  #title = expression(P[1]~"win") paste0("r =", parms["r"])) expression(α[1] == 0.35)+ #title = expression(α[1] == 0.35 ~","~ β[1] == 0.2 )
  geom_line(lwd = 1) +
  #geom_hline(yintercept = 0.28905636, color = "black", linetype = "dashed", size = 1) +
  scale_colour_manual(labels = c("H" = "Hyper", "P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "S" = "Host"),
                      values = c("H" = "#C03728", "P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", "S" = "#00AF66"))

#ggsave("P2 in m2 00536 Tless500.png", width = 20, height = 11.25, units = "cm", dpi = 1600)

### Plotting the per capita growth rate----
pop_size %>%
  as.data.frame() %>%
  filter(time %% 1 == 0) %>%
  #filter(time > 100) %>%
  #View()
  pivot_longer(cols = c("per_H", "per_P1H", "per_P2H", "per_P1", "per_P2", "per_S"), 
               #"H", "P1H", "P2H", "P1", "P2", "S" 
               #"per_H", "per_P1H", "per_P2H", "per_P1", "per_P2", "per_S"
               names_to = "species", values_to = "rate") %>%
  #filter(species == c("H","S")) %>%
  ggplot(mapping = aes(x = time, y = rate, color = species)) +
  labs(x = "Time", y = "Rate", title = expression(m[1] < m[2])) +  #paste0("r =", parms["r"]))
  geom_line(lwd = 1) +
  scale_colour_manual(labels = c("per_H" = expression(r[H]), "per_P1H" = expression(r[P[1/H]]),
                                 "per_P2H" = expression(r[P[2/H]]) , "per_P1" = expression(r[P[1]]),
                                 "per_P2" = expression(r[P[2]]), "per_S" = expression(r[S])),
                      values = c('per_H' = "#C03728", 'per_P1H' = "#82491E", 
                                 'per_P2H' = "#546E7A", 'per_P1' = "#BCAAA4", 
                                 'per_P2' = "#B0BEC5", 'per_S' = "#00AF66"))

View(pop_size)
sd(pop_size[(nrow(pop_size)-round(length(times)*0.35)):nrow(pop_size),"S"]) > 1e-8

## Theme setting
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


### Calculate the S* of each time step
pop_size =
  pop_size %>%
  as.data.frame() %>%
  mutate(
       S2 = with(as.list(parms), {
         ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)}),
       S1 = with(as.list(parms), {
         ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)})
       )

## Plotting
pop_size %>%
  as.data.frame() %>%
  filter(time %% 1 == 0) %>%
  filter(time < 5000) %>%
  #View()
  pivot_longer(cols = c("H", "P1", "P2", "S1", "S2" ), #"H", "P1H", "P2H", "P1", "P2", "S" 
               names_to = "species", values_to = "biomass") %>%
  #filter(species == c("H","S")) %>%
  ggplot(mapping = aes(x = time, y = biomass, color = species)) +
  labs(x = "Time", y = "Biomass") +  #, title = expression(P[1]~"win") paste0("r =", parms["r"])) expression(α[1] == 0.35)
  geom_line(lwd = 1) +
  #geom_hline(yintercept = 0.133791930, color = "black", linetype = "dashed", size = 1) +
  scale_colour_manual(labels = c("H" = "Hyper", "P1" = expression(P[1]), "P1H" = expression(P[1/H]), "P2" = expression(P[2]), "P2H" = expression(P[2/H]), "S" = "Host"),
                      values = c("H" = "#C03728", "P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", "S" = "#00AF66",
                                 "S1" = "#5e4a59", "S2" = "#a45ea1"))

tail(pop_size)
ggsave("m2 00536 sStar.png", width = 20, height = 11.25, units = "cm", dpi = 800)
