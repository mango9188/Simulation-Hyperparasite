library(deSolve)
library(tidyverse)

### Build up Mechanism 2 ----
M2_S <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1H) - (c1 * b1 * P1) * H - (d * H)
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dS_dt = r * S * (1-S/K) - (a1 * P1 + psi1 * a1 * P1H) * S 
    return(list(c(dH_dt, dP1H_dt, dP1_dt, dS_dt)))
  })
}

### Model parameters ----

times <- seq(0, 10000, by = 1)
#state <- c(H = 1, P1H = 0, P2H = 0, P1 = 2, P2 = 2, S = 5)
state <- c(H = 0.01, P1H = 0, P1 = 0.01, S = 0.2)
#Change alpha and beta--
parms <- c(r = 1, K = 10,
           a1 = 0.35, psi1 = 1, e1 = 0.5,
           b1 = 0.3, m1 = 0.06, e1H = 0.5,
           o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.03, DL = 0)


### Model application ----
pop_size = ode(func = M2_S, times = times, y = state, parms = parms)

tail(pop_size)
tail(mutate(as.data.frame(pop_size), P1total = P1 + P1H))
### Plot the result ----
## Plotting
pop_size %>%
  as.data.frame() %>%
  filter(time %% 1 == 0) %>%
  #filter(time > 1000, time < 5000) %>%
  #View()
  pivot_longer(cols = c("H", "P1H", "P1", "S"), #"H", "P1H", "P1", "S" 
               names_to = "species", values_to = "biomass") %>%
  #filter(species == c("H","S")) %>%
  ggplot(mapping = aes(x = time, y = biomass, color = species)) +
  labs(x = "Time", y = "Biomass", title = paste0("m = ", parms["m1"])) +  #paste0("r =", parms["r"]))β
  geom_line(lwd = 1) +
  scale_colour_manual(labels = c("H" = "Hyper", "P1" = expression(P[1]), 
                                 "P1H" = expression(P[1/H]), "S" = "Host"),
                      values = c("H" = "#C03728", "P1" = "#BCAAA4", 
                                 "P1H" = "#82491E", "S" = "#00AF66"))

## Plotting the per capita growth rate
pop_size %>%
  as.data.frame() %>%
  filter(time %% 1 == 0) %>%
  filter(time > 0, time < 500) %>%
  #View()
  pivot_longer(cols = c("per_H", "per_P1H", "per_P1", "per_S"),
               #"per_H", "per_P1H", "per_P1", "per_S"
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
sd(pop_size[(nrow(pop_size)-round(length(times)*0.2)):nrow(pop_size),"S"]) > 1e-8

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
