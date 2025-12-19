library(deSolve)
library(tidyverse)

### Build up tritrophic model ----
Model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1 * H) - (d * H)
    dP1_dt = e1 * a1 * P1 * S - o1 * P1 * H - m1 * P1
    dS_dt = r * S * (1-S/K) - a1 * P1 * S
    return(list(c(dH_dt, dP1_dt, dS_dt)))
  })
}

### Type-II functional response
Model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = h1 * (o1 * P1) / (Ks + P1) * H - d * H
    dP1_dt = e1 * (a1 * S) / (Ks + S) * P1  - (o1 * P1) / (Ks + P1) * H - m1 * P1
    dS_dt = r * S * (1-S/K) - (a1 * S) / (Ks + S) * P1
    return(list(c(dH_dt, dP1_dt, dS_dt)))
  })
}

### Type-II functional response
Model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = h1 * (o1 * P1) / (Ks + P1) * H - d * H
    dP1_dt = e1 * (a1 * S)* P1  - (o1 * P1) / (Ks + P1) * H - m1 * P1
    dS_dt = r * S * (1-S/K) - (a1 * S) * P1
    return(list(c(dH_dt, dP1_dt, dS_dt)))
  })
}

### Model parameters ----

times <- seq(0, 5000, by = 0.05)
#state <- c(H = 1, P1H = 0, P1 = 2, S = 5)
state <- c(H = 0.01, P1 = 0.01, S = 0.2)
#Change alpha and beta--
parms <- c(r = 1, K = 10,
           a1 = 0.5, psi1 = 1, e1 = 0.5,
           b1 = 0.45, m1 = 0.2, e1H = 0.5,
           o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.03, Ks = 0.8)


### Model application ----
pop_size = ode(func = Model, times = times, y = state, parms = parms)

tail(pop_size)
#tail(mutate(as.data.frame(pop_size), P1total = P1 + P1H))
### Plot the result ----
## Plotting
pop_size %>%
  as.data.frame() %>%
  filter(time %% 1 == 0) %>%
  #filter(time > 1000, time < 5000) %>%
  #View()
  pivot_longer(cols = c("H", "P1", "S"), #"H", "P1H", "P1", "S" 
               names_to = "species", values_to = "biomass") %>%
  #filter(species == c("H","S")) %>%
  ggplot(mapping = aes(x = time, y = biomass, color = species)) +
  labs(x = "Time", y = "Biomass", title = paste0("α = ", parms["a1"], ", β = ", parms["b1"])) +  #paste0("r =", parms["r"]))β
  geom_line(lwd = 1) +
  scale_colour_manual(labels = c("H" = "Hyper", "P1" = expression(P[1]), 
                                 "P1H" = expression(P[1/H]), "S" = "Host"),
                      values = c("H" = "#C03728", "P1" = "#BCAAA4", 
                                 "P1H" = "#82491E", "S" = "#00AF66"))

k = 0
avg = rep(0, 3)
for (j in colnames(pop_size)[2:4]) { #run a loop that can input each species
  peaks = findpeaks(pop_size[, j])[ ,2]
  periods = peaks[length(peaks)] - peaks[length(peaks) - 1]
  avg[k] = mean(pop_size[(length(times) - periods + 1):length(times), j]) #store the long term average in to the vector that just created
  k = k + 1 #increase k so that we can put the value of next species into the vector
}
avg

#Analytic Sol:

with(as.list(parms), {
  H = (e1*a1*K*(1-(a1*d)/(h1*o1*r))-m1)/o1
  P = d/(h1*o1)
  S = K*(1-(a1*d)/(h1*o1*r))
  return(c(H = H, P = P, S = S))
})
 