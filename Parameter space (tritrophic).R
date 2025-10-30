library(deSolve)
library(tidyverse)
library(pracma) 

#### ODE setup----
Model <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1 * H) - (d * H)
    dP1_dt = e1 * a1 * P1 * S - o1 * P1 * H - m1 * P1
    dS_dt = r * S * (1-S/K) - a1 * P1 * S 
    return(list(c(dH_dt, dP1_dt, dS_dt)))
  })
}

### Model parameters ----

times <- seq(0, 6000, by = 0.1)
#state <- c(H = 1, P1H = 0, P2H = 0, P1 = 2, P2 = 2, S = 5)
state <- c(H = 0.1, P1 = 0.1, S = 0.2)
#Change alpha and beta--
parms <- c(r = 1, K = 10,
           a1 = 0.5, psi1 = 1, e1 = 0.5,
           b1 = 0.45, m1 = 0.8, e1H = 0.5,
           o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.03, DL = 0)

comp_out = expand.grid(m1 = seq(1, 2.5, by = 0.05))

#### Create saving space for simulation output----
#combine different parameters and different state variables together
comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = length(state))))
names(comp_out) <- c("m1", "H", "P1", "S")


#### Simulate the ODE across the parameter space----
start_time <- Sys.time()
for(i in 1:dim(comp_out)[1]){
  
  temp_parms = parms
  temp_parms["m1"] <- comp_out[i, ]$m1
  
  temp_out <- ode(func = Model, 
                  times = times, 
                  y = state, 
                  parms = temp_parms)

  N_final <- nrow(temp_out)
  comp_out[i, 2:4] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
}
end_time <- Sys.time()
end_time - start_time

###Data analysis----
extinct_thres = 1e-7
comp_out$Outcome =
  ifelse(comp_out[, "H"] < extinct_thres, "H extincts", "Coexist")


####Plot the result----

###bifurcation plot
comp_out %>%
  select(c(m1, P1)) %>% #m1, H, P1, S
  #mutate(Total = P1H + P1 + H + S) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(m1)) %>%
  ggplot(aes(x = m1, y = Abundance, color = Species)) +
  geom_line(lwd = 1) + 
  labs(title = "Tritrophic model", x = "Mortality (m)", y = "Abundance", color = "Species")+
  scale_colour_manual(labels = 
                        c("P1" = "Pathogen", "P1H" = expression(P[1/H]),
                          "S" = "Host", "H" = "Hyperparasite", "Total" = "Total"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "S" = "#00AF66", "H" = "#C03728", "Total" = "blue"))


ggsave("Tritrophic model.png", width = 15, height = 11, units = "cm", dpi = 800)
##Save the simulation result
#saveRDS(comp_out, "SS2")
