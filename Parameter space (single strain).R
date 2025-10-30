library(deSolve)
library(tidyverse)
library(pracma) 

#### ODE setup----
M2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dH_dt = (h1 * o1 * P1H) - (c1 * b1 * P1) * H - (d * H)
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dS_dt = r * S * (1-S/K) - (a1 * P1 + psi1 * a1 * P1H) * S 
    return(list(c(dH_dt, dP1H_dt, dP1_dt, dS_dt)))
  })
}

### Model parameters ----

times <- seq(0, 10000, by = 0.1)
#state <- c(H = 1, P1H = 0, P2H = 0, P1 = 2, P2 = 2, S = 5)
state <- c(H = 0.01, P1H = 0, P1 = 0.01, S = 0.2)
#Change alpha and beta--
parms <- c(r = 1, K = 10,
           a1 = 0.5, psi1 = 1, e1 = 0.5,
           b1 = 0.35, m1 = 0.05, e1H = 0.5,
           o1 = 0.8, h1 = 1, c1 = 0.9, d = 0.03, DL = 0)

comp_out = expand.grid(b1 = seq(0, 1.5, by = 0.05))
#A = c(seq(0.05, 0.5, by = 0.05))
#B = c(rep(0.2, length(seq(0.05, 0.5, by = 0.05))))
#comp_out = matrix(c(A,B), ncol = 2, nrow = 10, byrow = F)

#### Create saving space for simulation output----
#combine different parameters and different state variables together
comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = length(state) + 2)))
names(comp_out) <- c("b1", "H", "P1H", "P1", "S", "sd", "Cycle")


#### Simulate the ODE across the parameter space----
start_time <- Sys.time()
for(i in 1:dim(comp_out)[1]){
  
  temp_parms <- parms
  temp_parms["b1"] <- comp_out[i, ]$b1
  
  temp_out <- ode(func = M2, 
                  times = times, 
                  y = state, 
                  parms = temp_parms)
  
  if (sd(temp_out[(nrow(temp_out)-round(length(times)*0.3)):nrow(temp_out),"S"]) > 1e-8){
    print(paste0("The system has not reached equilibrium, or there is a cycle when b1 = ", temp_parms["b1"]))
    
    N_semifinal <- nrow(na.omit(temp_out))
    
    temp_out <- ode(func = M2, 
                    times = times, 
                    y = temp_out[N_semifinal, -1],   #抓取最後一列的資料，排除時間
                    parms = temp_parms)
    
    N_final <- nrow(na.omit(temp_out))
    
    if (sd(temp_out[(nrow(temp_out)-round(length(times)*0.3)):nrow(temp_out),"S"]) > 1e-8){
      print(paste0("The system has not reached equilibrium after another ", times[length(times)], " steps. There should be a cycle when b1 = ", temp_parms["b1"]))
      
      k = 1
      avg = c(rep(NA, 4)) #create a vector to store the long term average of each species 
      for (j in names(comp_out)[2:5]) { #run a loop that can input each species
        peaks = findpeaks(temp_out[, j])[ ,2]
        periods = peaks[length(peaks)] - peaks[length(peaks) - 1]
        avg[k] = mean(temp_out[(length(times) - periods + 1):length(times), j]) #store the long term average in to the vector that just created
        k = k + 1 #increase k so that we can put the value of next species into the vector
      }
      k = 0 #reset k (it's unnecessary actually)
      comp_out[i, 2:5] = avg
      comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.3)):nrow(temp_out),"S"])
      comp_out[i, "Cycle"] = "T"
    }else{
      N_final <- nrow(na.omit(temp_out))
      comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.3)):nrow(temp_out),"S"])
      comp_out[i, 2:5] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
      comp_out[i, "Cycle"] = "F"
    }
  }else{
    N_final <- nrow(na.omit(temp_out))
    comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.3)):nrow(temp_out),"S"])
    comp_out[i, 2:5] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
    comp_out[i, "Cycle"] = "F"
  }
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
  select(c(b1, H, P1H, P1, S)) %>%
  mutate(Total = P1H + P1 + H + S) %>%
  pivot_longer(names_to = "Species", values_to = "Abundance", -c(b1)) %>%
  ggplot(aes(x = b1, y = Abundance, color = Species)) +
  geom_line(lwd = 1) + 
  labs(title = "α = 0.5", x = "β", y = "Abundance", color = "Species")+
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "S" = "Host", "H" = "Hyper", "Total" = "Total"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "S" = "#00AF66", "H" = "#C03728", "Total" = "blue"))


##Save the simulation result
saveRDS(comp_out, "SS2")
