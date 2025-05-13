library(deSolve)
library(tidyverse)
library(paletteer)
library(pracma) 

#### ODE setup----
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

times <- seq(0, 80000, by = 0.01)
state <- c(H = 0.1, P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.2)
parms <- c(r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

comp_out = expand.grid(m1 = seq(0.01, 0.1, by = 0.01), m2 = seq(0.01, 0.1, by = 0.01))

#### Create saving space for simulation output----
#combine different parameters and different state variables together
comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = length(state) + 2))) #2->include sd and cycle
names(comp_out) <- c("m1", "m2", "H", "P1H", "P2H", "P1", "P2", "S", "sd", "Cycle")


#### Simulate the ODE across the parameter space----
start_time <- Sys.time()
for(i in 1:dim(comp_out)[1]){
  
  temp_parms <- parms
  temp_parms["m1"] <- comp_out[i, ]$m1
  temp_parms["m2"] <- comp_out[i, ]$m2
  
  temp_out <- ode(func = M2, 
                  times = times, 
                  y = state, 
                  parms = temp_parms)
  
  if (sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"]) > 1e-8){
    print(paste0("The system has not reached equilibrium, or there is a cycle when m1 = ", temp_parms["m1"], ", m2 = ", temp_parms["m2"]))
    
    N_semifinal <- nrow(na.omit(temp_out))
    
    temp_out <- ode(func = M2, 
                    times = times, 
                    y = temp_out[N_semifinal, -1],   #抓取最後一列的資料，排除時間
                    parms = temp_parms)
    
    N_final <- nrow(na.omit(temp_out))
    
    if (sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"]) > 1e-8){
      print(paste0("The system has not reached equilibrium after another ", times[length(times)], " steps. There should be a cycle when m1 = ", temp_parms["m1"], ", m2 = ", temp_parms["m2"]))
      
      k = 1
      avg = c(rep(NA, 6)) #create a vector to store the long term average of each species 
      for (j in names(comp_out)[3:8]) { #run a loop that can input each species
        peaks = findpeaks(temp_out[, j])[ ,2]
        
        if (length(peaks) >= 3){
          periods = peaks[length(peaks)] - peaks[length(peaks) - 2]
          avg[k] = mean(temp_out[(length(times) - periods + 1):length(times), j]) #store the long term average in to the vector that just created
        }else{
          avg[k] = mean(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out), j])
        }
        
        k = k + 1 #increase k so that we can put the value of next species into the vector
      }
      k = 0 #reset k (it's unnecessary actually)
      comp_out[i, 3:8] = avg
      comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"])
      comp_out[i, "Cycle"] = "T"
    }else{
      N_final <- nrow(na.omit(temp_out))
      comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"])
      comp_out[i, 3:8] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
      comp_out[i, "Cycle"] = "F"
    }
  }else{
    N_final <- nrow(na.omit(temp_out))
    comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"])
    comp_out[i, 3:8] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
    comp_out[i, "Cycle"] = "F"
  }
}

end_time <- Sys.time()
end_time - start_time

###Data analysis----
extinct_thres = 1e-7
comp_out$Outcome =
  ifelse(comp_out[, "P1H"] + comp_out[, "P1"] < extinct_thres, "P2 win", 
         ifelse(comp_out[, "P2H"] + comp_out[, "P2"] < extinct_thres, "P1 win", "Coexist"))

comp_out$Outcome2 =
  apply(comp_out[, c("H", "P1H", "P2H", "P1", "P2")], 1, function(row){
    paste(ifelse(row > extinct_thres, "T", "F"), collapse = "")
  })

