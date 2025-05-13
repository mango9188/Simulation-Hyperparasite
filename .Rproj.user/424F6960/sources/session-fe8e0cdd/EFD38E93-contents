#### Create a parameter space for low H ----
library(deSolve)
library(tidyverse)

### Build up Mechanism 2 ----
Hcon_M2 <- function(times, state, parms) {
  with(as.list(c(state, parms)), {
    dP1H_dt = (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    dP2H_dt = (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
    dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
    dS_dt = r * S * (1-S/K) - (a1 * P1 + a2 * P2 + psi1 * a1 * P1H + psi2 * a2 * P2H) * S 
    return(list(c(dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
  })
}

### Model parameters ----
times <- seq(0, 80000, by = 0.01)
state <- c(P1H = 0, P2H = 0, P1 = 0.01, P2 = 0.01, S = 0.2)
parms <- c(H = 0.0239533,
           r = 1, K = 10,
           a1 = 0.5, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

###Create a parameter space and calculate the S* for each grids in parameter space(No PiH)----

comp_out = expand.grid(a1 = seq(0.05, 1, by = 0.05), b1 = seq(0.05, 1, by = 0.05))

#### Create saving space for output
#combine different parameters and different state variables together
comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = 2))) #Number = number of state variable
names(comp_out) <- c("a1", "b1", "S1", "S2")
#Because now H is small, it is assumed that all the competition outcome will follow exploitative competition outcome.
###This is for prediction, no simulation.
comp_out[, "S2"] = 
  with(as.list(parms), {
    S2 = m2/(e2*a2)
    return(S2)
    }) #S* of strain 2 is fixed now.


start_time <- Sys.time()
for(i in 1:dim(comp_out)[1]){
  
  temp_parms <- parms
  temp_parms["a1"] <- comp_out[i, ]$a1
  #temp_parms["b1"] <- comp_out[i, ]$b1 #Calculating S* will not need the value of b1.
  
  comp_out[i, "S1"] = 
    with(as.list(temp_parms), {
      S1 = m1/(e1*a1)
    })
  #calculate the S* for each grids in parameter space S* of strain 1
}
end_time <- Sys.time()
end_time - start_time

comp_out$Outcome <- 
  ifelse(comp_out[, 3] < comp_out[, 4], "P1+H", 
         ifelse(comp_out[, 3] > comp_out[, 4], "P2+H", "Coexist"))

comp_out

unique_outcomes = unique(comp_out$Outcome)
mycolor = c("Coexist" = "#BA6338FF",#AC
            "P1+H" = "#466983FF",#P1+H
            "P2+H" = "#749B58FF"#P2+H
            )

outcome_labels <- c(
  "Coexist" = "Coexistence",
  "P1+H" = expression(P[1] + P[1/H] ~ "win"),
  "P2+H" = expression(P[2] + P[2/H] ~ "win"))

ggplot(comp_out, aes(x = a1, y = b1, z = Outcome, fill = Outcome)) +
  geom_tile() +
  labs(title = expression(m[1] < m[2] ~","~ S^{"*"} ~ "rule"), x = expression(α[1]), y = expression(β[1]))+ #title = "δ = 0 (No vertical transmission)",  ~","~ r == 1.5
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = mycolor, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1)

###Create a parameter space with analytic way----

comp_out = expand.grid(a1 = seq(0.05, 1, by = 0.05), b1 = seq(0.05, 1, by = 0.05))


comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = 6))) #Number = number of state variable
names(comp_out) <- c("a1", "b1", "P1H", "P2H", "P1", "P2", "S1", "S2")

comp_out[, "S2"] = with((as.list(parms)), {
  S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
  return(S2)
})##The S* of strain 2 is fixed.


comp_out[, "S1"] = with((as.list(c(comp_out, parms))), {
  S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
  return(S1)
})##The S* of strain 1 depends on the value of parameter expanding.


### Or using for-loop (It will be a bit slower and messy in this way.)
for(i in 1:dim(comp_out)[1]){
  temp_parms <- parms
  temp_parms["a1"] <- comp_out[i, ]$a1
  temp_parms["b1"] <- comp_out[i, ]$b1
  comp_out[i, 7] = 
    with(as.list(temp_parms), {
      S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
      return(S1)
    })
  #calculate the S* for each grids in parameter space S* of strain 1
}



for(i in 1:dim(comp_out)[1]){
  temp_parms = parms
  temp_parms["a1"] = comp_out[i, "a1"]
  temp_parms["b1"] = comp_out[i, "b1"]
  
  if(comp_out[i, "S1"] > comp_out[i, "S2"]){
    comp_out[i, "P2"] = with((as.list(c(comp_out[i,], temp_parms))), {
      P2 = r * (1 - (S2/K)) * (m2 + o2) / (a2 * (m2 + o2 + psi2 * b2 * H))
      return(P2)
    })
    comp_out[i, "P2H"] = with((as.list(c(comp_out[i,], temp_parms))), {
      P2H = P2 * (b2 * H) / (m2 + o2)
      return(P2H)
    })
  }else{
    comp_out[i, "P1"] = with((as.list(c(comp_out[i,], temp_parms))), {
      P1 = r * (1 - (S1/K)) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
      return(P1)
    })
    comp_out[i, "P1H"] = with((as.list(c(comp_out[i,], temp_parms))), {
      P1H = P1 * (b1 * H) / (m1 + o1)
      return(P1H)
    })
  }
}


comp_out$Outcome <- 
  ifelse(comp_out[, "S1"] < comp_out[, "S2"], "P1+H", 
         ifelse(comp_out[, "S1"] > comp_out[, "S2"], "P2+H", "Coexist"))

comp_out

unique_outcomes = unique(comp_out$Outcome)
mycolor = c("Coexist" = "#BA6338FF",#AC
            "P1+H" = "#466983FF",#P1+H
            "P2+H" = "#749B58FF"#P2+H
)

outcome_labels <- c(
  "Coexist" = "Coexistence",
  "P1+H" = expression(P[1] + P[1/H] ~ "win"),
  "P2+H" = expression(P[2] + P[2/H] ~ "win"))

ggplot(comp_out, aes(x = a1, y = b1, z = Outcome, fill = Outcome)) +
  geom_tile() +
  labs(title = expression(m[1] < m[2] ~","~ S^{"*"} ~ "rule"), x = expression(α[1]), y = expression(β[1]))+ #title = "δ = 0 (No vertical transmission)",  ~","~ r == 1.5
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = mycolor, labels = outcome_labels) +
  #scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+#automatically choose color
  #scale_fill_manual(values = setNames(paletteer_d("ggsci::default_igv")[1:length(all_comb)], all_comb))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1)


###Create saving space for simulation output----
comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = length(state) + 2))) #2->include sd and cycle
names(comp_out) <- c("a1", "b1", "P1H", "P2H", "P1", "P2", "S", "sd", "Cycle")

###Simulate the ODE across the parameter space----
start_time <- Sys.time()
for(i in 1:dim(comp_out)[1]){
  
  temp_parms <- parms
  temp_parms["a1"] <- comp_out[i, ]$a1
  temp_parms["b1"] <- comp_out[i, ]$b1
  
  temp_out <- ode(func = Hcon_M2, 
                  times = times, 
                  y = state, 
                  parms = temp_parms)
  
  if (sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"]) > 1e-8){
    print(paste0("The system has not reached equilibrium, or there is a cycle when a1 = ", temp_parms["a1"], ", b1 = ", temp_parms["b1"]))
    
    N_semifinal <- nrow(na.omit(temp_out))
    
    temp_out <- ode(func = M2, 
                    times = times, 
                    y = temp_out[N_semifinal, -1],   #抓取最後一列的資料，排除時間
                    parms = temp_parms)
    
    N_final <- nrow(na.omit(temp_out))
    
    if (sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"]) > 1e-8){
      print(paste0("The system has not reached equilibrium after another ", times[length(times)], " steps. There should be a cycle when a1 = ", temp_parms["a1"], ", b1 = ", temp_parms["b1"]))
      
      k = 1
      avg = c(rep(NA, 5)) #create a vector to store the long term average of each species 
      for (j in names(comp_out)[3:7]) { #run a loop that can input each species
        peaks = findpeaks(temp_out[, j])[ ,2]
        periods = peaks[length(peaks)] - peaks[length(peaks) - 1]
        avg[k] = mean(temp_out[(length(times) - periods + 1):length(times), j]) #store the long term average in to the vector that just created
        k = k + 1 #increase k so that we can put the value of next species into the vector
      }
      k = 0 #reset k (it's unnecessary actually)
      comp_out[i, 3:7] = avg
      comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"])
      comp_out[i, "Cycle"] = "T"
    }else{
      N_final <- nrow(na.omit(temp_out))
      comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"])
      comp_out[i, 3:7] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
      comp_out[i, "Cycle"] = "F"
    }
  }else{
    N_final <- nrow(na.omit(temp_out))
    comp_out[i, "sd"] = sd(temp_out[(nrow(temp_out)-round(length(times)*0.2)):nrow(temp_out),"S"])
    comp_out[i, 3:7] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
    comp_out[i, "Cycle"] = "F"
  }
}

end_time <- Sys.time()
end_time - start_time

##Data analysis----
extinct_thres = 1e-7
comp_out$Outcome =
  ifelse(comp_out[, "P1H"] + comp_out[, "P1"] < extinct_thres, "P2 win", 
         ifelse(comp_out[, "P2H"] + comp_out[, "P2"] < extinct_thres, "P1 win", "Coexist"))

comp_out$Outcome2 =
  apply(comp_out[, c("H", "P1H", "P2H", "P1", "P2")], 1, function(row){
    paste(ifelse(row > extinct_thres, "T", "F"), collapse = "")
  })

