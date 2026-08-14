library(tidyverse)
library(deSolve)

M2_LVmodel = function(max_time = 1000, H_0 = 1, P1H_0 = 0.01, P2H_0 = 0.01, P1_0 = 1, P2_0 = 1, S_0 = 10, r = 0.8, K = 10, a1 = 0.5, a2 = 0.2, R1 = 1, R2 = 1, e1 = 0.3, e2 = 0.3, b1 = 0.5, b2 = 0.1, m1 = 0.1, m2 = 0.1, e1H = 0.3, e2H = 0.3, o1 = 0.1, o2 = 0.01, h1 = 2, h2 = 0.4, c1 = 0.01, c2 = 0.01, d = 0.3, DL = 1){

### Model specification
  M2 <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dH_dt = (h1 * o1 * P1H) + (h2 * o2 * P2H) - (c1 * b1 * P1 + c2 * b2 * P2) * H - (d * H)
      dP1H_dt = (b1 * P1 * H) + DL * (e1H * R1 * a1 * P1H * S) - (o1 + m1) * P1H
      dP2H_dt = (b2 * P2 * H) + DL * (e2H * R2 * a2 * P2H * S) - (o2 + m2) * P2H
      dP1_dt = (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * R1 * a1 * P1H * S) - m1 * P1
      dP2_dt = (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * R2 * a2 * P2H * S) - m2 * P2
      dS_dt = r * S * (1-S/K) - (a1 * P1 + a2 * P2 + R1 * a1 * P1H + R2 * a2 * P2H) * S 
      return(list(c(dH_dt, dP1H_dt, dP2H_dt, dP1_dt, dP2_dt, dS_dt)))
    })
  }
  
  ### Model parameters
  times <- seq(0, max_time, by = 0.1)
  state <- c(H = H_0, P1H = P1H_0, P2H = P2H_0, P1 = P1_0, P2 = P2_0, S = S_0)
  parms <- c(r = r, K = K, a1 = a1, a2 = a2, R1 = R1, R2 = R2, e1 = e1, e2 = e2, b1 = b1, b2 = b2, m1 = m1, m2 = m2, e1H = e1H, e2H = e2H, o1 = o1, o2 = o2, h1 = h1, h2 = h2, c1 = c1, c2 = c2, d = d, DL = DL)
  
  ### Model application
  pop_size = ode(func = M2, times = times, y = state, parms = parms)

  View(pop_size)
  ### Visualize the population dynamics
  #pop_size %>%
    #as.data.frame() %>%
    #pivot_longer(cols = c("H", "P1H", "P2H", "P1", "P2", "S"),
    #             names_to = "species", values_to = "biomass") %>%
    #ggplot(mapping = aes(x = time, y = biomass, color = species))+
    #labs(x = "Time", y = "Biomass")+
    #geom_line(lwd = 1)+
    #scale_colour_manual(labels = c("H", "P1", "P1/H", "P2", "P2/H", "S"),
                        #values = c("#96281BFF", "#019875FF", "#FF847CFF",
                                   #"#99B898FF", "#E84A5FFF", "#FECEA8FF"))+
    #theme_bw()+
    #theme(
      #axis.text.x = element_text(size = 15),
      #axis.text.y = element_text(size = 15),
      #axis.title.x = element_text(size = 20),
      #axis.title.y = element_text(size = 20))
}

####P1P2 same, coexistence, DL = 1 (Vertical transmission, virus-fungi)####
M2_LVmodel(max_time = 100000, H_0 = 4, P1H_0 = 1, P2H_0 = 1, P1_0 = 1, P2_0 = 1, S_0 = 10,
           r = 0.81, K = 20,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.691, o2 = 0.691, h1 = 1, h2 = 1, c1 = 1, c2 = 1, d = 0.2, DL = 1)

M2_LVmodel(max_time = 100000, H_0 = 15, P1H_0 = 0, P2H_0 = 0, P1_0 = 6, P2_0 = 6, S_0 = 15,
           r = 0.81, K = 20,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.691, o2 = 0.691, h1 = 1, h2 = 1, c1 = 1, c2 = 1, d = 0.2, DL = 1)

##increasing of o1o2 will make the densities of P1 and P2 further increase, but reduce S and H.
M2_LVmodel(max_time = 10000, H_0 = 15, P1H_0 = 0, P2H_0 = 0, P1_0 = 6, P2_0 = 6, S_0 = 15,
           r = 0.81, K = 20,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.7, o2 = 0.7, h1 = 1, h2 = 1, c1 = 1, c2 = 1, d = 0.2, DL = 1)

##increasing of a1a2 will decrease P1 and P2.
M2_LVmodel(max_time = 10000, H_0 = 15, P1H_0 = 0, P2H_0 = 0, P1_0 = 6, P2_0 = 6, S_0 = 15,
           r = 0.81, K = 20,
           a1 = 0.105, a2 = 0.105, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.7, o2 = 0.7, h1 = 1, h2 = 1, c1 = 1, c2 = 1, d = 0.2, DL = 1)

##increasing of h1h2
M2_LVmodel(max_time = 10000, H_0 = 3, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 0.8, K = 15,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.7, o2 = 0.7, h1 = 1.5, h2 = 1.5, c1 = 1, c2 = 1, d = 0.2, DL = 1)

##increasing of h1h2, decrease the ini of H, the density of whole population won't change
M2_LVmodel(max_time = 10000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 0.8, K = 15,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.7, o2 = 0.7, h1 = 1.5, h2 = 1.5, c1 = 1, c2 = 1, d = 0.2, DL = 1)

##decreasing of K will make P increase, P/H, H and S decrease.
#This one possess reasonable parameters 
M2_LVmodel(max_time = 10000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 0.8, K = 10,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.7, o2 = 0.7, h1 = 1.5, h2 = 1.5, c1 = 1, c2 = 1, d = 0.2, DL = 1)


####P1P2 different, coexistence, DL = 1 (Vertical transmission, virus-fungi)####
M2_LVmodel(max_time = 10000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 0.8, K = 10,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.7, o2 = 0.7, h1 = 1.5, h2 = 1.5, c1 = 1, c2 = 1, d = 0.2, DL = 1)



####, DL = 1 (Vertical transmission, virus-fungi)####
M2_LVmodel(max_time = 1000, H_0 = 4, P1H_0 = 1, P2H_0 = 1, P1_0 = 1, P2_0 = 1, S_0 = 10,
           r = 0.81, K = 20,
           a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.08, b2 = 0.08, m1 = 0.02, m2 = 0.02, e1H = 0.5, e2H = 0.5,
           o1 = 0.691, o2 = 0.691, h1 = 1, h2 = 1, c1 = 1, c2 = 1, d = 0.2, DL = 1)

start_time <- Sys.time()
#----------------------------------------
####P1P2 same, coexistence, DL = 0 (No vertical transmission, bacteria-nematodes)####
M2_LVmodel(max_time = 100000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 1, K = 10,##Host
           a1 = 0.5, a2 = 0.5, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

end_time <- Sys.time()
end_time - start_time
####P1P2 different, coexistence, Delta = 0 (No vertical transmission, bacteria-nematodes)####
##by changing omega
M2_LVmodel(max_time = 100000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 1, K = 10,
           a1 = 0.5, a2 = 0.25, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.01, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

##by changing beta
M2_LVmodel(max_time = 100000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 1, K = 10,
           a1 = 0.5, a2 = 0.25, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.05, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

##by changing beta and omega
M2_LVmodel(max_time = 100000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 1, K = 10,
           a1 = 0.5, a2 = 0.25, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.05, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.1, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

#-----------------TEST#########----------------------------------

M2_LVmodel(max_time = 10000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 1, K = 10,
           a1 = 0.95, a2 = 0.5, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.65, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)
#-----------------TEST#########----------------------------------


####P1P2 different, P1 P1H win, DL = 0 (No vertical transmission, bacteria-nematodes)####
##by changing omega
M2_LVmodel(max_time = 10000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 1, K = 10,
           a1 = 0.5, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.2, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

##by changing beta
M2_LVmodel(max_time = 10000, H_0 = 1, P1H_0 = 0, P2H_0 = 0, P1_0 = 2, P2_0 = 2, S_0 = 5,
           r = 1, K = 10,
           a1 = 0.5, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.2, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)


####P1 win, DL = 0 (No vertical transmission, bacteria-nematodes)####
M2_LVmodel(H_0 = 1, P1H_0 = 0.01, P2H_0 = 0.01, P1_0 = 1, P2_0 = 1, S_0 = 10,
           r = 0.8, K = 10, a1 = 0.5, a2 = 0.2, R1 = 1, R2 = 1, e1 = 0.3, e2 = 0.3, b1 = 0.5, b2 = 0.1, m1 = 0.1, m2 = 0.1, e1H = 0.3, e2H = 0.3, o1 = 0.1, o2 = 0.01, h1 = 2, h2 = 0.4, c1 = 0.01, c2 = 0.01, d = 0.3, DL = 0)

#### DL = 0 (No vertical transmission, bacteria-nematodes)####
M2_LVmodel(max_time = 400, H_0 = 1, P1H_0 = 0.1, P2H_0 = 0.1, P1_0 = 1, P2_0 = 1, S_0 = 10,
           r = 0.8, K = 20, a1 = 0.1, a2 = 0.1, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5, b1 = 0.1, b2 = 0.1, m1 = 0.08, m2 = 0.08, e1H = 0.8, e2H = 0.8, o1 = 0.1, o2 = 0.1, h1 = 0.1, h2 = 0.1, c1 = 1, c2 = 1, d = 0.3, DL = 0)

####Pwin#####
M2_LVmodel(H_0 = 1, P1H_0 = 0.01, P2H_0 = 0.01, P1_0 = 1, P2_0 = 1, S_0 = 10,
           r = 0.8, K = 50, a1 = 0.5, a2 = 0.2, R1 = 0.001, R2 = 1, e1 = 0.3, e2 = 0.3, b1 = 0.5, b2 = 0.1, m1 = 0.1, m2 = 0.1, e1H = 0.3, e2H = 0.3, o1 = 0.1, o2 = 0.01, h1 = 2, h2 = 0.4, c1 = 0.01, c2 = 0.01, d = 0.3, DL = 1)