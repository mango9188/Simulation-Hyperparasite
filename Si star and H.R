##H bifurcation plot
library(tidyverse)

parms = c(H = 0,
           r = 1, K = 10,
           a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.06, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

###Create a parameter space with analytic way----

comp_out = expand.grid(H = seq(0, 1, by = 0.01))


comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = 6))) #Number = number of state variable
names(comp_out) <- c("H", "P1H", "P2H", "P1", "P2", "S1", "S2")

comp_out[, "S2"] = with((as.list(c(comp_out, parms))), {
  S2 = ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
  return(S2)
})

comp_out[, "S1"] = with((as.list(c(comp_out, parms))), {
  S1 = ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
  return(S1)
})


for(i in 1:dim(comp_out)[1]){
  temp_parms = parms
  temp_parms["H"] = comp_out[i, "H"]
  
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

comp_out %>%
  select(c(S1, S2, H)) %>%
  pivot_longer(names_to = "Si", values_to = "Abundance", -c(H)) %>%

ggplot(comp_out, mapping = aes(x = H, y = Abundance, color = Si)) +
  geom_line(lwd = 1)+
  labs(title = expression(m[1] < m[2]), x = "H", y = expression(S[i]^{"*"}), color = "Si")+
  scale_y_continuous() +
  scale_x_continuous() +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S1" = "Strain 1", "S2" = "Strain 2", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S1" = "darkblue", "S2" = "darkgreen", "H" = "red"))
  
  #theme(axis.title.y.right = element_text(angle = 90))

ggsave("Si and H 00536.png", width = 20, height = 11.25, units = "cm", dpi = 1600)

####dPi_dt and S* (Fix H and expanding S)-------

parms = list(H = 0.722230302, S = 0,
             r = 1, K = 10,
             a1 = 0.35, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
             b1 = 0.2, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
             o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

comp_out = expand.grid(S = seq(0, 1, by = 0.01))

          

###Create a parameter space with analytic way----

comp_out = expand.grid(S = seq(0, 5, by = 0.01))


comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                                       ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = 4))) #Number = number of state variable
names(comp_out) <- c("S", "dP1_dt", "dP2_dt", "dP1H_dt", "dP2H_dt")


for(i in 1:dim(comp_out)[1]){
  temp_parms = parms
  temp_parms["S"] = comp_out[i, "S"]
  
  
  
  comp_out[i, "dP1_dt"] = 
    with(parms, {
      (e1 * a1 * P1) * S - (b1 * P1) * H + (1 - DL) * (e1H * psi1 * a1 * P1H * S) - m1 * P1
    })
  
  comp_out[i, "dP2_dt"] = 
    with(parms, {
      (e2 * a2 * P2) * S - (b2 * P2) * H + (1 - DL) * (e2H * psi2 * a2 * P2H * S) - m2 * P2
    })
  
  comp_out[i, "dP1H_dt"] =
    with(parms, {
      (b1 * P1 * H) + DL * (e1H * psi1 * a1 * P1H * S) - (o1 + m1) * P1H
    })
    
  comp_out[i, "dP2H_dt"] =
    with(parms, {
      (b2 * P2 * H) + DL * (e2H * psi2 * a2 * P2H * S) - (o2 + m2) * P2H
    })
  
}


comp_out

comp_out %>%
  #select(c(S, H)) %>%
  pivot_longer(names_to = "dPi_dt", values_to = "rate", -c(H)) %>%
  
  ggplot(comp_out, mapping = aes(x = S, y = rate, color = dPi_dt)) +
  geom_line(lwd = 1)+
  labs(title = expression(m[1] < m[2]), x = "S", y = expression(dP_{i}/dt), color = "dPi_dt")+
  scale_y_continuous() +
  scale_x_continuous() +
  scale_colour_manual(labels = 
                        c("P1" = expression(P[1]), "P1H" = expression(P[1/H]),
                          "P2" = expression(P[2]), "P2H" = expression(P[2/H]),
                          "S1" = "Strain 1", "S2" = "Strain 2", "H" = "Hyper"),
                      values = c("P1" = "#BCAAA4", "P1H" = "#82491E",
                                 "P2" = "#B0BEC5", "P2H" = "#546E7A", 
                                 "S1" = "darkblue", "S2" = "darkgreen", "H" = "red"))


