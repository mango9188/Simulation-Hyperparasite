##H bifurcation plot
library(tidyverse)

parms = c(H = 0.0239533,
           r = 1, K = 10,
           a1 = 0.5, a2 = 0.5, psi1 = 1, psi2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.75, b2 = 0.45, m1 = 0.05, m2 = 0.055, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)

###Create a parameter space with analytic way----

comp_out = expand.grid(H = seq(0, 3, by = 0.01))


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

####dPi_dt and S*-------

S1 =
  with(parms, {
    ((b1 * H + m1) * (m1 + o1))/ (e1 * a1 * (m1 + o1) + e1H * psi1 * a1 * b1 * H)
  })

S2 = 
  with(parms, {
    ((b2 * H + m2) * (m2 + o2))/ (e2 * a2 * (m2 + o2) + e2H * psi2 * a2 * b2 * H)
  })

P1 = 
  with(parms, {
    r * (1 - (S1/K)) * (m1 + o1) / (a1 * (m1 + o1 + psi1 * b1 * H))
  })

P1H = 
  with(parms, {
    P1 * (b1 * H) / (m1 + o1)
  })


P2 = 
  with(parms, {
    r * (1 - (S2/K)) * (m2 + o2) / (a2 * (m2 + o2 + psi2 * b2 * H))
  })

P2H =
  with(parms, {
    P2 * (b2 * H) / (m2 + o2)
  })