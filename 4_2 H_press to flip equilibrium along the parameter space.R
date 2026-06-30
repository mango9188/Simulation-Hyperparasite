## 4_2 H_press to flip equilibrium along the parameter space ----
### Read the data to further subset different initial states and parameters
parms <- list(r = 1, K = 10,
              a1 = 0.38, a2 = 0.51, psi1 = 0.8, psi2 = 0.8, e1 = 0.5, e2 = 0.5,
              b1 = 0.2, b2 = 0.42, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
              o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.02, DL = 0)

P.space = 
  readRDS("Pre4A1_d002_psi08_rev") %>%
  #readRDS("Pre4A1_forPSpace") %>%
  filter(Stable_E == "P2H,C")
P.space = P.space[,-c(9:10)]
P.space$H_press = 0
head(P.space)

### Model parameters
times <- c(0, 16000)
state <- c(H = 0, P1H = 0, P2H = 0, P1 = 0, P2 = 0, S = 0)
# parameter = expand_grid(P.space[, c("m1", "m2", "H", "P1H", "P2H", "P1", "P2", "S")])

#function that will output the competition outcome
f_outcome = function(state, parms){
  pop_size = ode(func = M2,
                 times = times,
                 y = state,
                 parms = parms,
                 maxsteps = 80000)
  Outcome =
    paste(ifelse(pop_size[nrow(pop_size), c("H", "P1H", "P2H", "P1", "P2")] > 1e-7, "T", "F"), collapse = "")
  return(Outcome)
}

#function of bisection method
Bisection_H = function(min_H, max_H, state, parms){
  ###The goal is to find a maximum value of H that make competition outcome == E2H
  ###Min
  Min_state = state
  Min_state["H"] = min_H
  min_E = f_outcome(state = Min_state, parms = parms)
  
  ###Max
  Max_state = state
  Max_state["H"] = max_H
  max_E = f_outcome(state = Max_state, parms = parms)
  
  if(min_E == max_E){
    print(paste0("The equilibrium of minimum H is equal to maximum H (", min_H, ")"))
    return(NA)
  }
  
  loop_count = 0
  while(max_H - min_H > 1e-9 && loop_count < 100){
    ###Mid
    mid_H = (min_H + max_H)/2
    Mid_state = state
    Mid_state["H"] = mid_H
    mid_E = f_outcome(state = Mid_state, parms = parms)
    
    ###determine the value of mid_H
    if(mid_E == min_E){
      min_H = mid_H
    }else if(mid_E == max_E){
      max_H = mid_H
    }else{
      return(NA)
    }
    loop_count = loop_count + 1
  }
  return((min_H + max_H)/2)
}

Start_time = Sys.time()
for (i in 1:dim(P.space)[1]) {
  #setting the initial condition
  temp_state = unlist(P.space[i,  c("H", "P1H", "P2H", "P1", "P2", "S")])
  
  #setting the parameters
  temp_parms = parms
  temp_parms[c("m1", "m2")] = P.space[i, c("m1", "m2")]
  
  #run with the bisection method
  P.space[i,]$H_press = Bisection_H(min_H = min(P.space[i, "H"]),
                                    max_H = 8,
                                    state = temp_state,
                                    parms = temp_parms)
  if(i %% 10 == 0) print(paste0("Total: ", nrow(P.space), ", Now: ", i))
}
Ending_time = Sys.time()
Ending_time - Start_time
#saveRDS(P.space, file = "")


### Plot the result----
P2 =
  readRDS("Result 4 (H_press to flip equilibrium)") %>%
  filter(round(m1, 3) == 0.06) %>%
  mutate(delta_H = H_press - H) %>%
  ggplot()+
  geom_line(mapping = aes(x = m2, y = delta_H), linewidth = 0.8, color = "darkgreen")+
  scale_x_continuous(breaks = c(seq(0.002, 0.020, by = 0.006))) +
  labs(x = expression("Intrinsic mortality rate of pathogen strain 2"~ (m[2])), y = expression(Delta~H))

P1 =
  readRDS("Result 4 (H_press to flip equilibrium)") %>%
  filter(round(m2, 3) == 0.02) %>%
  filter(round(m1, 3) > 0.059) %>%
  mutate(delta_H = H_press - H) %>%
  ggplot()+
  geom_line(mapping = aes(x = m1, y = delta_H), linewidth = 0.8, color = "darkgreen")+
  labs(x = expression("Intrinsic mortality rate of pathogen strain 1"~ (m[1])), y = expression(Delta~H))

P1 + labs(tag = "(A)") + P2 + labs(tag = "(B)") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) &
  theme(plot.tag = element_text(size = 13, face = "bold"))

ggsave("Result_4_2 H_press to flip equilibrium (LABEL).png", width = 27.5, height = 14.3, units = "cm", dpi = 800)

