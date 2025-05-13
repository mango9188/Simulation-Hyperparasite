library(deSolve)
library(tidyverse)
library(paletteer)

#### ODE setup
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



times <- seq(0, 100000, by = 0.1)
state <- c(H = 1, P1H = 0, P2H = 0, P1 = 2, P2 = 2, S = 5)
parms <- c(r = 1, K = 10,
           a1 = 0.5, a2 = 0.5, R1 = 1, R2 = 1, e1 = 0.5, e2 = 0.5,
           b1 = 0.45, b2 = 0.45, m1 = 0.05, m2 = 0.05, e1H = 0.5, e2H = 0.5,
           o1 = 0.8, o2 = 0.8, h1 = 1, h2 = 1, c1 = 0.9, c2 = 0.9, d = 0.03, DL = 0)


#### Set up parameter space of interest
#o1_min = 0.01
#o1_max = 0.8
#a1_min = 0.1
#a1_max = 1.5
#b1_min = 0.01
#b1_max = 1


#n <- 100
#R1_vec = seq(R1_min, R1_max, length.out = n)
#R2_vec = seq(R2_min, R2_max, length.out = n)
comp_out = expand.grid(a1 = seq(0, 1, by = 0.1), o1 = seq(0, 1, by = 0.1))


#### Create saving space for simulation output
#combine different parameters and different state variables together
comp_out <- as.data.frame(cbind(comp_out,
                                matrix(0, 
                                       nrow = dim(comp_out)[1],
                          ###dim(data)[1] is the number of row of data; [2] is col
                                       ncol = length(state))))
names(comp_out) <- c("a1", "o1", "H", "P1H", "P2H", "P1", "P2", "S")


#### Simulate the ODE across the parameter space
start_time <- Sys.time()
for(i in 1:dim(comp_out)[1]){
  
  temp_parms <- parms
  temp_parms["a1"] <- comp_out[i, ]$a1
  temp_parms["o1"] <- comp_out[i, ]$o1
  
  temp_out <- ode(func = M2, 
                  times = times, 
                  y = state, 
                  parms = temp_parms)
  
  N_final <- nrow(na.omit(temp_out)) 
  comp_out[i, 3:8] <- temp_out[N_final, -1] ##-1: exclude time
  
}
end_time <- Sys.time()
end_time - start_time

###Data analysis
extinct_thres = 1e-7
comp_out$Outcome =
  ifelse(comp_out[, "P1H"] + comp_out[, "P1"] < extinct_thres, "P2 win", 
         ifelse(comp_out[, "P2H"] + comp_out[, "P2"] < extinct_thres, "P1 win", "Coexist"))

comp_out$Outcome2 =
  apply(comp_out[, c("H", "P1H", "P2H", "P1", "P2")], 1, function(row){
    paste(ifelse(row > extinct_thres, "T", "F"), collapse = "")
  })

####Plot the result
#filter(filter(comp_out, round(a1, 3) >= 0.5, round(b1, 3) >= 0.45)

ggplot(comp_out, aes(x = a1, y = o1, z = Outcome, fill = Outcome2)) +
  geom_tile() +
  labs(title = "δ = 0", x = expression(α[1]), y = expression(ω[1]))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1)

##Save the simulation result
saveRDS(comp_out, "")
comp_out = readRDS("Pre1Sim3N")
######Read the data-------------------------------------
M2outcome = readRDS("M2outcomeDL1")

extinct_thres = 1e-7

M2outcome$Outcome =
  ifelse(M2outcome[, "P1H"] + M2outcome[, "P1"] < extinct_thres, "P2 win", 
         ifelse(M2outcome[, "P2H"] + M2outcome[, "P2"] < extinct_thres, "P1 win", "Coexist"))

M2outcome$Outcome2 =
  apply(M2outcome[, c("P1H", "P1", "P2H", "P2")], 1, function(row){
    paste(ifelse(row > extinct_thres, "T", "F"), collapse = "")
  })


####Fix omega-------------------------------------
o1_value = 0.3
####Plot the result
ggplot(filter(M2outcome, round(o1, 3) == o1_value), aes(x = a1, y = b1, z = Outcome2, fill = Outcome2)) +
  geom_tile() +
  labs(title = paste0("δ = 1, ω1 = ", o1_value), x = "α1", y = "β1" )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1)
####Fix beta-----------------------------------------
b1_value = 0.05
####Plot the result
ggplot(filter(M2outcome, b1 == b1_value), aes(x = a1, y = o1, z = Outcome2, fill = Outcome2)) +
  geom_tile() +
  labs(title = paste0("δ = 1,  β1 = ", b1_value), x = "α1", y = "ω1" )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1)


#--------Check which combination get most TTTT---------
M2outcome$o1 = as.factor(M2outcome$o1)
M2outcome$o1 = as.numeric(M2outcome$o1)

ggplot(M2outcome)+
  geom_bar(mapping = aes(x = Outcome2, fill = o1))

ggplot(M2outcome)+
  geom_bar(mapping = aes(x = Outcome, fill = Outcome2))+
  scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))

#Show the ratio
ggplot(M2outcome)+
  geom_bar(mapping = aes(x = Outcome, fill = Outcome2), position = "fill")+
  scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))


#----------------DL=0------------------------------------------------------------
M2outcome = readRDS("M2outcomeDL0")

extinct_thres = 1e-7

M2outcome$Outcome =
  ifelse(M2outcome[, "P1H"] + M2outcome[, "P1"] < extinct_thres, "P2 win", 
         ifelse(M2outcome[, "P2H"] + M2outcome[, "P2"] < extinct_thres, "P1 win", "Coexist"))

M2outcome$Outcome2 =
  apply(M2outcome[, c("P1H", "P1", "P2H", "P2")], 1, function(row){
    paste(ifelse(row > extinct_thres, "T", "F"), collapse = "")
  })

o1_value = 0.8
#### Visualize the competition outcomes in the grids
ggplot(filter(M2outcome, o1 == o1_value), aes(x = a1, y = b1, z = Outcome2, fill = Outcome2)) +
  geom_raster() +
  labs(title = paste0("δ = 0, ω1 = ", o1_value), x = "α1", y = "β1" )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) 

ggplot(filter(M2outcome, o1 == o1_value), aes(x = a1, y = b1, z = Outcome, fill = Outcome)) +
  geom_tile() +
  labs(title = paste0("δ = 0, ω1 = ", o1_value), x = "α1", y = "β1" )+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = as.character(paletteer_d("ggsci::default_igv")))+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed(ratio = 1) 

ggplot(filter(M2outcome, o1 == o1_value), aes(x = a1, y = b1, z = Outcome, fill = Outcome)) +
  geom_raster() +
  labs(title = paste0("ω1 = ", o1_value))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1") +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) 

#--------Check which combination get most TTTT---------
M2outcome$o1 = as.factor(M2outcome$o1)
M2outcome$o1 = as.numeric(M2outcome$o1)

ggplot(M2outcome)+
  geom_bar(mapping = aes(x = Outcome2, fill = o1))

ggplot(M2outcome)+
  geom_bar(mapping = aes(x = Outcome, fill = Outcome2))+
  scale_fill_brewer(palette = "Set1")

#Show the ratio
ggplot(M2outcome)+
  geom_bar(mapping = aes(x = Outcome, fill = Outcome2), position = "fill")+
  scale_fill_brewer(palette = "Set1")
