library(pracma)

# find time points when local peaks occur
peaks <- findpeaks(pop_size[, "H"])[ ,2]
peaks
tail(peaks)
# get period as time between peaks 
periods <- peaks[length(peaks)] - peaks[length(peaks) - 2] #3個 peak 之間的距離
# long-term average of H
avg_H <- mean(pop_size[(length(times) - periods + 1):length(times), "H"])
avg_H


###for the cycle situation
# for-loop to find all the long-term average of each species
avg = c(rep(NA, 6))
names(avg) = c("H", "P1H", "P2H", "P1", "P2", "S")
j = 1
species = c("H", "P1H", "P2H", "P1", "P2", "S")

sd(pop_size[(nrow(pop_size)-round(length(times)*0.3)):nrow(pop_size),"S"])

if (sd(pop_size[(nrow(pop_size)-round(length(times)*0.3)):nrow(pop_size),"S"]) > 1e-8){
  print("There is a cycle.")
  for (i in species) {
  peaks = findpeaks(pop_size[, i])[ ,2]
  periods = peaks[length(peaks)] - peaks[length(peaks) - 2]
  avg[j] = mean(pop_size[(length(times) - periods + 1):length(times), i])
  j = j + 1
  }
  print(avg)

}else{
  print("There is no cycle.")
}

avg

###for the cycle situation
#Check if the population is cycling
for(i in 1:dim(comp_out)[1]){
  
  temp_parms <- parms
  temp_parms["a1"] <- comp_out[i, ]$a1
  temp_parms["b1"] <- comp_out[i, ]$b1
  
  temp_out <- ode(func = M2, 
                  times = times, 
                  y = state, 
                  parms = temp_parms)
  
  N_final <- nrow(na.omit(temp_out))
  comp_out[i, 3:8] <- temp_out[N_final, -1] ##-1: except time
  
}



  if (sd(pop_size[(nrow(pop_size)-1000000):nrow(pop_size),"S"]) > 1e-8){
    #print(paste0("There is a cycle when a1 = ", temp_parms["a1"], ", b1 = ", temp_parms["b1"]))
    k = 1
    for (j in names(comp_out)[3:8]) {
      peaks = findpeaks(pop_size[, j])[ ,2]
      periods = peaks[length(peaks)] - peaks[length(peaks) - 1]
      avg[k] = mean(pop_size[(length(times) - periods + 1):length(times), j])
      k = k + 1 
    }
    k = 0
    comp_out[i, 3:8] = avg
  }else{
    N_final <- nrow(na.omit(temp_out)) #
    comp_out[i, 3:8] = temp_out[N_final, -1] #抓取最後一列的資料，排除時間
  }

avg

N_final <- nrow(na.omit(temp_out))
comp_out[i, 3:8] <- temp_out[N_final, -1] ##-1: except time

for (j in 9:10) {
  print(paste0("J = ", j))
  for (i in 1:6) {
    if (i > 4){
      print("STOP")
      break
    }else{
      print(i)
    }
  }
}


#For the case without P1H-----------------------------------------------------------------------------
avg = c(rep(NA, 5))
names(avg) = c("H", "P2H", "P1", "P2", "S")
j = 1
species = c("H", "P2H", "P1", "P2", "S")
for (i in species) {
  peaks = findpeaks(pop_size[, i])[ ,2]
  periods = peaks[length(peaks)] - peaks[length(peaks) - 1]
  avg[j] = mean(pop_size[(length(times) - periods + 1):length(times), i])
  j = j + 1
}
avg


#-----------------------------------------------------------------------------
# find time points when local peaks occur
peaks <- findpeaks(pop_size[, "N"])[ ,2]
peaks 

# get period as time between peaks 
periods <- peaks[length(peaks)] - peaks[length(peaks) - 1]
# long-term average of N
avg_N <- mean(pop_size[(length(times) - periods + 1):length(times), "N"])
avg_N

# long-term average of P
avg_P <- mean(pop_size[(length(times) - periods + 1):length(times), "P"])
avg_P

# equilibrium of N and P
E_np