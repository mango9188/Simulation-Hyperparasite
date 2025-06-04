#Stability analisis

####Monoculture equilibrium----

Ji = with(parms, {
  matrix(data = c(
    -c1*b1*P1*H-d , h1*o1 , -c1*b1*H , 0,
        b1*P1 ,    -(m1+o1) , b1*P1 ,  0,
          0   ,   e1H*psi1*a1*S , e1*a1*S-b1*H-m1 , e1*a1*P1 + e1H*psi1*a1*P1H,
          0   ,   -psi1*a1*S , -a1*P1 , r*(1-S/K) - a1*P1 - psi1*a1*P1H - S*r/K
          ), 
    byrow = T, nrow = 4, ncol = 4)

  })
